"""GPU benchmark of SAM+AF(31) and SAM+ESM2+AF(31) on the 363-protein v2 full-sequence split.

Completes the unfinished v4_train_sam_af (`both` 31-dim AlphaFold + 3D CA geometry arm)
and adds a regularized k=1 bottleneck ESM-2 35M (`facebook/esm2_t12_35M_UR50D`, 480-d)
branch across 3 seeds (7, 17, 31) on the firewall-clean 254-train / 109-dev DisProt v2 split.
Never reads or evaluates CAID3 or DP94 benchmark labels.
"""

from __future__ import annotations

import json
import math
import os
from pathlib import Path
import random
import sys
import tarfile
import time
import numpy as np
import torch
import torch.nn as nn
from transformers import AutoModel, AutoTokenizer

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SCRATCH = Path.home() / ".gemini" / "antigravity" / "brain" / "fe7cd636-644a-45ec-851c-0bf737a38f82" / "scratch" / "dfl_exp"
SCRATCH = Path(os.environ.get("DFL_EXP_SCRATCH", str(DEFAULT_SCRATCH)))
DEFAULT_TAR = Path.home() / "Documents" / "DFL-Scratch-Project" / "dfl_exp_code_and_results.tar.gz"
TAR_PATH = Path(os.environ.get("DFL_EXP_TAR", str(DEFAULT_TAR)))
OUT_JSON = ROOT / "reports" / "dfl_v2_esm2_af_classifier_benchmark_2026-10-26.json"


def ensure_scratch_extracted() -> None:
    SCRATCH.mkdir(parents=True, exist_ok=True)
    if not (SCRATCH / "v4_struct_features.py").exists():
        with tarfile.open(TAR_PATH, "r:gz") as tf:
            tf.extractall(SCRATCH, filter="fully_trusted")
    scratch_posix = str(SCRATCH).replace("\\", "/")
    for fpath in SCRATCH.glob("*.py"):
        txt = fpath.read_text(encoding="utf-8")
        if "/home/claude/dfl_exp" in txt or "/tmp/claude-0" in txt:
            txt = txt.replace("/home/claude/dfl_exp", scratch_posix)
            txt = txt.replace(
                "/tmp/claude-0/-home-claude/2e3c2d6e-3ccd-554d-bc19-532e1e14a294/scratchpad/sam_like",
                scratch_posix + "/scratch_sam",
            )
            fpath.write_text(txt, encoding="utf-8")


D, NQ, HEADS, FF, DROP = 64, 16, 4, 128, 0.1
SPAN_DIM, PRIOR_MIN, PRIOR_MAX = 8, 4, 96


def sinus_gpu(length: int, dim: int, dev: torch.device) -> torch.Tensor:
    pos = torch.arange(length, dtype=torch.float32, device=dev).unsqueeze(1)
    div = torch.exp(torch.arange(0, dim, 2, dtype=torch.float32, device=dev) * (-math.log(10000.0) / dim))
    pe = torch.zeros(length, dim, device=dev)
    pe[:, 0::2] = torch.sin(pos * div)
    pe[:, 1::2] = torch.cos(pos * div)
    return pe


class SAMAFGPU(nn.Module):
    def __init__(self, fdim: int, esm_dim: int = 0) -> None:
        super().__init__()
        self.esm_dim = esm_dim
        self.conv = nn.Conv1d(21, D, 5, padding=2)
        self.fb = nn.Sequential(nn.Linear(fdim, 32), nn.GELU(), nn.Linear(32, D))
        if esm_dim > 0:
            self.esm_branch = nn.Sequential(
                nn.LayerNorm(esm_dim),
                nn.Dropout(0.5),
                nn.Linear(esm_dim, 16),
                nn.GELU(),
                nn.Dropout(0.25),
                nn.Linear(16, D),
            )
            self.esm_scale = nn.Parameter(torch.tensor(0.25))
        self.enc = nn.TransformerEncoderLayer(D, HEADS, FF, DROP, batch_first=True, norm_first=True)
        self.q = nn.Parameter(torch.randn(NQ, D) * 0.02)
        x = np.log2(np.geomspace(PRIOR_MIN, PRIOR_MAX, NQ))[:, None]
        freqs = 2.0 ** np.arange(SPAN_DIM // 2) * 0.5
        feat = np.concatenate([np.sin(x * freqs), np.cos(x * freqs)], 1)
        self.register_buffer("span_feat", torch.tensor(feat, dtype=torch.float32))
        self.span = nn.Linear(SPAN_DIM, D)
        self.sa = nn.MultiheadAttention(D, HEADS, batch_first=True)
        self.ca = nn.MultiheadAttention(D, HEADS, batch_first=True)
        self.mlp = nn.Sequential(nn.Linear(D, D), nn.GELU(), nn.Linear(D, D))
        self.n1, self.n2, self.n3 = nn.LayerNorm(D), nn.LayerNorm(D), nn.LayerNorm(D)
        self.wq = nn.Linear(D, D, bias=False)
        self.wr = nn.Linear(D, D, bias=False)
        self.iou = nn.Linear(D, 1)
        self.bias = nn.Parameter(torch.zeros(1))

    def forward(
        self,
        x_oh: torch.Tensor,
        pad: torch.Tensor,
        f: torch.Tensor,
        x_esm: torch.Tensor | None = None,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        h = self.conv(x_oh.transpose(1, 2)).transpose(1, 2) + sinus_gpu(x_oh.shape[1], D, x_oh.device)
        h = h + self.fb(f) * f[..., -1:]
        if self.esm_dim > 0 and x_esm is not None:
            h = h + self.esm_scale * self.esm_branch(x_esm)
        r = self.enc(h, src_key_padding_mask=pad)
        bsz = x_oh.shape[0]
        q = (self.q + self.span(self.span_feat)).unsqueeze(0).expand(bsz, -1, -1)
        q = self.n1(q + self.sa(q, q, q, need_weights=False)[0])
        q = self.n2(q + self.ca(q, r, r, key_padding_mask=pad, need_weights=False)[0])
        q = self.n3(q + self.mlp(q))
        ml = torch.einsum("bkd,bld->bkl", self.wq(q), self.wr(r)) / math.sqrt(D) + self.bias
        iou = torch.sigmoid(self.iou(q)).squeeze(-1)
        sig = torch.sigmoid(ml)
        p = (iou.unsqueeze(-1) * sig).max(dim=1).values
        return p, sig


def sam_losses(
    p: torch.Tensor,
    sig: torch.Tensor,
    y: torch.Tensor,
    m: torch.Tensor,
    pw: torch.Tensor,
) -> tuple[torch.Tensor, torch.Tensor]:
    pc = p.clamp(1e-6, 1.0 - 1e-6)
    bce = -(pw * y * torch.log(pc) + (1.0 - y) * torch.log(1.0 - pc))
    l_bce = (bce * m).sum() / m.sum().clamp(min=1.0)
    ym = (y * m).unsqueeze(1)
    mm = m.unsqueeze(1)
    inter = (sig * ym).sum(-1)
    den = (sig * mm).sum(-1) + ym.sum(-1)
    dice = 1.0 - (2.0 * inter + 1.0) / (den + 1.0)
    ok = (m.sum(-1) > 0).float()
    l_dice = (dice.min(-1).values * ok).sum() / ok.sum().clamp(min=1.0)
    return l_bce, l_dice


def main() -> None:
    ensure_scratch_extracted()
    sys.path.insert(0, str(SCRATCH))
    import common  # type: ignore
    import v2_common as vc  # type: ignore
    import v4_struct_features as VF  # type: ignore
    from bilstm import make_batches  # type: ignore

    vc.use_policy("v2")
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}", flush=True)

    esm_cache = SCRATCH / "esm2_35m_features.npz"
    all_rows = [
        json.loads(line)
        for line in (SCRATCH / "v2_labels.jsonl").read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if esm_cache.exists():
        z = np.load(esm_cache)
        esm_feats = {k: z[k] for k in z.files}
        print(f"Loaded cached ESM-2 features: {len(esm_feats)}", flush=True)
    else:
        t0 = time.time()
        tok = AutoTokenizer.from_pretrained("facebook/esm2_t12_35M_UR50D", local_files_only=True)
        esm = AutoModel.from_pretrained("facebook/esm2_t12_35M_UR50D", local_files_only=True).to(device).eval()
        esm_feats = {}
        with torch.no_grad():
            for r in all_rows:
                acc, seq = r["accession"], r["sequence"]
                chunks = []
                for s in range(0, len(seq), 1022):
                    sub = seq[s : s + 1022]
                    inp = {k: v.to(device) for k, v in tok([sub], return_tensors="pt", add_special_tokens=True).items()}
                    with torch.amp.autocast("cuda", enabled=(device.type == "cuda")):
                        h = esm(**inp).last_hidden_state[0, 1 : 1 + len(sub)].float().cpu().numpy()
                    chunks.append(h)
                esm_feats[acc] = np.concatenate(chunks, axis=0).astype(np.float16)
        np.savez_compressed(esm_cache, **esm_feats)
        print(f"Cached ESM-2 features in {time.time() - t0:.1f}s", flush=True)

    stats = VF.get_stats()
    fdim = VF.dim("both")
    print(f"Built v4 struct features (both): fdim={fdim}", flush=True)

    def load_split_items(sp: str) -> list:
        out = []
        for a, s, l in vc.load("v2", sp):
            oh = common.onehot(s)
            y, m = common.targets(l)
            f_arr = VF.feats(a, s, "both", stats)
            e_arr = esm_feats[a].astype(np.float32)
            out.append((a, oh, y, m, f_arr, e_arr))
        return out

    def collate_gpu(its: list, idxs: list[int], use_esm: bool):
        max_len = max(its[i][1].shape[0] for i in idxs)
        bsz = len(idxs)
        fd = its[0][4].shape[1]
        x = torch.zeros(bsz, max_len, 21, device=device)
        y = torch.zeros(bsz, max_len, device=device)
        m = torch.zeros(bsz, max_len, device=device)
        f = torch.zeros(bsz, max_len, fd, device=device)
        xe = torch.zeros(bsz, max_len, 480, device=device) if use_esm else None
        lens = torch.zeros(bsz, dtype=torch.long, device=device)
        for b, i in enumerate(idxs):
            _, xo, yo, mo, fo, eo = its[i]
            n = xo.shape[0]
            x[b, :n] = torch.from_numpy(xo).to(device)
            y[b, :n] = torch.from_numpy(yo).to(device)
            m[b, :n] = torch.from_numpy(mo).to(device)
            f[b, :n] = torch.from_numpy(fo).to(device)
            if xe is not None:
                xe[b, :n] = torch.from_numpy(eo).to(device)
            lens[b] = n
        pad = torch.arange(max_len, device=device).unsqueeze(0) >= lens.unsqueeze(1)
        return x, y, m, f, xe, lens, pad

    pw, _, _ = vc.pos_weight("v2", "train")
    pwt = torch.tensor(float(pw), dtype=torch.float32, device=device)

    tr = load_split_items("train")
    dv = load_split_items("dev")

    results_dev: dict[str, dict[str, np.ndarray]] = {}
    summary_records: dict[str, dict] = {
        "v2_sam_af4_onehot_recorded_baseline": {
            "seed_ensemble_3": {"AP": 0.2770, "AUROC": 0.8390, "Fmax": 0.3670},
            "source": "DFL-Scratch-Project/HANDOFF.md",
        }
    }

    for name, use_esm in [("v4_sam_af31_onehot", False), ("v7_sam_af31_esm2_35m_bottleneck", True)]:
        seed_preds = []
        per_seed_metrics = {}
        for seed in [7, 17, 31]:
            t0 = time.time()
            random.seed(seed)
            np.random.seed(seed)
            torch.manual_seed(seed)
            model = SAMAFGPU(fdim, esm_dim=480 if use_esm else 0).to(device)
            if use_esm:
                esm_params = list(model.esm_branch.parameters()) + [model.esm_scale]
                esm_ids = {id(p) for p in esm_params}
                base_params = [p for p in model.parameters() if id(p) not in esm_ids]
                opt = torch.optim.AdamW(
                    [
                        {"params": base_params, "lr": 1e-3, "weight_decay": 0.0},
                        {"params": esm_params, "lr": 5e-4, "weight_decay": 2e-2},
                    ]
                )
            else:
                opt = torch.optim.Adam(model.parameters(), lr=1e-3)
            rng = random.Random(seed)
            best_ap = -1.0
            best_sc = None
            best_mt = None
            best_ep = 0
            for ep in range(1, 31):
                model.train()
                for idxs in make_batches(tr, rng, 60.0):
                    x, y, m, f, xe, lens, pad = collate_gpu(tr, idxs, use_esm)
                    if m.sum() == 0:
                        continue
                    p, sig = model(x, pad, f, xe)
                    lb, ld = sam_losses(p, sig, y, m, pwt)
                    loss = lb + ld
                    opt.zero_grad()
                    loss.backward()
                    opt.step()
                model.eval()
                sc = {}
                with torch.no_grad():
                    for idxs in make_batches(dv, random.Random(0), 0.0):
                        x, y, m, f, xe, lens, pad = collate_gpu(dv, idxs, use_esm)
                        p, _ = model(x, pad, f, xe)
                        p_np = p.cpu().numpy()
                        for b, i in enumerate(idxs):
                            sc[dv[i][0]] = p_np[b, : int(lens[b])].copy()
                yd, sd = vc.flatten("v2", "dev", sc)
                mt = vc.full_eval(yd, sd)
                if mt["AP"] > best_ap:
                    best_ap = mt["AP"]
                    best_sc = sc
                    best_mt = mt
                    best_ep = ep
            assert best_sc is not None and best_mt is not None
            seed_preds.append(best_sc)
            per_seed_metrics[str(seed)] = {
                "best_epoch": best_ep,
                "AP": round(float(best_mt["AP"]), 4),
                "AUROC": round(float(best_mt["AUROC"]), 4),
                "Fmax": round(float(best_mt["Fmax"]), 4),
                "seconds": round(time.time() - t0, 1),
            }
            print(
                f"{name} seed={seed} best_ep={best_ep} ({time.time() - t0:.1f}s): "
                f"AP={best_mt['AP']:.4f}, AUROC={best_mt['AUROC']:.4f}, Fmax={best_mt['Fmax']:.4f}",
                flush=True,
            )

        ens_sc = {a: np.mean([sp[a] for sp in seed_preds], axis=0) for a in seed_preds[0]}
        yd, sd = vc.flatten("v2", "dev", ens_sc)
        ens_mt = vc.full_eval(yd, sd)
        results_dev[name] = ens_sc
        summary_records[name] = {
            "per_seed": per_seed_metrics,
            "seed_ensemble_3": {
                "AP": round(float(ens_mt["AP"]), 4),
                "AUROC": round(float(ens_mt["AUROC"]), 4),
                "Fmax": round(float(ens_mt["Fmax"]), 4),
            },
        }
        print(
            f"==> {name} 3-SEED ENSEMBLE: "
            f"AP={ens_mt['AP']:.4f}, AUROC={ens_mt['AUROC']:.4f}, Fmax={ens_mt['Fmax']:.4f}",
            flush=True,
        )

    comb_sc = {
        a: 0.5 * results_dev["v4_sam_af31_onehot"][a]
        + 0.5 * results_dev["v7_sam_af31_esm2_35m_bottleneck"][a]
        for a in results_dev["v4_sam_af31_onehot"]
    }
    yd, sd_comb = vc.flatten("v2", "dev", comb_sc)
    comb_mt = vc.full_eval(yd, sd_comb)
    summary_records["combined_6model_sam_af31_plus_esm2_ensemble"] = {
        "AP": round(float(comb_mt["AP"]), 4),
        "AUROC": round(float(comb_mt["AUROC"]), 4),
        "Fmax": round(float(comb_mt["Fmax"]), 4),
    }
    print(
        f"==> COMBINED 6-MODEL ENSEMBLE (v4_sam_af31 + v7_sam_af31_esm2_bottleneck): "
        f"AP={comb_mt['AP']:.4f}, AUROC={comb_mt['AUROC']:.4f}, Fmax={comb_mt['Fmax']:.4f}",
        flush=True,
    )

    OUT_JSON.write_text(
        json.dumps(
            {
                "date": "2026-10-26",
                "split": "v2_full_sequence_disprot_254train_109dev",
                "dev_proteins": 109,
                "dev_residues": int(len(yd)),
                "dev_positives": int(yd.sum()),
                "dev_positive_rate": round(float(yd.mean()), 6),
                "benchmark_seal": "CAID3 and DP94 never read or evaluated",
                "models": summary_records,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
