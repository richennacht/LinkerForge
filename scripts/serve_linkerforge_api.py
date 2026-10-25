"""Small JSON API for the LinkerForge web interface.

Endpoints
---------
GET  /api/status   dataset tier counts (from the hash-pinned manifest) and result-file provenance
GET  /api/demo     the precomputed sfGFP (2B3P) - linker - mCherry (2H5Q) demonstration
POST /api/score    score one user-supplied linker sequence at the demo junction
                   body: {"sequence": "...", "rigidity": "flexible"|"rigid", "cleavable": bool}

Nothing is computed from invented inputs: the junction geometry, port accessibilities and flank
context come from `reports/part4_end_to_end_ranker_and_sfgfp_mcherry_2026-10-19.json`, which is
produced by `scripts/evaluate_part4_end_to_end_pipeline.py`.  The service never returns a
fusion-success probability.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from http.server import BaseHTTPRequestHandler, HTTPServer
from pathlib import Path
import sys
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from linkerforge.pipeline import PairwiseCompatibilityRanker  # noqa: E402

RESULT = ROOT / "reports" / "part4_end_to_end_ranker_and_sfgfp_mcherry_2026-10-19.json"
MANIFEST = ROOT / "data" / "processed" / "mined_and_natural_linkers_v1" / "manifest.json"
DISCLAIMER = (
    "Deterministic geometric / sequence rules and a relative score on one schematic junction. "
    "Not a prediction of expression, folding or activity."
)
ALLOWED = set("ACDEFGHIKLMNPQRSTVWY")


class Service:
    def __init__(self) -> None:
        self.result = json.loads(RESULT.read_text(encoding="utf-8"))
        self.result_sha = hashlib.sha256(RESULT.read_bytes()).hexdigest()
        self.ranker = PairwiseCompatibilityRanker(
            feature_names=tuple(self.result["ranker"]["feature_names"]),
            weights=tuple(self.result["ranker"]["weights"]),
        )
        self.manifest = json.loads(MANIFEST.read_text(encoding="utf-8")) if MANIFEST.exists() else {}

    def status(self) -> dict[str, Any]:
        return {
            "dataset_id": self.result.get("dataset_id"),
            "dataset_manifest_sha256": self.result.get("dataset_manifest_sha256"),
            "manifest_tier_counts": self.manifest.get("tier_counts") or self.manifest.get("record_counts"),
            "result_file": RESULT.name,
            "result_file_sha256": self.result_sha,
            "disclaimer": DISCLAIMER,
        }

    def score(self, body: dict[str, Any]) -> tuple[int, dict[str, Any]]:
        seq = str(body.get("sequence", "")).upper().strip()
        if len(seq) > 80 or (seq and not set(seq) <= ALLOWED):
            return 400, {"error": "sequence must be <= 80 residues of the 20 standard amino acids"}
        rigidity = body.get("rigidity", "flexible")
        if rigidity not in ("flexible", "rigid"):
            return 400, {"error": "rigidity must be 'flexible' or 'rigid'"}
        ctx = self.result["fusion_demo"]["context"]
        target = int(body.get("target_length") or max(len(seq), 5))
        score, card, feats = self.ranker.score_candidate(
            seq,
            target_length=target,
            required_path_length_A=ctx["path_A"],
            requested_rigidity=rigidity,
            requested_cleavable=bool(body.get("cleavable", False)),
            left_flank=ctx["left_flank"],
            right_flank=ctx["right_flank"],
            port_a_sasa=ctx["rsa_a"],
            port_b_sasa=ctx["rsa_b"],
        )
        return 200, {
            "sequence": seq,
            "length": len(seq),
            "junction_path_A": ctx["path_A"],
            "relative_score": round(score, 4),
            "scorecard": card.to_dict(),
            "features": {k: round(v, 4) for k, v in feats.items()},
            "disclaimer": DISCLAIMER,
        }


def make_handler(svc: Service):
    class Handler(BaseHTTPRequestHandler):
        def _send(self, code: int, payload: dict[str, Any]) -> None:
            body = json.dumps(payload).encode("utf-8")
            self.send_response(code)
            self.send_header("Content-Type", "application/json; charset=utf-8")
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Access-Control-Allow-Origin", "http://127.0.0.1:5173")
            self.send_header("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
            self.send_header("Access-Control-Allow-Headers", "Content-Type")
            self.end_headers()
            self.wfile.write(body)

        def do_OPTIONS(self) -> None:  # noqa: N802
            self._send(200, {"ok": True})

        def do_GET(self) -> None:  # noqa: N802
            if self.path.startswith("/api/status"):
                self._send(200, svc.status())
            elif self.path.startswith("/api/demo"):
                self._send(200, svc.result["fusion_demo"])
            else:
                self._send(404, {"error": "not_found"})

        def do_POST(self) -> None:  # noqa: N802
            if not self.path.startswith("/api/score"):
                self._send(404, {"error": "not_found"})
                return
            try:
                length = int(self.headers.get("Content-Length", "0"))
                body = json.loads(self.rfile.read(length).decode("utf-8") or "{}")
            except (ValueError, json.JSONDecodeError):
                self._send(400, {"error": "invalid JSON body"})
                return
            code, payload = svc.score(body)
            self._send(code, payload)

        def log_message(self, format: str, *args: Any) -> None:  # noqa: A002
            return

    return Handler


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=8765)
    args = ap.parse_args()
    svc = Service()
    print(f"LinkerForge API on http://{args.host}:{args.port} (result sha256 {svc.result_sha[:12]}...)")
    HTTPServer((args.host, args.port), make_handler(svc)).serve_forever()


if __name__ == "__main__":
    main()
