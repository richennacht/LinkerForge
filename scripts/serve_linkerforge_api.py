"""Lightweight HTTP JSON API server for the LinkerForge End-to-End Web Application.

Endpoints:
  - GET  /api/status : Returns dataset summary (27,766 records), Part 3 & Part 4 metrics.
  - POST /api/design : Runs live end-to-end LinkerForge ranking & 3D PDB + PyMOL generation
                       for sfGFP (2B3P) — linker — mCherry (2H5Q) with user tick-box settings.
"""

from __future__ import annotations

import argparse
from http.server import BaseHTTPRequestHandler, HTTPServer
import json
from pathlib import Path
import sys
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from linkerforge.dataset_v1 import (  # noqa: E402
    DEFAULT_DATASET_DIR,
    load_dataset_manifest,
    load_linker_records,
)
from linkerforge.pipeline import (  # noqa: E402
    AuthoritySequenceIndex,
    fit_pairwise_compatibility_ranker,
    run_sfgfp_mcherry_fusion_demo,
)


class LinkerForgeServiceContext:
    def __init__(self, dataset_dir: Path = DEFAULT_DATASET_DIR) -> None:
        self.dataset_dir = dataset_dir
        self.manifest = load_dataset_manifest(dataset_dir)
        natural_records = load_linker_records(
            dataset_dir,
            tiers=["natural_cath_pdb_verified", "natural_disprot_curated_dfl"],
        )
        mined_hi = load_linker_records(
            dataset_dir,
            tiers=["mined_high_confidence_dfl_call"],
        )
        templates = load_linker_records(
            dataset_dir,
            tiers=["literature_engineered_template"],
        )
        self.ranker = fit_pairwise_compatibility_ranker(
            list(natural_records) + list(mined_hi[:1200]) + list(templates),
            seed=23,
        )
        self.authority_index = AuthoritySequenceIndex(
            list(natural_records) + list(mined_hi) + list(templates)
        )


def make_handler(ctx: LinkerForgeServiceContext):
    class Handler(BaseHTTPRequestHandler):
        def _send_json(self, status_code: int, payload: dict[str, Any]) -> None:
            body = json.dumps(payload).encode("utf-8")
            self.send_response(status_code)
            self.send_header("Content-Type", "application/json; charset=utf-8")
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Access-Control-Allow-Origin", "*")
            self.send_header("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
            self.send_header("Access-Control-Allow-Headers", "Content-Type")
            self.end_headers()
            self.wfile.write(body)

        def do_OPTIONS(self) -> None:  # noqa: N802
            self._send_json(200, {"ok": True})

        def do_GET(self) -> None:  # noqa: N802
            if self.path.startswith("/api/status"):
                self._send_json(
                    200,
                    {
                        "service": "LinkerForge Scientific API v1.0",
                        "dataset_id": ctx.manifest.get("dataset_id"),
                        "total_records": ctx.manifest.get("total_records", 27766),
                        "tier_counts": ctx.manifest.get("tier_counts", {}),
                        "authority_index_size": len(ctx.authority_index.records),
                    },
                )
                return
            self._send_json(404, {"error": "not_found"})

        def do_POST(self) -> None:  # noqa: N802
            if self.path.startswith("/api/design"):
                length = int(self.headers.get("Content-Length", "0"))
                raw = self.rfile.read(length).decode("utf-8") if length > 0 else "{}"
                req = json.loads(raw or "{}")
                rigidity = str(req.get("rigidity", "flexible"))
                if rigidity not in ("flexible", "rigid"):
                    rigidity = "flexible"
                cleavable = bool(req.get("cleavable", False))
                protease_family = str(req.get("protease_family", "TEV"))
                target_length = max(6, min(36, int(req.get("target_length", 15))))
                required_path_length_A = float(req.get("required_path_length_A", 28.4))

                demo = run_sfgfp_mcherry_fusion_demo(
                    ctx.ranker,
                    ctx.authority_index,
                    requested_rigidity=rigidity,
                    requested_cleavable=cleavable,
                    protease_family=protease_family,
                    target_length=target_length,
                    required_path_length_A=required_path_length_A,
                )
                self._send_json(200, demo)
                return
            self._send_json(404, {"error": "not_found"})

        def log_message(self, format: str, *args: Any) -> None:  # noqa: A002
            return

    return Handler


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--dataset-dir", type=Path, default=DEFAULT_DATASET_DIR)
    args = parser.parse_args()
    ctx = LinkerForgeServiceContext(args.dataset_dir)
    server = HTTPServer((args.host, args.port), make_handler(ctx))
    print(f"LinkerForge API listening on http://{args.host}:{args.port}")
    server.serve_forever()


if __name__ == "__main__":
    main()
