"""Public multi-pair summary provenance tests; not biological validation."""

from scripts.summarize_public_multipair_colabfold_panel import existing_rows


def test_existing_rows_retain_exact_linker_sequences():
    existing = {
        "records": [
            {
                "candidate_id": "ggggs_10",
                "linker_length": 10,
                "mean_plddt": 80.0,
                "domain_a_mean_plddt": 90.0,
                "linker_mean_plddt": 40.0,
                "domain_b_mean_plddt": 75.0,
                "ptm": 0.5,
                "mean_cross_domain_pae_angstrom": 20.0,
                "pdb": "rank_001_alphafold2_ptm_model_1_seed_000.pdb",
                "pdb_sha256": "0" * 64,
                "score_json": "scores_rank_001_alphafold2_ptm_model_1_seed_000.json",
                "score_json_sha256": "1" * 64,
            }
        ]
    }

    rows = existing_rows(existing, {"ggggs_10": "GGGGSGGGGS"})

    assert rows[0]["linker_sequence"] == "GGGGSGGGGS"
    assert rows[0]["construct_id"] == "1UBQ_2GB1__ggggs_10"
