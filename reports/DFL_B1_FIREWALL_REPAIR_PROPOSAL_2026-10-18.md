# B1 firewall repair: sequence-only null diagnostic and resolution path

Status: **proposal, not a change to the binding preflight**. B1 remains open and
no training job is authorized. This report uses only training-candidate
sequences and DP94's sealed sequence-only registry. It does not open DP94 or
CAID3 labels, prediction scores, or per-protein errors.

## Finding

The current hard-exclusion rule accepts a BLOSUM62 local alignment with at
least 50 aligned residues, 25% identity and an alignment length at least half
the shorter protein. Its implementation divides `hit.length` by the shorter
protein length; that quantity can exceed 1 when gaps extend the alignment, so
it is not a bounded sequence-coverage fraction. Moreover, local alignments
between compositionally similar disordered proteins may satisfy the rule
without being the same protein family.

`scripts/audit_dfl_firewall_null.py` tests this without touching benchmark
labels. From the pinned 205 DisProt DFL-positive-candidate records, it samples
40 candidates and 40 DP94 sequences per seed. Each observed pair is compared
with a whole-residue shuffled candidate and a 12-residue-block shuffled
candidate. Shuffles preserve composition, and block shuffles additionally
retain short local runs. The declared screen gave:

| Seed | Original hits / 1,600 | Residue-shuffle hits / 1,600 | Block-shuffle hits / 1,600 |
| --- | ---: | ---: | ---: |
| 173 | 30 | 23 | 23 |
| 271 | 31 | 18 | 9 |
| 811 | 57 | 37 | 27 |
| **Pooled descriptive count** | **118 / 4,800** | **78 / 4,800** | **59 / 4,800** |

The pairs overlap across seeds, so the pooled numbers are not independent
confidence intervals. The shuffles are not known nonhomologous biological
families. Nevertheless, a large fraction of screen-positive events survives
destruction of long-range sequence order. This is concrete evidence that the
25% local rule alone has poor specificity as a *family* definition. It does
not prove which individual candidate–DP94 connections are false, nor does it
permit us to ignore an existing exclusion.

Source JSONL SHA-256 `2bc6726545f47256bcec70376bd51d689db0b5c80a64f45a221b9d88bb0c4c8d`;
candidate manifest SHA-256 `f4fad95c86f969b4659345b7c0d38afc0dbe66d5f8ba9d3a3bcdcc07f6b04ea0`;
DP94 sequence-only manifest SHA-256 `fe64e7c688cfac48ed0dac73bcf09daae5cb0adfd0e03e9d530824d3d92c0051`.
Audit executable SHA-256 `9965dc1ec7fb29f9400502b4f905b9e6be73855bd76f28bbb7cd197c309ac3af`.

Reproduce from a checkout that has the ignored pinned source files:

```powershell
C:\Users\Perthro\Documents\ChatGPT\Folddiffusion\.venv-gpu\Scripts\python.exe scripts\audit_dfl_firewall_null.py --source-root C:\Users\Perthro\Documents\ChatGPT\Folddiffusion --sample-size 40 --seed 173
C:\Users\Perthro\Documents\ChatGPT\Folddiffusion\.venv-gpu\Scripts\python.exe scripts\audit_dfl_firewall_null.py --source-root C:\Users\Perthro\Documents\ChatGPT\Folddiffusion --sample-size 40 --seed 271
C:\Users\Perthro\Documents\ChatGPT\Folddiffusion\.venv-gpu\Scripts\python.exe scripts\audit_dfl_firewall_null.py --source-root C:\Users\Perthro\Documents\ChatGPT\Folddiffusion --sample-size 40 --seed 811
```

## Proposed resolution, in order

1. **Keep the benchmark seal and exact exclusions.** DP94's authors document a
   DFL annotation track, so it cannot be discarded merely to increase the
   training count. Keep exact accession, exact sequence, containment and
   strong full-length homolog exclusions hard. CAID3 remains a spent
   retrospective registry, not a tuning set.
2. **Pre-register a biologically grounded family rule using sequences only.**
   Test a composition-aware full-sequence clustering method (for example,
   MMseqs2 with low-complexity masking and reciprocal coverage) and an
   independent Pfam/InterPro domain-family audit. Choose and document tool
   versions, thresholds and a shuffle-null false-connection rate *before*
   opening any held-out labels. Investigate borderline pairs manually from
   sequence/domain evidence; unresolved pairs stay excluded. Do not replace
   the present rule in the launcher until this evidence and a contract review
   exist.
3. **Rebuild source-specific masks.** Use the curated DisProt flexible-linker
   regions and published DFL training files as positive/unknown supervision
   only where source semantics support it. Unannotated residues remain `x`,
   not `0`; structural DLD and the 884 UniProt weak candidates stay separate
   tiers. Reserve a family-disjoint development set and an independently
   adjudicated no-linker set for specificity.
4. **Re-run the capacity audit under both policies.** Report retained proteins,
   spans and family clusters under the old rule and the proposed validated
   rule, plus every reclassified edge and its reason. If the new rule still
   leaves too little labelled data, stop: obtain more curated DFL annotations
   or commission prospective adjudication. Never fill the gap with inferred
   domain boundaries.
5. **Only then create B1.** Version the split and all benchmark sequence
   registries, update the binding contract and both preflights together after
   review, and run them before any GPU allocation. DP94 benchmark fitness and
   competitor training exposure remain separate unresolved questions; a
   publishable superiority claim still needs a genuinely untouched common-
   target test with independently justified labels and coverage.

This is the shortest defensible path through the current *methodological*
deadlock. It may recover data lost to non-specific screening, but it makes no
promise of enough DFL labels or superior model performance.
