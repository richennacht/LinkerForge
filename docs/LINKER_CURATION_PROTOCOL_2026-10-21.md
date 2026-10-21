# LinkerForge mined-corpus curation protocol

Date: 21 October 2026 (Asia/Calcutta)
Scope: linker identity and source evidence; **not** fusion compatibility or biological success

## Decision

The 27,766-record `mined_and_natural_linkers_v1` composite remains a valid,
provenance-stratified discovery, retrieval and weak-supervision corpus. It is
not renamed or promoted to a fusion-compatibility dataset.

Its evidence classes remain separate:

| Evidence class | Records | Permitted use |
| --- | ---: | --- |
| Classifier-mined high-confidence calls | 6,556 | retrieval, weak supervision, review sampling |
| Classifier-mined provisional candidates | 21,039 | discovery and review sampling |
| CATH/PDB natural structural observations | 91 | natural-junction authority and geometry reference |
| DisProt curator-validated DFL annotations | 57 | independent DFL authority |
| Published engineered templates | 23 | cited template retrieval, not general compatibility proof |

## Review packet

`scripts/build_linker_curation_review_packet.py` verifies every source hash and
record count before selecting a deterministic, stratified review packet.
Configuration is frozen in `configs/linker-curation-review-v1.json`.

The initial packet contains 971 records:

- 360 high-confidence classifier calls;
- 440 provisional classifier candidates;
- all 91 CATH/PDB natural observations;
- all 57 DisProt curated DFL annotations;
- all 23 published templates.

Classifier-mined rows are stratified by source corpus, length, authority tier
and classifier score band. Homology-cluster IDs are preserved for later
family-level splitting.

## Independent review

The packet creates separate `curator_a.jsonl` and `curator_b.jsonl` files.
Curators must work independently until both files are locked. For each record
they assess:

1. coordinate correctness;
2. whether the span is an inter-domain connection;
3. direct, indirect, absent or uncertain linker-behaviour evidence;
4. flexible/disordered, structured, mixed, unresolved or uncertain state;
5. whether the cited source supports the record;
6. include, exclude or uncertain disposition.

Disagreements require a separately versioned adjudication record. Reviewer
identity, review time and source identifiers must be retained. Agreement and
stratum-specific precision are reported only after the two reviews are locked.

## Claim firewall

Every review record fixes `compatibility_outcome` to `not_assessed`. A reviewed
natural linker is not thereby evidence that the same sequence will express,
fold or function between an arbitrary engineered domain pair.

The compatibility ranker therefore remains blocked until a separate
literature-derived engineered-fusion outcome dataset exists with exact domain
sequences, linker sequences, construct context, measured outcomes, citations
and family-disjoint splits. New wet-lab work is outside the requested scope;
published measurements may be curated without claiming they are newly
generated evidence.

## Training gate

The active DFL detector firewall is independent of this curation packet and
remains fail-closed. Neither CAID3's spent labels nor DP94's sealed labels may
be used for training, filtering, threshold selection or feature design. No new
detector or compatibility-ranker training is authorized merely because this
review packet exists.
