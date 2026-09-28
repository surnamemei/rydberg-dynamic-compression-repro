# Decision-rule records: provenance

The study's decision rules were fixed in time-stamped project records before the corresponding runs. They are internal records: they were **not** deposited with an external registry before the runs. This archive makes them available so that the rules and their timing can be audited after the fact. The earlier files keep their historical names (`*preregistration*.json`, `00_preregistered_criteria.json`); the article calls the tests they govern *prespecified*.

The 13 records remain at their original paths under `results/`, next to the runs they govern; this directory holds byte-identical copies (same SHA-256). Each file carries the time at which it was written, as a field inside the file.

| File | In-file timestamp | First committed to the private project history | Committed separately, before its runs? | SHA-256 |
|---|---|---|---|---|
| `results/stage06_dwell_physics/00_preregistered_criteria.json` (copy: `stage06_dwell_physics__00_preregistered_criteria.json`) | `written_before_results`: 2026-09-26T12:06:36+10:00 | 2026-09-26T16:40:25+10:00 | no (same commit as its results) | `0c5eec2c35c1bc72…` |
| `results/stage06_dwell_physics/11_followup_preregistration.json` (copy: `stage06_dwell_physics__11_followup_preregistration.json`) | `written_before_results`: 2026-09-26T16:46:05+10:00 | 2026-09-26T19:09:10+10:00 | no (same commit as its results) | `f5fc95a3d08972fa…` |
| `results/stage06_dwell_physics/phase_mechanism/00_preregistration.json` (copy: `stage06_dwell_physics__phase_mechanism__00_preregistration.json`) | `written_before_any_run`: 2026-09-26T19:44:20+10:00 | 2026-09-26T19:55:50+10:00 | no (same commit as its results) | `f9a8068dab8a7dee…` |
| `results/final_validation/regen_stage05/00_protocol.json` (copy: `final_validation__regen_stage05__00_protocol.json`) | `written_before_any_run`: 2026-09-26T20:35:28+10:00 | 2026-09-26T20:35:28+10:00 | yes | `1fbfd35db94a9de1…` |
| `results/final_validation/replication/00_preregistration.json` (copy: `final_validation__replication__00_preregistration.json`) | `written_before_any_run`: 2026-09-26T20:35:28+10:00 | 2026-09-26T20:35:28+10:00 | yes | `be5e96f4006476b9…` |
| `results/revision/00_revision_preregistration.json` (copy: `revision__00_revision_preregistration.json`) | `written_before_any_run`: 2026-09-27T14:54:00+10:00 | 2026-09-27T14:54:05+10:00 | yes | `4ce45766e0a0eb98…` |
| `results/revision2/00_decision_rules.json` (copy: `revision2__00_decision_rules.json`) | `written_before_any_run`: 2026-09-27T23:53:30+10:00 | 2026-09-27T23:54:10+10:00 | yes | `5a762c0fe7fb8bf6…` |
| `results/revision2/01_rule_clarifications.json` (copy: `revision2__01_rule_clarifications.json`) | `written_before_any_run`: 2026-09-28T00:02:36+10:00 | 2026-09-28T00:03:04+10:00 | yes | `d3bbe581931fa807…` |
| `results/revision2/02_deviations.json` (copy: `revision2__02_deviations.json`) | `written`: 2026-09-28T00:15:34+10:00 | 2026-09-28T00:15:48+10:00 | separately, during execution; appended in 3 commits (each entry states which results existed when it was written) | `11827a7de3e8c17f…` |
| `results/tqe_viability/00_decision_rules.json` (copy: `tqe_viability__00_decision_rules.json`) | `written_before_any_run`: 2026-09-28T01:55:49+10:00 | 2026-09-28T01:56:34+10:00 | yes | `2ab479b01b3af851…` |
| `results/tqe_viability/01_deviations.json` (copy: `tqe_viability__01_deviations.json`) | `written (first entry)`: 2026-09-28T02:06:27+10:00 | 2026-09-28T02:28:02+10:00 | separately, during execution; appended in 1 commits (each entry states which results existed when it was written) | `75730ac299bd066a…` |
| `results/tqe_viability_final/00_plan.md` (copy: `tqe_viability_final__00_plan.md`) | `Written`: 2026-09-28, before any run of this campaign (the in-file time of the commit that adds this file is the registration time; internal time-stamped project record, not an external registration). | 2026-09-28T02:37:07+10:00 | yes | `2a7e98aae6142872…` |
| `results/tqe_viability_final/00_deviations.md` (copy: `tqe_viability_final__00_deviations.md`) | `Written`: 2026-09-28T03:18:46+10:00 | 2026-09-28T02:37:07+10:00 | separately, during execution; appended in 2 commits (each entry states which results existed when it was written) | `d398e8f2ca5ef2d7…` |

**What this record can and cannot show.**
- The in-file timestamps are self-reported.
- Files marked 'yes' were committed on their own, before any of their runs.
- The deviation record of the final checks (`results/revision2/02_deviations.json`) was written during execution, by design; each entry states which results existed when it was written.
- Files marked 'no' were first committed together with their results. For those, the project history does not establish the order independently.
- The private history is not part of this repository.

**Test codes used inside the files, and the names used in the article:**
- replication `H_gen`: magnitude criterion;
- `H_shape`: transition-shape criterion;
- `H_detune`: detuning criterion.
