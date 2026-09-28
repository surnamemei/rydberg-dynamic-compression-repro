# Deviations from 00_plan.md

Each entry: time stamp, the state when written (what had been run or seen), the deviation and its reason.

## D1 — Part II provisional flag (implementation fix, no rule change)

- **Written:** 2026-09-28T03:18:46+10:00
- **State:** Part I complete and classified; the Part II coarse scan (31 IFs x 3 LO fields) analysed once; the refinement rule fired for A_LO 0.35 and 0.5 V/m (8 IFs); no refinement run existed yet.
- **Deviation:** on its first call, `tvf_p2_analyze.py` decided the refinement *after* checking for missing runs, so it labelled its coarse-scan output `final: true`. The class from that call (SUGGESTIVE ONLY) was treated as provisional, as the plan requires. The script now recomputes the missing-run check after deciding the refinement. Decision rules, metrics and the refinement decision are unchanged; the class is final only after `tvf_run.py IIr` and a second call.

## D2 — Part III selection-step crash (implementation fix, no rule change)

- **Written:** 2026-09-28T03:20:42+10:00
- **State:** Part III training Gram matrix, the 30 fits and the validation pass had completed. The script then stopped with a KeyError before writing any output, so no test result had been computed or seen. The only information visible was the key of the validation-selected grid point, printed in the error message.
- **Deviation:** the grid stored the depths 1, 2 and 4 µs as integers, and the selected key was rebuilt with floats. The grid now uses float depths, and the selected key is looked up by value. Grid, split, selection rule and metrics are unchanged. The script was rerun from the start.
