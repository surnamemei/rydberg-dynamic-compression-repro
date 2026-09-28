# Manuscript Changes Required by the Revision Results

The preregistration (`00_revision_preregistration.json`, committed before any run) requires every change to be recorded here. Each change was made in both the Markdown record (`manuscript/main_draft.md`, `manuscript/supplement_draft.md`, with SRC comments pointing to `RV_*` ledger rows) and the submission LaTeX (`submission/tqe/main.tex`, `supplement.tex`).

## Changes forced by a preregistered outcome

### P1: the dwell contrast's surrogate excess does not survive D_ref (preregistered verdict: DOES_NOT_SURVIVE)

**Abstract.** The sentence "Reordering identical field-amplitude samples changes in-band distortion by up to 0.39 (95% CI ±0.02), which the surrogates do not reproduce" now reads "…lowers the large-signal gain from 0.470 to 0.360 of its small-signal value, which the surrogates do not reproduce". The gain decomposition is post hoc.

**Section 1, finding 2.** "In-band nonlinear distortion depends…" becomes "The large-signal gain and in-band distortion depend…".

**Section 2.3.** Added a FIR-span robustness sentence and the definition of D_ref (the metric was already defined post hoc in Section 6).

**Section 4.2.** Added the "What the change in D consists of" paragraph:
- FIR span does not matter (the D change stays at 0.379–0.389);
- under D_ref the LTI+static surrogate's residual rises more (excess −0.007 [−0.009, −0.005]), so the criterion failed;
- the unmatched part is the fitted gain (−0.110 [−0.118, −0.102]);
- declustering passes every variant.

The section title became "Temporal order changes gain and distortion", and the "What this shows" paragraph was rewritten.

**Fig. 5 caption.** States that D combines gain and residual changes, and points to Supplementary Fig. S5.

**Section 5.3: time scale.** The growth in D between ξ ≈ 0.25 and 4 is attributed to the residual component. The gain falls mainly by ξ = 0.25 (post hoc).

**Section 5.3: single-slow-state surrogate.** New sentence: it reproduces the long-dwell gain drop approximately, but not its onset (post hoc).

**Section 5.3: carrier ablation.** The statement "The strong dwell dependence of Section 4 is therefore tied to the phase carrier…; power-driven memory alone … not that dependence" was withdrawn. It held only for the self-normalized D. In the fitted gain, the dwell dependence is larger with the unmodulated carrier (−0.38 [−0.42, −0.35]; post hoc).

**Section 5 "What this shows".** The ordering sensitivity now "arises at sub-microsecond to microsecond dwell"; added that the gain changes at shorter dwell than a single τ_atom memory predicts.

**Section 8.4.** The single-state surrogate "captures clustering and, approximately, the long-dwell gain drop, but not its sub-microsecond onset or the phase effect".

**Section 8.6.** Added the numerical lesson that self-normalized D can register gain compression as distortion.

**Limitations.** New bullet "Metric-dependent ordering result"; "Distortion, not information rate" became "Gain and distortion, not information rate".

**Conclusion, second finding.** Rewritten around gain compression and declustering-dependent distortion.

**Supplement S9.** Table S7, Fig. S5 (new).

### P2: six sequences per phase case

**OP1 reference case.** In Section 7.1 the value X = 0.520 is now reported per sequence, with the six-sequence mean 0.494 [0.431, 0.557] and SD 0.06. The earlier two-sequence agreement (0.520/0.523) understated the sequence-to-sequence spread.

**Abstract.** "0.43 to 0.21" becomes "0.43 to 0.22 (six sequences)", and "0.16 vs 0.52" becomes "0.16 vs 0.49".

**Section 7.2, shape at OP1.** Steps add +0.06 [+0.01, +0.11]; the 900 ns change is unresolved. The claim "insensitivity to transition shape does not hold at OP2" now reads "transition shape matters much more at OP2".

**Section 7.3.** Six-sequence rate values were added. The OP1 turnover is resolved (−0.07 [−0.13, −0.01]). The deterministic sequence lies above all six random sequences.

**Section 7.5.**
- OP2 six-sequence values: X = 0.165 [0.134, 0.196], and OP1 exceeds OP2 in every case.
- 900 ns suppression is stronger at OP2 than at OP1 (difference +0.10 [+0.03, +0.16]).
- "rate dependence saturates" became "does not turn over". FAST−REF is +0.05 at OP2, against −0.07 at OP1 (difference −0.12 [−0.17, −0.07]).

**Fig. 8.** Panels (b) and (e) now show the six-sequence means with 95% CIs; the caption was updated.

**Limitations.** The Numerics bullet now gives the sequence SDs.

**Supplement.** Table S8 (new).

### P3: amplitude sweep

**Onset.** X becomes material near P1dB: from 0 dB at OP1 and from −3 dB at OP2.
- Added to Section 7.1 and the abstract ("emerging near P1dB").
- Also in Supplement S9 and Fig. S6(a).

### P4: offset sweep and a third operating point (classification: MONOTONIC_IN_LO)

**Section 7.2.** Finer sweep: R = 0.51 at −0.25 MHz and 1.17 at +0.125 MHz at OP1; at OP2, R stays within 0.98–1.10 for |δf| ≤ 0.75 MHz.

**Section 7.5.** New "A third operating point" paragraph. OP1 is not isolated: X and the selectivity S5 grow monotonically with the LO field. The "What this shows" paragraph now covers three LO fields.

**Limitations.** "Two LO operating points only" became "Three LO operating points, the third with a minimal contrast only". The OP1 selectivity statement was narrowed to the near-IF feature.

**Conclusion.** The phase effect "emerged near P1dB and persisted at two further LO fields, with magnitude growing with the LO field".

**Supplement.** Fig. S6(b), (c).

### P5: frequency-dependent static surrogate (verdict: FAILS)

**Where added.** A new paragraph in Section 7.2; one clause each in the abstract, Section 8.3, Section 8.4 and the Conclusion.

**Result.** The surrogate predicts X of about 0 (median). A post-hoc mean-based check gives the opposite sign.

### P6: numerical closure

**Section 2.6.** Double-precision reruns gave X = 0.5199 vs 0.5200 and changed the long-dwell D_ref gap by 0.007 percentage points. Archived jobs rerun bit-identically.

**Limitations.** The precision statement was updated. The 200 µs persistence criterion failed (X of 0.45–0.57 across windows, no decay), and this is now stated.

**Supplement S9.** Numerical-closure list, including wall times.

## Editorial

**Supplement.** New Section S9 (Tables S7, S8; Figs. S5, S6).

**Ledger and audits.**
- `manuscript/numbers_ledger.csv` gained 109 `RV_*` rows.
- `check_draft.py`: 0 unmatched literals.
- `make_all.py`: ALL PASS for 14 figures (201/201 checks).
