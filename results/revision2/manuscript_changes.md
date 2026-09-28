# Manuscript changes: final receiver-physics validation pass (2026-09-28)

Decision rules: `00_decision_rules.json` (committed 6f73843, before any run), `01_rule_clarifications.json` (c535691, before any run),
`02_deviations.json` (D1–D2 2940fae, before any 15 MHz result; D3 52d5c51, before any weak-probe run; D4–D5 da38640, after all runs).
Every number added to the text is a row of `manuscript/numbers_ledger.csv` (IDs `R2_*`, built by `revision2()` in
`manuscript/build_numbers_ledger.py`). Files changed: `manuscript/main_draft.md`, `manuscript/supplement_draft.md`,
`submission/tqe/main.tex`, `submission/tqe/supplement.tex`, figure scripts and PDFs.

## Required by the results

| Result | Change | Where |
|---|---|---|
| P1: IF dependence classified DOES_NOT_SURVIVE (phase effect absent at 10 and 15 MHz; at 10 MHz the static surrogate reproduces most of the dwell gain change) | New subsection "Dependence on IF and probe strength" (VII-F / 7.6); phase claims scoped to the 5 MHz IF in the abstract, the contribution list, VII-E "What this shows", Discussion (P1dB incompleteness, descriptor, regimes, behavioural modelling, design), Limitations (configuration bullet rewritten; new "matched levels across configurations" bullet; ordering-result bullet) and Conclusion | main |
| P1: dwell gain contrast at other IFs | "Other intermediate frequencies" paragraph in IV-B; "Not established" list in IV-D | main |
| P1: harmonic-zone leakage | Numerical-lessons sentence (VIII-F); leakage paragraph in S10 (leakage cannot account for X, because the constant-amplitude envelope carries no baseline modulation) | main, supplement |
| P2: terminology rule | X renamed "excess gain reduction" (the "extra compression" name failed the rule: X vs X_mag differ by up to 0.11 at OP1); X_coh may be called compression (phase share at most 16%) | main (metrics, VII-A, VII-C, VII-D, Fig. 8 caption and axis label), supplement (S6, Table S8 caption) |
| P2: decomposition | New "Magnitude or phase?" paragraph in VII-A; Table S10 and the terminology paragraph in S10 | main, supplement |
| P3: weaker probe | Result in VII-F; slow-component sentence in V-A; S10 weak-probe paragraph; Table S9 column | main, supplement |
| P4: internal trace | "Operational signature" subsection renamed "Operational Signature (Post Hoc) and Internal State"; the sentence "Internal atomic coherences were not recorded" replaced by the prespecified internal-state paragraph (wording "associated with"; no mechanism inferred); Limitations "Operational, not microscopic, mechanism" bullet updated; Discussion VIII-C updated | main, supplement S10 |
| Configuration changes | "The only change is the LO field" replaced by a statement that only configuration parameters change (LO field; IF and probe Rabi frequency in VII-F) | main II-A |
| New supplementary material | Section S10, Tables S9 and S10, Fig. S7; Acknowledgment now lists Supplementary Figs. S1–S7 | supplement, main |

## Requested restructuring (P5)

- Section IV now reads: A. Controlled test signals (unchanged); B. Dwell ordering primarily changes the fitted gain (gain contrast first; states explicitly that the prespecified D-based dwell robustness test failed); C. Declustering robustly changes residual distortion (D_ref first, FIR-span robustness, D second); D. What this establishes, and what it does not.
- Figures swapped so that the gain decomposition leads: the former Supplementary Fig. S5 (fitted gain and D_ref vs dwell) is now Fig. 5, with a new panel (d) showing the declustering D_ref change vs FIR span (+3 and +6 dB); the former Fig. 5 (self-normalized D vs dwell, single-slow-state check) is now Supplementary Fig. S5, with the D-based details moved to S9 ("The self-normalized dwell result"). Cross-references in Sections V and S9 updated. Fig. 4(d) caption now names the metric (self-normalized D) and points to Fig. 5(d).
- The old self-normalized D dwell result no longer leads Section IV.

## Requested wording audit (P6)

Every visible use of "preregistered"/"preregistration" was replaced (main.tex: 34 lines; supplement.tex: 15 lines; figure artwork: 7 labels in Figs. 2, 7, 8 and S1):
- "Evidence labels" paragraph rewritten: decision rules fixed in time-stamped project records before the runs ("prespecified"); not deposited with an external registry before the runs; the public archive makes them auditable after the fact.
- Abstract: "with most decision rules fixed before the runs"; "prespecified distortion criterion".
- Code Availability (both variants): "time-stamped decision-rule records".
- All other occurrences: "prespecified".
- Remaining occurrences are archived file names inside hidden SRC comments of the Markdown drafts only.
- The figure audit now fails on "preregist" or "extra compression" in any figure text.
- Reproducibility package: the `preregistration/` folder is now `decision_records/` (nine records, including the three revision2 records), and README, PROVENANCE and REPRODUCE use the same wording.

## Checks after the changes

- `check_draft.py`: 0 unresolved ledger IDs, 0 unmatched numerical literals, 0 figure/table reference errors; abstract 246 words.
- Markdown vs LaTeX: all decimal values agree (differences are section numbers only).
- `make_all.py`: 15 figures, 215 checks, ALL PASS.
- `main.pdf` 17 pages, `supplement.pdf` 9 pages; no LaTeX errors or undefined references; no new overfull boxes; remaining placeholders: repository URL (once) and DOI (twice).

# Approved rescope (2026-09-28)

The author approved the RESCOPE classification. No new simulations or analyses were run; the science is frozen. Every number in the rewritten text is a row of the numbers ledger (new rows: `R2_gcoh_rel_OP1`, `R2_Xcoh_IF`, `R2_criteria`).

| Item | Change |
|---|---|
| Title | "Configuration-Dependent Dynamic Compression in a Simulated Rydberg Atomic Receiver" (alternates prepared: "Dynamic Compression and Temporal Memory in a Simulated Rydberg Atomic Receiver"; "Configuration Sensitivity of Large-Signal Dynamics in a Simulated Rydberg Atomic Receiver") |
| Abstract | Rewritten from scratch (202 words): reference configuration; surrogate failure; dwell → gain; declustering → D_ref; failed prespecified dwell-distortion robustness test stated; 2.5–3.1 µs memory; coherent gain 0.47 of CW; phase effect absent at 10/15 MHz and with a weaker probe at matched levels; memory persists; configuration dependent; no mechanism |
| Introduction | Configuration question added to the problem statement; contributions replaced by C1–C5 |
| Methods | "Reference configuration" defined (5 MHz IF, standard probe; OP1–OP3 differ in the LO field) |
| Section III | Retitled "CW-Derived Surrogates Mispredict Waveform Degradation in the Reference Configuration" |
| Section IV | IV-D renamed "Interpretation" (different observables; no single temporal-order law) |
| Section V | τ_atom renamed τ_ref (archived name kept for traceability); "the D-based dwell dependence is specific to the phase carrier ... follows from the power-driven memory" withdrawn and replaced by "occurs with and without a phase carrier; consistent with a power-driven memory" |
| Section VI | Quasi-static limit explicitly restricted to unmodulated excitation |
| Section VII | Retitled "Configuration-Specific Phase-Modulated Gain Reduction"; A reference configuration, B amplitude vs phase, C operating-point dependence at 5 MHz, D IF dependence, E weaker-probe check, F internal-coherence observation |
| Section VIII | Rewritten: evidence-summary Table 5; A CW characterization does not transfer; B different observables; C configuration dependence; D what remains robust; E what does not generalize; F implications |
| Section IX | Limitations rewritten; literature-currency bullet removed |
| Section X | Conclusion rewritten; ends with the call to identify the physical origin of the configuration dependence and test it experimentally |
| Figures | New Fig. 9 (configuration dependence, from existing P1/P3 results); Fig. S7 reduced to the P4 internal traces; captions use "reference configuration"; OP comparisons labelled "5 MHz IF" (Figs. 8, S6); τ label τ_ref in Figs. 1, 4–7, S3, S5 |
| Final skeptical review | 0 blocking, 3 important, 5 minor comments; all fixed in text. Important: matched-level qualifier added wherever the absence of the phase effect is stated (abstract, C4, VII, VIII-E, conclusion, Fig. 8 and 9 captions); declustering scoped to the reference configuration in VIII-D; C2 scoped. Minor: "bit-identically (dwell gains to 1e-14)"; weak-probe and configuration-dependence wording; post-hoc gain contrast's later prespecified confirmation noted. Residual minor points that are not text defects: 18-page length; float placement; X_coh not computed for the weak probe (stated in Table S9); internal Markdown table numbering differs from the LaTeX numbering (drafts only) |
