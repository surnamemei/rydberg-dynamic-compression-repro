"""Write the follow-up job plan consumed by progress06.py (FU-A..FU-D, 11_followup_preregistration.json)."""
import json
from common import OUT

plan = []
for f in (0.8, 1.0, 1.2):
    plan.append({"tag": "fu_cw_steps", "job_id": f"tau_{f:.1f}xE1", "weight": 0.14, "family": "FU-A tau step", "realization": 0, "xi": f"step@{f:.1f}E1", "power_db": 0})
for p in (0, 3):
    for c in ("cw", "qpsk"):
        plan.append({"tag": "fu_cw_steps", "job_id": f"plateau_p{p:+d}_{c}", "weight": 0.21, "family": f"FU-B plateau {c}", "realization": 0, "xi": "60us", "power_db": p})
std = ["xi_0.05", "xi_0.1", "xi_0.25", "xi_0.5", "xi_1", "xi_2", "xi_4"]
longv = ["xi_0.05", "xi_1", "xi_4", "xi_8", "xi_16", "xi_32"]
for r in range(4):
    for p in (3, 0):
        for v in std:
            plan.append({"tag": "fu_cwcarrier", "job_id": f"dwell_r{r}_{v}_p{p:+d}_Nd4001_dt1", "weight": 1.0, "family": "FU-C dwell CW-carrier", "realization": r, "xi": v[3:], "power_db": p})
for tag, fam in (("fu_long_qpsk", "FU-D long QPSK-carrier"), ("fu_long_cw", "FU-D long CW-carrier")):
    for r in range(4):
        for v in longv:
            plan.append({"tag": tag, "job_id": f"dwell_r{r}_{v}_p+3_Nd4001_dt1", "weight": 2.9, "family": fam, "realization": r, "xi": v[3:], "power_db": 3})
(OUT / "plan_followup.json").write_text(json.dumps(plan, indent=1))
print(len(plan), "jobs; weighted main-job equivalents:", round(sum(j["weight"] for j in plan), 1))
