import csv
from pathlib import Path
import p1db_comm_stage0 as m
ROOT=Path(__file__).resolve().parent
rows=[];stats=[]
for rep in range(4):
 st,rr=m.run("OFDM",20261004+100*rep,[-15,-10,-6,-3,0,3,6])
 rows.extend(rr);stats.append(st)
with open(ROOT/"p1db_comm_ofdm_fair.csv","w",newline="") as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
with open(ROOT/"waveform_statistics_ofdm_fair.csv","w",newline="") as f:
 w=csv.DictWriter(f,fieldnames=list(stats[0]));w.writeheader();w.writerows(stats)

