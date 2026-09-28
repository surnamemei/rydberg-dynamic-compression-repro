"""Add the upstream RydbergComms simulator (commit 4095724) to this directory without overwriting any file.

The simulator is not part of this repository (it is the prior work of J. Zhu and L. Dai, MIT License). This script
clones https://github.com/johnzja/RydbergComms, exports commit 4095724 with `git archive`, and extracts every file
that does not already exist here. Files of this repository always take precedence.
"""
import io
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

URL = "https://github.com/johnzja/RydbergComms.git"
COMMIT = "4095724"
HERE = Path(__file__).resolve().parent

with tempfile.TemporaryDirectory() as tmp:
    subprocess.run(["git", "clone", "--quiet", URL, tmp], check=True)
    full = subprocess.run(["git", "-C", tmp, "rev-parse", "--verify", COMMIT + "^{commit}"],
                          check=True, capture_output=True, text=True).stdout.strip()
    data = subprocess.run(["git", "-C", tmp, "archive", "--format=tar", full], check=True, capture_output=True).stdout
added = skipped = 0
with tarfile.open(fileobj=io.BytesIO(data)) as tar:
    for m in tar.getmembers():
        dst = HERE / m.name
        if m.isdir():
            dst.mkdir(parents=True, exist_ok=True)
            continue
        if dst.exists() or dst.with_name(dst.name.lower()).exists() or dst.with_name(dst.name.upper()).exists():
            skipped += 1
            continue
        dst.parent.mkdir(parents=True, exist_ok=True)
        tar.extract(m, HERE, set_attrs=False, filter="data") if sys.version_info >= (3, 12) else tar.extract(m, HERE)
        added += 1
print(f"RydbergComms {full[:12]}: added {added} upstream files, kept {skipped} existing files of this repository")
