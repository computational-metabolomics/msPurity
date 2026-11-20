"""Regenerate reference.json with Python and a checkout of biorack/blink."""
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

revision = "dfdf0c13672b9f5e5bbf59187ea4f6cd04627c64"
repo = Path(sys.argv[1]).resolve()
if subprocess.check_output(["git", "-C", repo, "rev-parse", "HEAD"], text=True).strip() != revision:
    raise ValueError("Check out the pinned BLINK revision " + revision)
sys.path.insert(0, str(repo))
from blink import discretize_spectra, score_sparse_spectra

query = [[50.0625, 50.1875, 50.20, 60.0, 90], [4, 1, 9, 16, 0]]
references = [
    [[49.75, 49.875, 50.0, 50.25, 50.375, 60.125], [1, 4, 1, 9, 1, 4]],
    [[50.0625, 50.1875, 50.20, 60.0], [4, 1, 9, 16]],
    [[75.0, 90.0], [2, 1]],
    [[49.75, 50.5], [1, 4]],
]
cases = []
for width, tolerance, power in [(0.001, 0.01, 0.5), (0.125, 0.375, 0.5),
                                (0.125, 0.3125, 1), (0.125, 0.125, 0)]:
    data = discretize_spectra([np.array(query, dtype=float)],
                              [np.array(r, dtype=float) for r in references],
                              [200], [200] * len(references), tolerance=tolerance,
                              bin_width=width, intensity_power=power)
    scores = score_sparse_spectra(data, gpu=False)
    cases.append(dict(binWidth=width, blurTolerance=tolerance, intensityPower=power,
                      cosine=scores["mzi"].toarray()[0].tolist(),
                      count=scores["mzc"].toarray()[0].tolist()))
out = dict(revision=revision, query=query, references=references, cases=cases)
Path(__file__).with_name("reference.json").write_text(json.dumps(out, indent=2) + "\n")
