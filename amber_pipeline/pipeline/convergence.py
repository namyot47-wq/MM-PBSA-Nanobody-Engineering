#pipeline/convergence.py
import subprocess
import numpy as np
import re

def parse_sander_out(out_file: str, key: str) -> list[float]:
    values = []
    pattern = re.compile(rf"{key}\s*=\s*(-?\d+\.\d+)")
    with open(out_file) as f:
        for line in f:
            m = pattern.search(line)
            if m:
                values.append(float(m.group(1)))
    return values
#Statistical test to determine tail convergence through standard deviation of tail / mean value of tail
def tail_stable(values: list[float], tol_pct: float = 1.0, tail_frac: float = 0.2) -> bool:
    arr = np.array(values)
    tail = arr[int(len(arr) * (1 - tail_frac)):]
    pct_stdev = 100 * tail.std() / abs(tail.mean())
    return pct_stdev < tol_pct