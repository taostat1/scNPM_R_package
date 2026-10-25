"""Record the upstream AJGM stopping state without changing its iterations."""
import json
import math
import sys


def checked_ajgm(fn, path, **kwargs):
    record = {}

    def profile(frame, event, arg):
        if event == "return" and frame.f_code is fn.__code__:
            loc = frame.f_locals
            gap = float(loc["diff_1"])
            threshold = float(loc["EMthreshold"])
            record.update(iterations=int(loc["z"]), objective_change=gap,
                          threshold=threshold, max_iterations=int(loc["MAX_iter"]),
                          converged=math.isfinite(gap) and gap <= threshold)

    previous = sys.getprofile()
    try:
        sys.setprofile(profile)
        result = fn(**kwargs)
    finally:
        sys.setprofile(previous)
    if not record:
        raise RuntimeError("AJGM stopping diagnostics unavailable; check upstream version")
    with open(path, "w", encoding="utf-8") as stream:
        json.dump(record, stream, indent=2)
    return result
