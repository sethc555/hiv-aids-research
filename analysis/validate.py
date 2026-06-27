#!/usr/bin/env python3
"""Standing VALIDATOR for the within-host TIP model assumptions. Reads assumptions.json and prints:
 (1) DASHBOARD        -- counts by status / role / controversy.
 (2) DIG-HERE QUEUE   -- load-bearing/structural assumptions NOT 'supported' (prioritized by kills*).
 (3) BLIND-SPOT SURFACER -- scans each model's parameter dict / module constants and flags parameters
     with NO catalogued assumption: a written-but-uncatalogued modelling choice = a candidate blind spot.
 (4) HONEST UNKNOWNS  -- the explicitly-acknowledged gaps + the conceptual-frame limit the tool can't see.
 (5) --run -- also execute the LIGHTWEIGHT verify_*.py (verify_analytic / verify_loo) under the 4 GB cap.
     (verify_claims.py is the heavy stochastic aggregate -- run it separately, capped.)

Pattern borrowed from the sibling T1D repo. Run after every model change. Exit 0 iff no blind spots and
no unparsed sources.   python3 validate.py [--run]
"""
import argparse
import glob
import json
import os
import re
import subprocess
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
CATALOG = os.path.join(HERE, "assumptions.json")
# parameter SOURCES: (file, kind, dict-name)
SOURCES = [("tip_model.py", "dict", "P"),
           ("tip_model_p13_wm.py", "dict", "QS"),
           ("tip_model_p4_reservoir.py", "consts", None)]
HEAVY = {"verify_claims.py"}     # stochastic aggregate -- never auto-run here (OOM-safety)


def _read(fn):
    p = os.path.join(HERE, fn)
    return open(p).read() if os.path.exists(p) else None


def source_params(fn, kind, name):
    """Return (params, parsed_ok). parsed_ok False = file missing or dict not found (FLAG; no silent skip)."""
    src = _read(fn)
    if src is None:
        return [], False
    if kind == "dict":
        m = re.search(rf'\b{name}\s*=\s*dict\((.*?)\)', src, re.S)
        if not m:
            return [], False
        block = re.sub(r'#.*', '', m.group(1))
        return sorted(set(re.findall(r'(\w+)\s*=', block))), True
    return sorted(set(re.findall(r'^([A-Z][A-Z0-9_]*)\s*=\s*[-+\d.]', src, re.M))), True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="store_true", help="also run the lightweight verify_*.py (4 GB cap)")
    args = ap.parse_args()
    cat = json.load(open(CATALOG))
    A = cat["assumptions"]

    print("=" * 78)
    print(f"ASSUMPTION VALIDATOR — {len(A)} catalogued assumptions (HIV within-host TIP model)")
    print("=" * 78)
    print("\n(1) DASHBOARD")
    print("    status     :", dict(Counter(a["status"] for a in A)))
    print("    role       :", dict(Counter(a["role"] for a in A)))
    print("    controversy:", dict(Counter(a["controversy"] for a in A)))

    print("\n(2) DIG-HERE QUEUE  (load-bearing/structural AND not 'supported';  * = kills the result if wrong)")
    queue = [a for a in A if a["role"] in ("load-bearing", "structural") and a["status"] != "supported"]
    queue.sort(key=lambda a: (not a["kills_result_if_wrong"], a["status"]))
    if not queue:
        print("    (none — every load-bearing/structural assumption is 'supported')")
    for a in queue:
        star = "*" if a["kills_result_if_wrong"] else " "
        print(f"   {star}[{a['id']:>3}] {a['status']:>9}/{a['controversy']:<6} {a['statement'][:62]}")

    covered = set().union(*[set(a.get("params", [])) for a in A]) if A else set()
    print("\n(3) BLIND-SPOT SURFACER  (model parameters with NO catalogued assumption)")
    any_blind = False
    for fn, kind, name in SOURCES:
        params, ok = source_params(fn, kind, name)
        if not ok:
            print(f"    {fn:30s} WARNING: not found / no {name or 'constants'} parsed — NOT gap-checked")
            any_blind = True
            continue
        blind = [p for p in params if p not in covered]
        any_blind = any_blind or bool(blind)
        print(f"    {fn:30s} {len(params)-len(blind):>2}/{len(params):<2} covered | "
              f"uncovered: {', '.join(blind) if blind else '(all covered)'}")
    if any_blind:
        print("    -> each uncovered parameter is a modelling choice with no recorded justification: catalog it")
        print("       (even 'standard kinetic constant, low controversy') or remove it. That is how a blind spot dies.")
    else:
        print("    -> every model parameter maps to a catalogued assumption. No blind spots.")

    print("\n(4) HONEST UNKNOWNS / OPEN GAPS  (contested / shaky / unchecked)")
    for a in [a for a in A if a["status"] in ("contested", "shaky", "unchecked")]:
        print(f"    [{a['id']:>3}] {a['statement'][:60]}  ({a['status']}/{a['controversy']})")
    print("\n    CONCEPTUAL-FRAME LIMIT (the unknown-unknowns this tool CANNOT see):")
    print("   ", cat["_meta"]["honest_limit"])

    if args.run:
        print("\n(5) VERIFY RUN — lightweight only (4 GB cap; verify_claims.py is heavy, run separately)")
        for v in sorted(glob.glob(os.path.join(HERE, "verify_*.py"))):
            base = os.path.basename(v)
            if base in HEAVY:
                print(f"    {base:24s} SKIPPED (heavy stochastic aggregate — run capped on its own)")
                continue
            r = subprocess.run(["bash", "-c", f"ulimit -v 4194304; timeout 595 python3 {base}"],
                               capture_output=True, text=True, cwd=HERE)
            tail = (r.stdout.strip().splitlines() or ["(no output)"])[-1]
            print(f"    {base:24s} exit={r.returncode}  {tail[:44]}")

    print(f"\n==> {'BLIND SPOTS / unparsed sources present — fix above' if any_blind else 'clean: every param catalogued, no silent skips'}")
    return 1 if any_blind else 0


if __name__ == "__main__":
    raise SystemExit(main())
