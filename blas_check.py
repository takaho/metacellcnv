#!/usr/bin/env python3
"""Standalone diagnostic script for BLAS thread count and effective dense matmul performance.

Usage:
    python3 blas_check.py            # measure as-is
    OMP_NUM_THREADS=1 python3 blas_check.py
    OMP_NUM_THREADS=8 python3 blas_check.py

SEACells' A/B optimization cost is essentially a (k x k) @ (k x n) dense
matmul, so the GFLOP/s measured here directly determines SEACells speed.
"""
import os
import sys
import time

import numpy as np

K_DIM = int(os.environ.get("BC_K", 365))      # number of metacells
N_DIM = int(os.environ.get("BC_N", 27411))    # number of cells


def main() -> None:
    print(f"python  : {sys.version.split()[0]}")
    print(f"numpy   : {np.__version__}")
    try:
        import scipy
        print(f"scipy   : {scipy.__version__}")
    except Exception:
        pass
    print(f"CPU     : os.cpu_count()={os.cpu_count()}", end="")
    try:
        print(f" / sched_getaffinity={len(os.sched_getaffinity(0))}")
    except Exception:
        print()
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        print(f"  {var}={os.environ.get(var, '(unset)')}")

    print("\n-- BLAS linked by numpy --")
    try:
        cfg = np.show_config(mode="dicts")
        blas = cfg.get("Build Dependencies", {}).get("blas", {})
        for key in ("name", "found", "version", "detection method", "pc file directory"):
            if key in blas:
                print(f"  {key}: {blas[key]}")
    except Exception as exc:
        print(f"  np.show_config(dicts) unavailable ({exc}) -> falling back:")
        np.show_config()

    print("\n-- actual thread count seen by threadpoolctl --")
    try:
        from threadpoolctl import threadpool_info
        info = threadpool_info()
        if not info:
            print("  (nothing detected)")
        for d in info:
            print(f"  {d.get('user_api')}/{d.get('internal_api')} "
                  f"threads={d.get('num_threads')} "
                  f"prefix={d.get('prefix', '')} "
                  f"path={str(d.get('filepath', ''))[-50:]}")
    except ImportError:
        print("  threadpoolctl not installed (pip install threadpoolctl)")

    k, n = K_DIM, N_DIM
    print(f"\n-- dense matmul ({k}x{k}) @ ({k}x{n:,}) float64 --")
    rng = np.random.default_rng(0)
    t1 = np.ascontiguousarray(rng.random((k, k)))
    a = np.ascontiguousarray(rng.random((k, n)))
    flops = 2.0 * k * k * n
    t1 @ a                                  # warm-up
    ts = []
    for _ in range(5):
        s = time.perf_counter()
        t1 @ a
        ts.append(time.perf_counter() - s)
    med = float(np.median(ts))
    print(f"  {med:.3f}s  ->  {flops / med / 1e9:.1f} GFLOP/s "
          f"({flops / 1e9:.1f} GFLOP, fastest {min(ts):.3f}s)")

    print(f"\n-- reference: square float64 dgemm (2048^3) --")
    x = np.ascontiguousarray(rng.random((2048, 2048)))
    x @ x
    ts = []
    for _ in range(3):
        s = time.perf_counter()
        x @ x
        ts.append(time.perf_counter() - s)
    med2 = float(np.median(ts))
    print(f"  {med2:.3f}s  ->  {2.0 * 2048 ** 3 / med2 / 1e9:.1f} GFLOP/s")

    g = flops / med / 1e9
    print()
    if g < 20:
        print(f"[verdict] {g:.1f} GFLOP/s is low. BLAS is likely running single-threaded")
        print("       or on a reference implementation. Check the thread counts and BLAS name above.")
    elif g < 60:
        print(f"[verdict] {g:.1f} GFLOP/s. Several cores are being used, but there")
        print("       seems to be headroom relative to core count. Try raising OMP_NUM_THREADS.")
    else:
        print(f"[verdict] {g:.1f} GFLOP/s. BLAS is well utilized. "
              "Further gains would come from float32 / the line-search side.")


if __name__ == "__main__":
    main()
