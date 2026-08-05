"""
Plot Python-vs-C++ spot centroid position differences (dx=xc_py-xc_cpp,
dy=yc_py-yc_cpp) as a function of xc/yc, for a single bundle's own-selection
spot lists (pyspots.txt vs cppspots_pass4.txt -- both fiber,wave,xc,yc).

Spots are matched by (fiber, nearest wave within --wave-tol). Spots present
on only one side (selection differs slightly between the two pipelines) are
plotted at delta=0 in a distinct color/marker, positioned at their own
xc/yc, rather than dropped -- this is the point of the plot: show both the
position agreement on matched spots AND the extent of the selection
mismatch itself.

Usage:
  python testing/plot_spot_position_diff.py --py py-z8-b5.pyspots.txt \
      --cpp cpp-z8-b5_05.cppspots_pass4.txt --out spot_diff.png
"""
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load_spots(path):
    fiber, wave, xc, yc = [], [], [], []
    with open(path) as f:
        for line in f:
            p = line.strip().split(',')
            if len(p) < 4:
                continue
            fiber.append(int(float(p[0])))
            wave.append(float(p[1]))
            xc.append(float(p[2]))
            yc.append(float(p[3]))
    return np.array(fiber), np.array(wave), np.array(xc), np.array(yc)


def match(fiber_a, wave_a, xc_a, yc_a, fiber_b, wave_b, xc_b, yc_b, wave_tol):
    """Match spots in a against b by (fiber, nearest wave). Returns index
    arrays into a and b for matches, plus boolean masks of unmatched a/b."""
    matched_a, matched_b = [], []
    used_b = np.zeros(len(fiber_b), dtype=bool)
    for i in range(len(fiber_a)):
        cand = np.where((fiber_b == fiber_a[i]) & (~used_b))[0]
        if len(cand) == 0:
            continue
        j = cand[np.argmin(np.abs(wave_b[cand] - wave_a[i]))]
        if abs(wave_b[j] - wave_a[i]) <= wave_tol:
            matched_a.append(i)
            matched_b.append(j)
            used_b[j] = True
    matched_a = np.array(matched_a, dtype=int)
    matched_b = np.array(matched_b, dtype=int)
    unmatched_a = np.setdiff1d(np.arange(len(fiber_a)), matched_a)
    unmatched_b = np.setdiff1d(np.arange(len(fiber_b)), matched_b)
    return matched_a, matched_b, unmatched_a, unmatched_b


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--py", required=True, help="Python pyspots.txt (fiber,wave,xc,yc)")
    ap.add_argument("--cpp", required=True, help="C++ cppspots_pass4.txt (fiber,wave,xc,yc)")
    ap.add_argument("--out", default="spot_position_diff.png")
    ap.add_argument("--wave-tol", type=float, default=0.05, help="Angstrom tolerance for matching")
    ap.add_argument("--title", default=None)
    args = ap.parse_args()

    fp, wp, xp, yp = load_spots(args.py)
    fc, wc, xc, yc = load_spots(args.cpp)

    ma, mb, ua, ub = match(fp, wp, xp, yp, fc, wc, xc, yc, args.wave_tol)

    dx = xp[ma] - xc[mb]
    dy = yp[ma] - yc[mb]
    x_matched = xp[ma]
    y_matched = yp[ma]

    x_py_only = xp[ua]
    y_py_only = yp[ua]
    x_cpp_only = xc[ub]
    y_cpp_only = yc[ub]

    n_matched = len(ma)
    n_py_only = len(ua)
    n_cpp_only = len(ub)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

    ax = axes[0]
    ax.axhline(0, color='0.6', lw=0.8, zorder=1)
    ax.scatter(x_matched, dx, s=10, alpha=0.6, color='#2166ac',
               label=f'matched (n={n_matched})', zorder=3)
    ax.scatter(x_py_only, np.zeros_like(x_py_only), s=28, marker='^',
               color='#d6604d', label=f'python-only, {{$\\Delta$=0}} (n={n_py_only})', zorder=4)
    ax.scatter(x_cpp_only, np.zeros_like(x_cpp_only), s=28, marker='v',
               color='#4d9221', label=f'C++-only, {{$\\Delta$=0}} (n={n_cpp_only})', zorder=4)
    ax.set_xlabel('xc (pixels)')
    ax.set_ylabel(r'$\Delta x_c$ = python $-$ C++ (pixels)')
    ax.set_title('X centroid difference vs xc')
    ax.legend(fontsize=8, loc='best')
    ax.grid(alpha=0.25)

    ax = axes[1]
    ax.axhline(0, color='0.6', lw=0.8, zorder=1)
    ax.scatter(y_matched, dy, s=10, alpha=0.6, color='#2166ac',
               label=f'matched (n={n_matched})', zorder=3)
    ax.scatter(y_py_only, np.zeros_like(y_py_only), s=28, marker='^',
               color='#d6604d', label=f'python-only, {{$\\Delta$=0}} (n={n_py_only})', zorder=4)
    ax.scatter(y_cpp_only, np.zeros_like(y_cpp_only), s=28, marker='v',
               color='#4d9221', label=f'C++-only, {{$\\Delta$=0}} (n={n_cpp_only})', zorder=4)
    ax.set_xlabel('yc (pixels)')
    ax.set_ylabel(r'$\Delta y_c$ = python $-$ C++ (pixels)')
    ax.set_title('Y centroid difference vs yc')
    ax.legend(fontsize=8, loc='best')
    ax.grid(alpha=0.25)

    title = args.title or f'{args.py} vs {args.cpp}'
    fig.suptitle(f'Python$-$C++ spot centroid differences: {title}\n'
                 f'matched={n_matched}  python-only={n_py_only}  C++-only={n_cpp_only}  '
                 f'(dx rms={dx.std():.4f}px, dy rms={dy.std():.4f}px)', fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(args.out, dpi=150)
    print(f"Wrote {args.out}")
    print(f"matched={n_matched} python_only={n_py_only} cpp_only={n_cpp_only}")
    print(f"dx: mean={dx.mean():.4f} std={dx.std():.4f} max_abs={np.abs(dx).max():.4f}")
    print(f"dy: mean={dy.mean():.4f} std={dy.std():.4f} max_abs={np.abs(dy).max():.4f}")


if __name__ == "__main__":
    main()
