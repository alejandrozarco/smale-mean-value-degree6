#!/usr/bin/env python3
"""Regenerate the figures in figures/ from the data in this repository (light and dark variants).
usage: python3 figures/make_figures.py        (needs matplotlib, numpy)"""
import itertools, json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "figures")
THEMES = {
    "light": dict(surface="#fcfcfb", ink="#0b0b0b", ink2="#52514e", muted="#898781", grid="#e1e0d9",
                  axis="#c3c2b7", mark="#2a78d6", mark2="#d6792a", faint="#2a78d6"),
    "dark": dict(surface="#1a1a19", ink="#ffffff", ink2="#c3c2b7", muted="#898781", grid="#2c2c2a",
                 axis="#383835", mark="#3987e5", mark2="#e58a39", faint="#3987e5"),
}
plt.rcParams.update({"font.size": 9, "svg.fonttype": "none", "svg.hashsalt": "smale-d6"})   # deterministic SVG ids


def style(ax, t, grid_y=True):
    ax.set_facecolor(t["surface"])
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(t["axis"])
    ax.tick_params(colors=t["muted"], labelsize=9, length=0)
    ax.xaxis.label.set_color(t["ink2"]); ax.yaxis.label.set_color(t["ink2"])
    if grid_y:
        ax.grid(axis="y", color=t["grid"], linewidth=0.6); ax.set_axisbelow(True)


def title(ax, t, s):
    ax.set_title(s, color=t["ink"], fontsize=10, loc="left")


def read_done(path):
    rec = {}
    for l in open(os.path.join(ROOT, path)):
        f = l.split()
        if f[0] == "R":
            rec[int(f[1])] = dict(zip(["processed", "F", "E", "L", "outside", "sym", "excl", "unres", "maxdepth"],
                                      map(int, f[2:11])))
    return rec


def read_tasks(path):
    lines = open(os.path.join(ROOT, path)).read().split("\n")[1:]
    return {int(l.split()[0]): np.array([float.fromhex(x) for x in l.split()[1:]]) for l in lines if l.strip()}


def equality_points(d):
    n = d - 1
    pts = []
    for p in itertools.permutations(range(1, n)):
        z = np.conj(np.exp(2j * np.pi * np.array(p) / n))
        pts.append(np.ravel(np.column_stack([z.real, z.imag])))
    return np.array(pts)


def save(fig, t, name):
    fig.savefig(os.path.join(OUT, name), facecolor=t["surface"], metadata={"Date": None})
    plt.close(fig)


def fig_tasks(t, name):
    """Per-task number of boxes processed, d = 5 and d = 6."""
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.6))
    fig.patch.set_facecolor(t["surface"])
    for ax, (d, run) in zip(axes, [(5, "runs/d5_v2/d5.done"), (6, "runs/d6_v2/d6.done")]):
        p = np.array([r["processed"] for r in read_done(run).values()])
        bins = np.logspace(0, np.ceil(np.log10(p.max())) + 0.01, 40)
        ax.hist(p, bins=bins, color=t["mark"], linewidth=0)
        ax.set_xscale("log"); ax.set_yscale("log")
        style(ax, t)
        ax.set_xlabel("boxes processed in one top-level task"); ax.set_ylabel("tasks")
        title(ax, t, f"d = {d}: {len(p)} tasks, {p.sum():,} boxes")
    fig.tight_layout(); save(fig, t, name)


def fig_distance(t, name):
    """d = 6: boxes per task against the distance from the task centre to the nearest equality point."""
    rec, task = read_done("runs/d6_v2/d6.done"), read_tasks("runs/d6_v2/d6.tasks")
    E = equality_points(6)
    ids = sorted(rec)
    dist = np.array([np.min(np.linalg.norm(E - task[i], axis=1)) for i in ids])
    p = np.array([rec[i]["processed"] for i in ids])
    fig, ax = plt.subplots(figsize=(10, 3.8))
    fig.patch.set_facecolor(t["surface"])
    jitter = np.random.default_rng(0).uniform(-0.012, 0.012, len(dist))
    ax.scatter(dist + jitter, p, s=3, color=t["mark"], alpha=0.35, linewidths=0)
    ax.set_yscale("log")
    style(ax, t)
    ax.set_xlabel("Euclidean distance (in R^8) from the task-box centre to the nearest of the 24 equality points "
                  "(small horizontal jitter)")
    ax.set_ylabel("boxes processed")
    title(ax, t, "d = 6: work per top-level task (half-width 1/4)")
    fig.tight_layout(); save(fig, t, name)


def fig_outcomes(t, name):
    """d = 6: how the leaves of the subdivision were closed."""
    tot = {k: 0 for k in ["F", "E", "L", "outside", "sym", "excl", "unres"]}
    for r in read_done("runs/d6_v2/d6.done").values():
        for k in tot:
            tot[k] += r[k]
    labels = [("F", "F test"), ("E", "E test"), ("L", "L test"), ("sym", "symmetry discard"),
              ("outside", "outside chart"), ("excl", "inside an excluded ball"), ("unres", "unresolved")]
    fig, ax = plt.subplots(figsize=(10, 3.2))
    fig.patch.set_facecolor(t["surface"])
    y = np.arange(len(labels))[::-1]
    vals = [tot[k] for k, _ in labels]
    ax.barh(y, [v / 1e6 for v in vals], color=t["mark"], height=0.62)
    for yi, v in zip(y, vals):
        ax.text(v / 1e6, yi, f"  {v:,}", va="center", color=t["ink2"], fontsize=9)
    ax.set_yticks(y); ax.set_yticklabels([l for _, l in labels])
    ax.set_xlim(0, max(vals) * 1.25 / 1e6)
    style(ax, t, grid_y=False); ax.grid(axis="x", color=t["grid"], linewidth=0.6); ax.set_axisbelow(True)
    ax.set_xlabel("leaf boxes (millions)")
    title(ax, t, f"d = 6: leaf boxes by outcome ({sum(vals):,} leaves)")
    fig.tight_layout(); save(fig, t, name)


def fig_local(t, name):
    """Deficit (log c - log F)/|x|^2 on E near the equality point: numerical samples and certified lower bounds."""
    S = np.loadtxt(os.path.join(ROOT, "runs/e_profile_d6.tsv"))
    cert = []
    for f in sorted(os.listdir(os.path.join(ROOT, "runs"))):
        if f.startswith("local_cert_v3_d6_T") and f.endswith(".log"):
            r = json.load(open(os.path.join(ROOT, "runs", f)))
            if r.get("certified"):
                cert.append((float(r["T"]), -r["g_upper"]))
    cert.sort()
    fig, ax = plt.subplots(figsize=(10, 4.4))
    fig.patch.set_facecolor(t["surface"])
    jit = np.random.default_rng(0).uniform(0.97, 1.03, len(S))
    ts = sorted(set(S[:, 0]))
    ax.plot(ts, [S[S[:, 0] == x, 1].min() for x in ts], color=t["mark"], linewidth=1.0,
            label="numerical: minimum over the 400 directions")
    ax.scatter(S[:, 0] * jit, S[:, 1], s=4, color=t["mark"], alpha=0.3, linewidths=0,
               label="numerical: 400 random directions per radius (local/e_profile.py, float64)")
    for i, (T, lb) in enumerate(cert):
        ax.hlines(lb, 0.004, T, color=t["mark2"], linewidth=1.6,
                  label="certified lower bound, valid for |x| ≤ T (local/local_cert.py)" if i == 0 else None)
        ax.plot([T], [lb], "o", color=t["mark2"], ms=3)
    r = 0.05
    ax.axvline(r / (1 - r), color=t["muted"], linewidth=0.8, linestyle=(0, (3, 3)))
    ax.text(r / (1 - r) * 1.03, 0.0052, "r/(1−r), r = 0.05", color=t["muted"], fontsize=8)
    ax.set_xscale("log"); ax.set_xlim(0.004, 0.25); ax.set_yscale("log"); ax.set_ylim(0.005, 0.3)
    style(ax, t)
    ax.set_xlabel("|x| (distance from the equality point in the ε-chart)")
    ax.set_ylabel("(log c − log F) / |x|²  on E")
    title(ax, t, "d = 6, near z − z⁶/6")
    leg = ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2, frameon=False, fontsize=8)
    for tx in leg.get_texts():
        tx.set_color(t["ink2"])
    fig.tight_layout(); save(fig, t, name)


if __name__ == "__main__":
    for theme, t in THEMES.items():
        fig_tasks(t, f"tasks_{theme}.svg")
        fig_distance(t, f"distance_{theme}.svg")
        fig_outcomes(t, f"outcomes_{theme}.svg")
        fig_local(t, f"local_{theme}.svg")
    print("figures written to", OUT)
