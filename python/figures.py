#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
figures.py -- Publication figures for the SVT / CMB manuscript.

Produces three two-panel figures, each written at 300 dpi in three formats
(PNG, TIFF/LZW, PDF-vector) into the ``figures`` sub-folder.

    Figure 1  Homogenization and the noise floor
    Figure 2  Direction of entangled adjustment
    Figure 3  Power collapse

Dependencies: Python 3, numpy, pandas, matplotlib (+ Pillow, used by
matplotlib for TIFF output).  No seaborn.

Design notes
------------
* Figure width is fixed at 7.09 in (180 mm, standard double-column width).
  Figures are therefore saved *without* ``bbox_inches='tight'`` so the
  physical dimensions submitted to the journal are exact.
* Every series is separable in grayscale: colour is redundant with
  linestyle *and* marker shape.
* Fonts are embedded as TrueType (``pdf.fonttype = 42``) so the vector PDF
  remains editable at the typesetter.

Run:  python figures.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


# --------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------

DATA_DIR = Path(r"C:\Users\Salim\Desktop\makaleler\buse fidan turkon stabilizator")
OUT_DIR = DATA_DIR / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

DPI = 300
FIGSIZE = (7.09, 3.15)          # inches -- 180 mm double-column width
FORMATS = ("png", "tif", "pdf")


# --------------------------------------------------------------------------
# Palette  (Okabe-Ito colourblind-safe pair + black for reference elements)
# --------------------------------------------------------------------------

C1 = "#0072B2"      # blue   -- series 1
C2 = "#D55E00"      # vermil -- series 2
CK = "#000000"      # black  -- derived / dual-criterion series
CREF = "#8C8C8C"    # gray   -- reference lines and in-plot annotation
CGRID = "#DCDCDC"   # gray   -- recessive grid

LW = 1.4            # data line width
MS = 4.5            # marker size
MEW = 0.5           # white marker rim: keeps overlapping points separable
FS_LABEL = 10
FS_TICK = 9
FS_LEGEND = 8
FS_PANEL = 11
FS_ANNOT = 7


# --------------------------------------------------------------------------
# Global style
# --------------------------------------------------------------------------

def set_style() -> None:
    plt.rcParams.update({
        # -- typography -----------------------------------------------------
        "font.family": "serif",
        "font.serif": ["Times New Roman", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "axes.labelsize": FS_LABEL,
        "axes.titlesize": FS_LABEL,
        "xtick.labelsize": FS_TICK,
        "ytick.labelsize": FS_TICK,
        "legend.fontsize": FS_LEGEND,
        "legend.frameon": False,
        "legend.handlelength": 2.2,
        "legend.handletextpad": 0.5,
        "legend.labelspacing": 0.35,
        "legend.columnspacing": 1.1,
        "legend.borderaxespad": 0.3,
        # -- marks ----------------------------------------------------------
        "lines.linewidth": LW,
        "lines.markersize": MS,
        "lines.markeredgewidth": MEW,
        "lines.solid_capstyle": "round",
        # -- axes -----------------------------------------------------------
        "axes.linewidth": 0.7,
        "axes.edgecolor": "#333333",
        "axes.labelcolor": "black",
        "axes.axisbelow": True,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "xtick.major.size": 3.0,
        "ytick.major.size": 3.0,
        "xtick.color": "#333333",
        "ytick.color": "#333333",
        # -- output ---------------------------------------------------------
        "figure.dpi": DPI,
        "savefig.dpi": DPI,
        "pdf.fonttype": 42,     # embed TrueType, keep text selectable
        "ps.fonttype": 42,
        "pdf.compression": 6,
    })


def style_axes(ax, grid: bool = True, keep_right: bool = False) -> None:
    """Recessive frame: no top/right spine, faint horizontal grid only."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(keep_right)
    if grid:
        ax.grid(True, axis="y", color=CGRID, linewidth=0.5, linestyle="-")
    ax.tick_params(labelcolor="black")


def panel_letter(ax, letter: str) -> None:
    ax.text(0.02, 0.98, letter, transform=ax.transAxes,
            fontsize=FS_PANEL, fontweight="bold", va="top", ha="left",
            zorder=10)


def series_kw(color, marker, linestyle="-", **extra):
    """Standard keyword bundle for a data series."""
    kw = dict(color=color, marker=marker, linestyle=linestyle,
              linewidth=LW, markersize=MS,
              markerfacecolor=color, markeredgecolor="white",
              markeredgewidth=MEW, clip_on=False, zorder=3)
    kw.update(extra)
    return kw


def hline_ref(ax, y, text, x_frac=0.985, ha="right", linestyle=":", va="bottom"):
    """Gray reference line with a small in-plot label."""
    ax.axhline(y, color=CREF, linestyle=linestyle, linewidth=0.9, zorder=1)
    ax.text(x_frac, y, text, transform=ax.get_yaxis_transform(),
            fontsize=FS_ANNOT, color=CREF, ha=ha, va=va,
            bbox=dict(boxstyle="square,pad=0.12", fc="white", ec="none",
                      alpha=0.85),
            zorder=2)


# --------------------------------------------------------------------------
# Data loading + column contract
# --------------------------------------------------------------------------

EXPECTED: dict[str, list[str]] = {
    "mc_fazA_ozet.csv": ["sd_het", "rho_m", "cv0"],
    "mc_fazE.csv":      ["q", "r_xy", "cv0_med", "cv0_q25", "cv0_q75"],
    "mc_fazB.csv":      ["rho_m", "lam_zx", "boot", "cift_or", "p_poz"],
    "mc_fazB2.csv":     ["rho_m", "phi", "r_xz", "boot", "cift_or", "p_poz"],
    "mc_fazC1r.csv":    ["rho_m", "boot", "binom", "cift"],
    "mc_fazC2r.csv":    ["rho_m", "boot", "sign", "cift_sign"],
}


def load_all() -> dict[str, pd.DataFrame]:
    """Load every input file, report its shape/columns, and enforce the
    column contract.  Any missing file or column aborts the run."""
    frames: dict[str, pd.DataFrame] = {}
    problems: list[str] = []

    print("=" * 74)
    print("INPUT INSPECTION")
    print("=" * 74)

    for fname, needed in EXPECTED.items():
        path = DATA_DIR / fname
        if not path.exists():
            problems.append(f"{fname}: FILE NOT FOUND at {path}")
            print(f"\n{fname}\n  -> FILE NOT FOUND")
            continue

        df = pd.read_csv(path)
        frames[fname] = df
        missing = [c for c in needed if c not in df.columns]

        print(f"\n{fname}")
        print(f"  shape   : {df.shape}")
        print(f"  columns : {list(df.columns)}")
        print(f"  required: {needed}")
        print(f"  status  : {'OK' if not missing else 'MISSING ' + str(missing)}")

        if missing:
            problems.append(f"{fname}: missing column(s) {missing}; "
                            f"present = {list(df.columns)}")

    if problems:
        print("\n" + "=" * 74)
        print("ABORTED -- the column contract is not satisfied.")
        print("No figure was drawn; nothing was guessed or substituted.")
        print("=" * 74)
        for p in problems:
            print("  * " + p)
        sys.exit(1)

    print("\nAll expected columns present in all six files.\n")
    return frames


# --------------------------------------------------------------------------
# Saving
# --------------------------------------------------------------------------

def save_figure(fig, stem: str) -> list[Path]:
    """Save one figure as PNG, TIFF (LZW) and vector PDF.

    The two raster formats get an explicit ``dpi=300``; the PDF is written
    without a dpi argument so it stays a true vector file.
    NOTE: no ``bbox_inches='tight'`` anywhere -- that would re-crop the
    canvas and lose the exact 7.09 x 3.15 in trim size.
    """
    written: list[Path] = []

    png_path = OUT_DIR / f"{stem}.png"
    fig.savefig(png_path, dpi=300, facecolor="white")
    written.append(png_path)

    tif_path = OUT_DIR / f"{stem}.tif"
    fig.savefig(tif_path, dpi=300, facecolor="white",
                pil_kwargs={"compression": "tiff_lzw"})   # lossless
    written.append(tif_path)

    pdf_path = OUT_DIR / f"{stem}.pdf"
    fig.savefig(pdf_path, facecolor="white")              # vector, no dpi
    written.append(pdf_path)

    plt.close(fig)
    return written


# ==========================================================================
# FIGURE 1 -- Homogenization and the noise floor
# ==========================================================================

def figure1(frames) -> list[Path]:
    a = frames["mc_fazA_ozet.csv"]
    e = frames["mc_fazE.csv"]

    fig, (axa, axb) = plt.subplots(1, 2, figsize=FIGSIZE,
                                   constrained_layout=True)

    # ---- panel (a): total CV, sampling floor, extracted signal ------------
    het = (a[np.isclose(a["sd_het"], 0.18)][["rho_m", "cv0"]]
           .sort_values("rho_m"))
    flr = (a[np.isclose(a["sd_het"], 0.00)][["rho_m", "cv0"]]
           .sort_values("rho_m"))

    # align the two conditions on rho_m before differencing
    m = het.merge(flr, on="rho_m", suffixes=("_het", "_flr"))
    m["signal"] = np.sqrt(np.maximum(m["cv0_het"] ** 2 - m["cv0_flr"] ** 2, 0.0))

    axa.plot(m["rho_m"], m["cv0_het"],
             label=r"Total CV ($\sigma_h = 0.18$)",
             **series_kw(C1, "o", "-"))
    axa.plot(m["rho_m"], m["cv0_flr"],
             label=r"Sampling floor ($\sigma_h = 0$)",
             **series_kw(C2, "s", "-"))
    axa.plot(m["rho_m"], m["signal"],
             label="Extracted signal",
             **series_kw(CK, "^", "--"))

    axa.set_xlabel(r"Method-variance share $\rho_M$")
    axa.set_ylabel("Baseline CV of group coefficients (%)")
    axa.set_xlim(-0.02, 0.97)
    axa.set_ylim(0, 17.5)
    axa.set_xticks([0.0, 0.15, 0.30, 0.45, 0.60, 0.75, 0.90])
    style_axes(axa)
    axa.legend(loc="upper right")
    panel_letter(axa, "(a)")

    # ---- panel (b): straightlining -- CV shrinks while r(X,Y) inflates ----
    e = e.sort_values("q")

    axb.fill_between(e["q"], e["cv0_q25"], e["cv0_q75"],
                     color=C1, alpha=0.15, linewidth=0, zorder=1)
    l_cv, = axb.plot(e["q"], e["cv0_med"],
                     label="Median baseline CV",
                     **series_kw(C1, "o", "-"))
    p_iqr = Patch(facecolor=C1, alpha=0.15, edgecolor="none",
                  label="Interquartile range")

    axb.set_xlabel("Fraction of straightlining respondents $q$")
    axb.set_ylabel("Baseline CV (%)")
    axb.set_xlim(-0.02, 0.52)
    axb.set_ylim(0, 24)
    axb.set_xticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
    style_axes(axb, keep_right=True)

    axr = axb.twinx()
    l_r, = axr.plot(e["q"], e["r_xy"],
                    label=r"Pooled $r(X,Y)$",
                    **series_kw(C2, "s", "--"))
    axr.set_ylabel(r"Pooled correlation $r(X,Y)$")
    axr.set_ylim(0.45, 0.75)
    axr.spines["top"].set_visible(False)
    axr.tick_params(labelcolor="black")

    # one merged legend drawn from both axes
    axb.legend(handles=[l_cv, p_iqr, l_r], loc="lower left")
    panel_letter(axb, "(b)")

    return save_figure(fig, "Figure1")


# ==========================================================================
# FIGURE 2 -- Direction of entangled adjustment
# ==========================================================================

def figure2(frames) -> list[Path]:
    b = frames["mc_fazB.csv"]
    b2 = frames["mc_fazB2.csv"]

    fig, (axa, axb) = plt.subplots(1, 2, figsize=FIGSIZE,
                                   constrained_layout=True)

    rho_styles = {0.3: "-", 0.8: "--"}

    # ---- panel (a): rate vs. trait-proxy loading -------------------------
    handles_a: list[Line2D] = []
    for rho, ls in rho_styles.items():
        sub = b[np.isclose(b["rho_m"], rho)].sort_values("lam_zx")
        rho_tag = f"{rho:.1f}".lstrip("0")          # .3 / .8
        for col, colour, marker, name in (
            ("p_poz",   C1, "o", r"$p_+$"),
            ("boot",    C2, "s", "Criterion 1"),
            ("cift_or", CK, "^", "Dual criterion"),
        ):
            ln, = axa.plot(sub["lam_zx"], sub[col],
                           label=rf"{name}, $\rho_M={rho_tag}$",
                           **series_kw(colour, marker, ls))
            handles_a.append(ln)

    hline_ref(axa, 0.05, r"$\alpha=.05$", x_frac=0.985, ha="right")
    # label hangs *below* the p_+ = 1/2 line: both p_+ curves run above it
    hline_ref(axa, 0.50, r"$p_+ = 1/2$", x_frac=0.985, ha="right", va="top")

    # C3 admissibility bound -- stopped just above the p_+ curves so the
    # rule does not run through the legend text
    axa.axvline(0.359, ymin=0.0, ymax=0.824,
                color=CREF, linestyle="-.", linewidth=0.9, zorder=1)
    axa.text(0.359 - 0.012, 0.315, r"C3 bound ($\rho_M=.3$)",
             rotation=90, rotation_mode="anchor",
             fontsize=FS_ANNOT, color=CREF, ha="center", va="bottom",
             zorder=2)

    axa.set_xlabel(r"Trait-proxy loading $\lambda$")
    axa.set_ylabel("Rate / proportion")
    axa.set_xlim(-0.015, 0.615)
    axa.set_ylim(0, 0.68)
    axa.set_xticks([0.0, 0.15, 0.30, 0.45, 0.60])
    style_axes(axa)
    # Column-major ordering: left column rho_M = .3, right column rho_M = .8.
    # Anchored just right of the panel letter so the two never collide.
    axa.legend(handles=handles_a, loc="upper left",
               bbox_to_anchor=(0.075, 1.005), ncol=2, fontsize=7,
               handlelength=1.9, columnspacing=0.8)
    panel_letter(axa, "(a)")

    # ---- panel (b): rate vs. candidate entanglement r(X,Z) ---------------
    # Per-design-point label placement: the phi = .35 point is labelled above
    # its marker, the phi = .60 point below-right, so neither label lands on
    # the descending simulation line or on the field-gradient overlay.
    #                (rho_M, phi): (dx, dy, ha, va)
    phi_offsets = {
        (0.3, 0.35): (0.000,  0.030, "center", "bottom"),
        (0.3, 0.60): (0.012, -0.028, "left",   "top"),
        (0.8, 0.35): (0.000,  0.024, "center", "bottom"),
        (0.8, 0.60): (0.014, -0.028, "left",   "top"),
    }

    handles_b: list = []
    for rho, ls in rho_styles.items():
        sub = b2[np.isclose(b2["rho_m"], rho)].sort_values("r_xz")
        rho_tag = f"{rho:.1f}".lstrip("0")
        ln_p, = axb.plot(sub["r_xz"], sub["p_poz"],
                         label=rf"$p_+$, $\rho_M={rho_tag}$",
                         **series_kw(C1, "o", ls))
        ln_d, = axb.plot(sub["r_xz"], sub["cift_or"],
                         label=rf"Dual criterion, $\rho_M={rho_tag}$",
                         **series_kw(CK, "^", ls))
        handles_b += [ln_p, ln_d]

        # label each simulated design point next to its p_+ marker
        for _, row in sub.iterrows():
            phi_tag = f"{row['phi']:.2f}".lstrip("0")
            dx, dy, ha, va = phi_offsets[(round(rho, 2), round(row["phi"], 2))]
            axb.text(row["r_xz"] + dx, row["p_poz"] + dy,
                     rf"$\varphi={phi_tag}$",
                     fontsize=FS_ANNOT, color="#444444",
                     ha=ha, va=va, zorder=4)

    # field gradient overlay (Table 2, Panel C)
    x_field = [0.05, 0.15, 0.25, 0.40, 0.75]
    y_field = [0.440, 0.465, 0.423, 0.245, 0.162]
    ln_f, = axb.plot(x_field, y_field,
                     color=CK, linestyle=":", linewidth=0.9,
                     marker="D", markersize=MS, markerfacecolor=CK,
                     markeredgecolor="white", markeredgewidth=MEW,
                     label=r"Empirical gradient ($N$ = 505)",
                     clip_on=False, zorder=3)
    handles_b.append(ln_f)

    hline_ref(axb, 0.05, r"$\alpha=.05$", x_frac=0.985, ha="right")
    hline_ref(axb, 0.50, r"$p_+ = 1/2$", x_frac=0.015, ha="left")

    axb.set_xlabel(r"Candidate entanglement $r(X,Z)$")
    axb.set_ylabel("Rate / proportion")
    axb.set_xlim(0, 0.82)
    axb.set_ylim(0, 0.68)
    axb.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8])
    style_axes(axb)
    axb.legend(handles=handles_b, loc="upper right", fontsize=7)
    panel_letter(axb, "(b)")

    return save_figure(fig, "Figure2")


# ==========================================================================
# FIGURE 3 -- Power collapse
# ==========================================================================

def figure3(frames) -> list[Path]:
    c1 = frames["mc_fazC1r.csv"].sort_values("rho_m")
    c2 = frames["mc_fazC2r.csv"].sort_values("rho_m")

    fig, (axa, axb) = plt.subplots(1, 2, figsize=FIGSIZE,
                                   constrained_layout=True)

    panels = (
        (axa, "(a)", c1, "binom", "cift"),
        (axb, "(b)", c2, "sign", "cift_sign"),
    )

    for ax, letter, df, crit2_col, dual_col in panels:
        # The dual criterion is drawn first and heaviest: where it coincides
        # with Criterion 2 it reads as a halo instead of hiding it.
        ln_dual, = ax.plot(df["rho_m"], df[dual_col], label="Dual criterion",
                           **series_kw(CK, "^", "-", linewidth=2.0, zorder=3))
        ln_c1, = ax.plot(df["rho_m"], df["boot"], label="Criterion 1",
                         **series_kw(C2, "s", "-", zorder=4))
        ln_c2, = ax.plot(df["rho_m"], df[crit2_col], label="Criterion 2",
                         **series_kw(C1, "o", "-", zorder=4))

        # left edge: the right edge is where every curve converges on alpha
        hline_ref(ax, 0.05, r"$\alpha=.05$", x_frac=0.015, ha="left",
                  linestyle="--")

        ax.set_xlabel(r"Method-variance share $\rho_M$")
        ax.set_ylabel("Rejection rate")
        ax.set_xlim(-0.02, 0.92)
        ax.set_ylim(0, 0.42)            # identical in both panels
        ax.set_xticks([0.0, 0.15, 0.30, 0.45, 0.60, 0.75, 0.90])
        style_axes(ax)
        ax.legend(handles=[ln_c1, ln_c2, ln_dual], loc="upper right")
        panel_letter(ax, letter)

    return save_figure(fig, "Figure3")


# --------------------------------------------------------------------------
# DPI / integrity verification
# --------------------------------------------------------------------------

EXPECTED_PX = (round(FIGSIZE[0] * DPI), round(FIGSIZE[1] * DPI))   # 2127 x 945
PX_TOL = 2


def verify(paths: list[Path]) -> list[str]:
    """Re-open every raster with PIL, print its pixel size and info dpi, and
    check both against the target.  Returns a list of failures."""
    from PIL import Image

    failures: list[str] = []

    print("=" * 74)
    print(f"OUTPUT VERIFICATION  (target raster: "
          f"{EXPECTED_PX[0]} x {EXPECTED_PX[1]} px @ (300, 300) dpi)")
    print("=" * 74)

    for p in paths:
        size_kb = p.stat().st_size / 1024

        if not p.exists():
            failures.append(f"{p.name}: file missing")
            print(f"  {p}\n      !! FILE MISSING")
            continue

        if p.suffix.lower() in (".png", ".tif", ".tiff"):
            with Image.open(p) as im:
                size = im.size
                info_dpi = im.info.get("dpi")
                mode = im.mode

            dpi_ok = (info_dpi is not None
                      and round(float(info_dpi[0])) == DPI
                      and round(float(info_dpi[1])) == DPI)
            size_ok = (abs(size[0] - EXPECTED_PX[0]) <= PX_TOL
                       and abs(size[1] - EXPECTED_PX[1]) <= PX_TOL)

            if not dpi_ok:
                failures.append(f"{p.name}: info dpi = {info_dpi}, expected (300, 300)")
            if not size_ok:
                failures.append(f"{p.name}: size = {size}, expected {EXPECTED_PX}")

            status = "OK" if (dpi_ok and size_ok) else "!! MISMATCH"
            detail = (f"size = {size} px, info dpi = {info_dpi}, "
                      f"mode = {mode}, {status}")
        else:
            with open(p, "rb") as fh:
                head = fh.read(5)
            if head != b"%PDF-":
                failures.append(f"{p.name}: not a valid PDF")
            detail = ("vector PDF (resolution-independent), header "
                      f"{'OK' if head == b'%PDF-' else 'INVALID'}")

        print(f"  {p}")
        print(f"      {detail}  [{size_kb:,.0f} KB]")

    print()
    return failures


def main() -> None:
    set_style()
    frames = load_all()

    paths: list[Path] = []
    for name, fn in (("Figure 1", figure1),
                     ("Figure 2", figure2),
                     ("Figure 3", figure3)):
        made = fn(frames)
        print(f"{name}: wrote {', '.join(p.name for p in made)}")
        paths += made
    print()

    failures = verify(paths)

    if failures:
        print("!! VERIFICATION FAILED -- do not submit these files:")
        for f in failures:
            print("   * " + f)
        sys.exit(2)

    missing = [p for p in paths if not p.exists()]
    print(f"All {len(paths)} files verified"
          f"{'' if not missing else ' WITH MISSING FILES'} in {OUT_DIR}")


if __name__ == "__main__":
    main()
