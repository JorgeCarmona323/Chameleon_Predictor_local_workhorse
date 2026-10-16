"""Shared publication plot style for Chameleon_Predictor figures.

Import and call `apply()` at the top of any plotting script so every figure in the
project shares one typographic + color system. Single source of truth -> consistent,
publication-ready output across the pipeline.

    import plot_style; plot_style.apply()
    from plot_style import C          # role/semantic colors
"""
import matplotlib as mpl

# semantic palette (roles + neutrals); colorblind-safe-ish, legible on white
C = dict(
    hit   = "#3B6EA5",   # ordinary data point
    pc    = "#2E7D5B",   # positive control / favorable
    nc    = "#A1222A",   # negative control / unfavorable
    accent= "#C75917",   # highlight
    ink   = "#1A2B34",   # text / axes
    mut   = "#5A6B73",   # muted text / captions
    grid  = "#E3EAE9",   # faint gridlines
    band  = "#EEF3F2",   # reference band fill
    stem  = "#C6D1CF",   # lollipop stems
)

def apply():
    mpl.rcParams.update({
        # typography
        "font.family":        "sans-serif",
        "font.sans-serif":    ["Arial", "Helvetica Neue", "Helvetica", "DejaVu Sans"],
        "font.size":          11,
        "axes.titlesize":     12.5,
        "axes.titleweight":   "semibold",
        "axes.labelsize":     11,
        "axes.labelcolor":    C["ink"],
        "xtick.labelsize":    9.5,
        "ytick.labelsize":    9.5,
        "legend.fontsize":    9,
        "figure.titlesize":   13.5,
        # color of text/axes
        "text.color":         C["ink"],
        "xtick.color":        C["ink"],
        "ytick.color":        C["ink"],
        "axes.edgecolor":     C["ink"],
        # spines / lines
        "axes.linewidth":     0.9,
        "axes.spines.top":    False,
        "axes.spines.right":  False,
        "xtick.major.width":  0.9,
        "ytick.major.width":  0.9,
        "xtick.major.size":   4,
        "ytick.major.size":   4,
        # grid (call ax.grid(axis=...) to enable per-plot)
        "grid.color":         C["grid"],
        "grid.linewidth":     0.8,
        # figure / export
        "figure.facecolor":   "white",
        "axes.facecolor":     "white",
        "savefig.facecolor":  "white",
        "figure.dpi":         120,
        "savefig.dpi":        300,
        "savefig.bbox":       "tight",
        "savefig.pad_inches": 0.25,
    })
