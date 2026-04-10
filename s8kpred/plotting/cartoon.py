"""
s8kpred/plotting/cartoon.py
----------------------------
Draw a secondary-structure cartoon plot using biotite.
biotite is an optional dependency — if it is not installed the function
raises ImportError with a helpful message.
"""
from __future__ import annotations

from pathlib import Path


def plot_secondary_structure_cartoon(
    sec_str: str,
    output_path: Path,
    symbols_per_line: int = 60,
    dpi: int = 300,
) -> None:
    """
    Render a helix/sheet cartoon from a secondary-structure string.

    Parameters
    ----------
    sec_str          : secondary structure string using H (helix) and E (sheet);
                       all other characters are drawn as a plain line.
    output_path      : where to save the PNG (parent directory must exist).
    symbols_per_line : residues per row in the plot.
    dpi              : image resolution.
    """
    try:
        import biotite
        import biotite.sequence as seq
        import biotite.sequence.graphics as graphics
        import matplotlib
        matplotlib.use("Agg")          # headless rendering
        import matplotlib.pyplot as plt
        import numpy as np
        from matplotlib.patches import Rectangle
    except ImportError as exc:
        raise ImportError(
            "Cartoon plotting requires biotite and matplotlib.\n"
            "Install them with:  pip install s8kpred[plot]"
        ) from exc

    # ── Custom feature plotters ───────────────────────────────────────────

    class HelixPlotter(graphics.FeaturePlotter):
        def matches(self, feature):
            return (
                feature.key == "SecStr"
                and feature.qual.get("sec_str_type") == "helix"
            )

        def draw(self, axes, feature, bbox, loc, style_param):
            n_turns = np.ceil((loc.last - loc.first + 1) / 3.6)
            x_val = np.linspace(0, n_turns * 2 * np.pi, 100)
            y_val = (-0.4 * np.sin(x_val) + 1) / 2
            x_val = x_val * bbox.width / (n_turns * 2 * np.pi) + bbox.x0
            y_val = y_val * bbox.height + bbox.y0
            bg = Rectangle(bbox.p0, bbox.width, bbox.height, color="white", linewidth=0)
            axes.add_patch(bg)
            axes.plot(x_val, y_val, linewidth=2, color="#ff4d6d")

    class SheetPlotter(graphics.FeaturePlotter):
        def __init__(self, head_width=0.8, tail_width=0.5):
            self._head_width = head_width
            self._tail_width = tail_width

        def matches(self, feature):
            return (
                feature.key == "SecStr"
                and feature.qual.get("sec_str_type") == "sheet"
            )

        def draw(self, axes, feature, bbox, loc, style_param):
            draw_head = not bool(loc.defect & seq.Location.Defect.MISS_RIGHT)
            axes.add_patch(
                biotite.AdaptiveFancyArrow(
                    bbox.x0,
                    bbox.y0 + bbox.height / 2,
                    bbox.width, 0,
                    self._tail_width * bbox.height,
                    self._head_width * bbox.height,
                    head_ratio=0.5,
                    draw_head=draw_head,
                    color="#ffc600",
                    linewidth=0,
                )
            )

    # ── Convert string → biotite Annotation ──────────────────────────────

    def _to_annotation(ss: str) -> seq.Annotation:
        features = []
        i = 0
        while i < len(ss):
            start = i
            if ss[i] == "H":
                while i < len(ss) and ss[i] == "H":
                    i += 1
                features.append(
                    seq.Feature("SecStr", [seq.Location(start, i)], {"sec_str_type": "helix"})
                )
            elif ss[i] == "E":
                while i < len(ss) and ss[i] == "E":
                    i += 1
                features.append(
                    seq.Feature("SecStr", [seq.Location(start, i)], {"sec_str_type": "sheet"})
                )
            else:
                i += 1
        return seq.Annotation(features)

    annotation = _to_annotation(sec_str)
    fig = plt.figure(figsize=(8.0, 3.0))
    ax  = fig.add_subplot(111)
    graphics.plot_feature_map(
        ax, annotation,
        symbols_per_line=symbols_per_line,
        show_numbers=True,
        show_line_position=True,
        loc_range=(1, len(sec_str) + 1),
        feature_plotters=[HelixPlotter(), SheetPlotter()],
    )
    fig.tight_layout()
    plt.savefig(str(output_path), dpi=dpi)
    plt.close(fig)
