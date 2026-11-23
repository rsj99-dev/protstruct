"""二级结构一致性图（matplotlib 为可选依赖）。

图形内容：
    上：逐列主要二级结构的一致比例；
    中：全部序列按比对列排布的二级结构色块（H/E/T/C，gap 为白色）；
    右：各序列相对第一条序列的总分；
    图上圆点/方块标记二硫键、盐桥涉及的残基。

配色采用“海盐珊瑚”色系，中文使用宋体，英文与数字使用 Times New Roman。
"""

from __future__ import annotations

from typing import List, Tuple

from .analysis import AnalysisResult
from .consistency import ConsistencyReport

SS_COLORS = {
    "H": "#F4A7AE",
    "E": "#A7C7E7",
    "T": "#A6D0C6",
    "C": "#DADDE3",
    "-": "#F7FAFC",
}
SS_ORDER = ["H", "E", "T", "C"]
SS_LABELS = {
    "H": "α-螺旋 (H)",
    "E": "β-折叠 (E)",
    "T": "β-转角 (T)",
    "C": "无规卷曲 (C)",
}
_CODE = {"-": 0, "H": 1, "E": 2, "T": 3, "C": 4}


def _import_matplotlib():
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib import font_manager
        from matplotlib.colors import BoundaryNorm, ListedColormap
        from matplotlib.lines import Line2D
        from matplotlib.patches import Patch
    except ImportError as exc:  # pragma: no cover - 取决于运行环境
        raise RuntimeError(
            "生成一致性图需要 matplotlib，请先执行：pip install matplotlib"
        ) from exc
    return plt, font_manager, BoundaryNorm, ListedColormap, Line2D, Patch


def _configure_fonts(plt, font_manager) -> None:
    available = {font.name for font in font_manager.fontManager.ttflist}
    preferred = [
        "Times New Roman",
        "SimSun",
        "Songti SC",
        "Noto Serif CJK SC",
        "DejaVu Serif",
    ]
    families = [name for name in preferred if name in available] or ["DejaVu Serif"]
    plt.rcParams["font.family"] = families
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["text.color"] = "#334155"
    plt.rcParams["axes.labelcolor"] = "#334155"
    plt.rcParams["xtick.color"] = "#334155"
    plt.rcParams["ytick.color"] = "#334155"


def _ruler_step(width: int) -> int:
    if width <= 60:
        return 10
    if width <= 200:
        return 20
    if width <= 600:
        return 50
    return 100


def _shorten(label: str, limit: int = 22) -> str:
    return label if len(label) <= limit else label[: limit - 1] + "…"


def plot_consistency(
    result: AnalysisResult,
    report: ConsistencyReport,
    path: str,
    dpi: int = 200,
    figure_scale: float = 1.0,
) -> str:
    """绘制一致性图并保存为 png/svg/pdf，返回文件路径。

    figure_scale < 1 可在保持字号不变的情况下缩小画布，便于论文排版时保持
    图片内文字的可读性。
    """
    plt, font_manager, BoundaryNorm, ListedColormap, Line2D, Patch = _import_matplotlib()
    _configure_fonts(plt, font_manager)

    width = max(1, report.width)
    rows = max(1, len(report.aligned_ss))
    figure_width = min(34.0, max(11.0, 3.0 + width * 0.055 + 2.6)) * figure_scale
    figure_height = min(26.0, max(4.2, 2.4 + rows * 0.52)) * figure_scale

    fig = plt.figure(figsize=(figure_width, figure_height), facecolor="white")
    grid = fig.add_gridspec(
        2,
        2,
        height_ratios=[1.0, rows],
        width_ratios=[width, max(14.0, width * 0.22)],
        hspace=0.06,
        wspace=0.05,
    )
    ax_agreement = fig.add_subplot(grid[0, 0])
    ax_map = fig.add_subplot(grid[1, 0], sharex=ax_agreement)
    ax_score = fig.add_subplot(grid[1, 1], sharey=ax_map)

    # 上：逐列一致比例
    agreement = list(report.agreement) + [0.0] * (width - len(report.agreement))
    ax_agreement.bar(range(width), agreement, width=1.0, color="#A7C7E7", linewidth=0)
    ax_agreement.set_ylim(0, 1.0)
    ax_agreement.set_yticks([0.0, 0.5, 1.0])
    ax_agreement.set_yticklabels(["0%", "50%", "100%"], fontsize=8)
    ax_agreement.set_ylabel("逐列一致性", fontsize=8)
    ax_agreement.tick_params(axis="x", labelbottom=False, length=0)
    for spine in ("top", "right"):
        ax_agreement.spines[spine].set_visible(False)

    # 中：二级结构色块
    matrix: List[List[int]] = []
    for ss_row in report.aligned_ss:
        matrix.append([_CODE.get(char, 4) for char in ss_row])
    colormap = ListedColormap(
        [SS_COLORS["-"], SS_COLORS["H"], SS_COLORS["E"], SS_COLORS["T"], SS_COLORS["C"]]
    )
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5, 4.5], colormap.N)
    ax_map.imshow(matrix, aspect="auto", interpolation="nearest", cmap=colormap, norm=norm)

    identifiers = []
    for index, score in enumerate(report.scores):
        label = _shorten(score.identifier)
        identifiers.append(label + "（参考）" if index == 0 else label)
    ax_map.set_yticks(range(rows))
    ax_map.set_yticklabels(identifiers, fontsize=9)
    ax_map.tick_params(axis="y", length=0)
    ax_map.set_xlim(-0.5, width - 0.5)
    ax_map.set_ylim(rows - 0.5, -0.5)

    step = _ruler_step(width)
    ticks = list(range(0, width, step))
    ax_map.set_xticks(ticks)
    ax_map.set_xticklabels([str(tick + 1) for tick in ticks], fontsize=8)
    ax_map.set_xlabel("比对位置（列）", fontsize=9)
    if width <= 400:
        for column in range(10, width, 10):
            ax_map.axvline(column - 0.5, color="white", linewidth=0.4, alpha=0.65)
    ax_map.grid(axis="y", color="white", linewidth=0.8)
    for spine in ax_map.spines.values():
        spine.set_color("#A7C7E7")

    # 二硫键 / 盐桥残基标记
    disulfide_points: List[Tuple[int, int]] = []
    salt_points: List[Tuple[int, int]] = []
    for row, score in enumerate(report.scores):
        for column_a, column_b in score.disulfide_columns:
            disulfide_points.extend([(column_a, row), (column_b, row)])
        for column_a, column_b in score.salt_bridge_columns:
            salt_points.extend([(column_a, row), (column_b, row)])
    if disulfide_points:
        ax_map.scatter(
            [point[0] for point in disulfide_points],
            [point[1] for point in disulfide_points],
            s=9, c="#5878A6", marker="o", linewidths=0, zorder=3,
        )
    if salt_points:
        ax_map.scatter(
            [point[0] for point in salt_points],
            [point[1] for point in salt_points],
            s=12, c="#FBD9C9", marker="s", linewidths=0.5,
            edgecolors="#5878A6", zorder=3, alpha=0.95,
        )

    # 右：总分条形图
    totals = [score.total for score in report.scores]
    colors = ["#5878A6" if score.is_reference else "#A7C7E7" for score in report.scores]
    ax_score.barh(range(rows), totals, color=colors, height=0.62)
    ax_score.set_xlim(0, report.max_total * 1.28 if report.max_total else 1.0)
    ax_score.set_xlabel("总分（满分 %.1f）" % report.max_total, fontsize=8, labelpad=3)
    ax_score.tick_params(axis="y", labelleft=False, length=0)
    ax_score.tick_params(axis="x", labelsize=8)
    ax_score.set_ylim(rows - 0.5, -0.5)
    for row, total in enumerate(totals):
        ax_score.text(total + report.max_total * 0.02, row, "%.1f" % total, va="center", fontsize=8)
    for spine in ("top", "right"):
        ax_score.spines[spine].set_visible(False)
    ax_score.grid(axis="x", color="#DADDE3", linewidth=0.8)
    ax_score.set_axisbelow(True)

    legend_items = [
        Patch(facecolor=SS_COLORS[kind], edgecolor="#A7C7E7", label=SS_LABELS[kind])
        for kind in SS_ORDER
    ]
    legend_items.append(
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#5878A6",
               markersize=4, label="二硫键残基")
    )
    legend_items.append(
        Line2D([0], [0], marker="s", color="none", markerfacecolor="#FBD9C9",
               markeredgecolor="#5878A6", markersize=4, label="盐桥残基")
    )
    fig.legend(
        handles=legend_items,
        loc="lower center",
        ncol=6,
        frameon=False,
        fontsize=9,
        bbox_to_anchor=(0.5, 0.055),
    )

    fig.suptitle(
        "多序列二级结构一致性图（参考序列：%s）" % report.reference,
        fontsize=13,
        y=0.995,
    )
    fig.text(
        0.5,
        0.932,
        "逐列比较：二级结构相同 +1，不同 0；二硫键/盐桥对应 +0.5；参考序列为第一条序列",
        ha="center",
        fontsize=8,
        color="#5878A6",
    )
    fig.subplots_adjust(left=0.09, right=0.985, top=0.89, bottom=0.21)

    fig.savefig(path, dpi=dpi, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    return path
