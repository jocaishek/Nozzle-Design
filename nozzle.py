"""
Interactive nozzle dashboard — matplotlib
Run with: python nozzle_dashboard.py
"""

import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.widgets import Button
from matplotlib.patches import FancyBboxPatch

# ── theme ──────────────────────────────────────────────────────────────────
BG       = "#0d0f12"
SURFACE  = "#13161b"
SURFACE2 = "#1a1e25"
BORDER   = "#2a2e38"
TEXT     = "#e8eaf0"
MUTED    = "#6b7280"
ACCENT   = "#4fffb0"

PALETTE = [
    "#4fffb0","#00d4ff","#ff6b6b","#ffd166","#c77dff",
    "#06d6a0","#ef476f","#ffd60a","#4cc9f0","#f77f00",
    "#80b918","#e63946","#48cae4","#fb8500","#b5e48c",
]

plt.rcParams.update({
    "figure.facecolor":  BG,
    "axes.facecolor":    SURFACE,
    "axes.edgecolor":    BORDER,
    "axes.labelcolor":   MUTED,
    "axes.titlecolor":   TEXT,
    "xtick.color":       MUTED,
    "ytick.color":       MUTED,
    "text.color":        TEXT,
    "grid.color":        BORDER,
    "grid.linewidth":    0.5,
    "font.family":       "monospace",
    "font.size":         9,
    "axes.titlesize":    10,
    "axes.labelsize":    9,
})

# ── physics ────────────────────────────────────────────────────────────────
class Nozzle:
    def __init__(self, name, env, pc, tc, gamma, R, r_profile, pa, L=0.5):
        self.name      = name
        self.env       = env
        self.pc        = pc
        self.tc        = tc
        self.gamma     = gamma
        self.R         = R
        self.r_profile = np.array(r_profile)
        self.pa        = pa
        self.L         = L
        self._solve()

    def _area_mach_eq(self, M):
        g = self.gamma
        return (1/M) * ((2/(g+1)) * (1 + (g-1)/2 * M**2)) ** ((g+1)/(2*(g-1))) - self.epsilon

    def _solve(self):
        n          = len(self.r_profile)
        self.x     = np.linspace(0, self.L, n)
        A          = np.pi * self.r_profile**2
        self.at    = A[0]
        self.ae    = A[-1]
        self.epsilon = self.ae / self.at

        self.me = fsolve(self._area_mach_eq, 2.5)[0]

        g          = self.gamma
        self.pe    = self.pc * (1 + (g-1)/2 * self.me**2) ** (-g/(g-1))
        self.te    = self.tc / (1 + (g-1)/2 * self.me**2)
        self.ve    = self.me * np.sqrt(g * self.R * self.te)
        self.mdot  = (self.at * self.pc * np.sqrt(g / (self.R * self.tc))
                      * (2/(g+1)) ** ((g+1)/(2*(g-1))))
        half_angle     = np.arctan((self.r_profile[-1] - self.r_profile[0]) / self.L)
        self.eta_div   = np.cos(half_angle)
        self.half_angle_deg = np.degrees(half_angle)
        self.thrust    = self.eta_div * self.mdot * self.ve + (self.pe - self.pa) * self.ae
        self.isp       = self.thrust / (self.mdot * 9.81)
        self.p_ratio   = self.pe / self.pa if self.pa > 0 else None
        raw_score      = self.thrust / 1e6
        if self.pa > 0:
            raw_score /= (1 + abs(self.pe - self.pa) / self.pa)
        self.raw_score = raw_score

        # mach along profile
        self.mach_profile = np.ones(len(self.r_profile))
        for i in range(1, len(self.r_profile)):
            eps_local = (np.pi * self.r_profile[i]**2) / self.at
            self.mach_profile[i] = fsolve(
                lambda M: (1/M)*((2/(g+1))*(1+(g-1)/2*M**2))**((g+1)/(2*(g-1))) - eps_local,
                2.0
            )[0]

    def pressure_label(self):
        if self.p_ratio is None:
            return "vacuum", ACCENT
        if self.p_ratio < 0.9:
            return f"overexpanded  pe/pa={self.p_ratio:.2f}", "#ff6b6b"
        if self.p_ratio > 1.1:
            return f"underexpanded  pe/pa={self.p_ratio:.2f}", "#ffd166"
        return f"matched  pe/pa={self.p_ratio:.2f}", ACCENT


# ── engine definitions ──────────────────────────────────────────────────────
def ls(a, b, n=10):
    return np.linspace(a, b, n)

ENGINES = [
    Nozzle("Generic",       "SL",  3e6,   3500, 1.4,  287, ls(0.056,0.11),  101325),
    Nozzle("NASA RS-25",    "SL",  2.1e7, 3550, 1.22, 360, ls(0.056,0.15),  101325),
    Nozzle("Merlin 1D SL",  "SL",  9.7e6, 3400, 1.22, 310, ls(0.056,0.08),  101325),
    Nozzle("Merlin 1D Vac", "VAC", 9.7e6, 3400, 1.22, 310, ls(0.056,0.2),   0),
    Nozzle("Saturn F-1",    "SL",  7e6,   3300, 1.23, 300, ls(0.056,0.08),  101325),
    Nozzle("Raptor Vac",    "VAC", 3e7,   3500, 1.2,  370, ls(0.056,0.22),  0),
    Nozzle("RL-10B-2",      "VAC", 4.4e6, 3400, 1.22, 444, ls(0.056,0.28),  0),
    Nozzle("Vulcain 2",     "SL",  11.7e6,3500, 1.22, 360, ls(0.056,0.17),  101325),
    Nozzle("RD-180",        "SL",  25.7e6,3676, 1.2,  390, ls(0.056,0.13),  101325),
    Nozzle("RD-170",        "SL",  24.5e6,3680, 1.2,  390, ls(0.056,0.13),  101325),
    Nozzle("Vikas",         "SL",  5.85e6,3200, 1.23, 330, ls(0.056,0.1),   101325),
    Nozzle("BE-3",          "SL",  8.7e6, 3200, 1.26, 460, ls(0.056,0.12),  101325),
    Nozzle("Aestus",        "VAC", 1.1e6, 3050, 1.25, 300, ls(0.056,0.14),  0),
    Nozzle("Rutherford",    "SL",  8.3e6, 3300, 1.22, 320, ls(0.056,0.085), 101325),
]

max_score = max(e.raw_score for e in ENGINES)
for e in ENGINES:
    e.norm_score = e.raw_score / max_score


# ── dashboard ──────────────────────────────────────────────────────────────
class Dashboard:
    def __init__(self):
        self.selected = 0
        self.tab      = "overview"   # overview | contour | comparison | scores

        self.fig = plt.figure(figsize=(16, 9), facecolor=BG)
        self.fig.canvas.manager.set_window_title("Nozzle Dashboard")
        self._build_layout()
        self._draw_sidebar()
        self._draw_tabs()
        self._render()

    # ── layout ──────────────────────────────────────────────────────────────
    def _build_layout(self):
        # outer grid: sidebar | content
        outer = gridspec.GridSpec(1, 2, figure=self.fig,
                                  left=0.01, right=0.99,
                                  top=0.96, bottom=0.04,
                                  wspace=0.12,
                                  width_ratios=[0.18, 0.82])

        self.ax_sidebar = self.fig.add_subplot(outer[0])
        self.ax_sidebar.set_facecolor(SURFACE)
        for spine in self.ax_sidebar.spines.values():
            spine.set_edgecolor(BORDER)
        self.ax_sidebar.set_xticks([])
        self.ax_sidebar.set_yticks([])

        # content area subdivided later per tab
        self.content_gs = gridspec.GridSpecFromSubplotSpec(
            10, 2, subplot_spec=outer[1], hspace=0.55, wspace=0.45
        )

        # tab buttons along the top of content
        self.tab_axes = []
        tab_gs = gridspec.GridSpecFromSubplotSpec(
            1, 4, subplot_spec=self.content_gs[0, :], wspace=0.05
        )
        for i, label in enumerate(["Overview", "Contour", "Comparison", "Scores"]):
            ax = self.fig.add_subplot(tab_gs[i])
            ax.set_facecolor(SURFACE2)
            for sp in ax.spines.values():
                sp.set_edgecolor(BORDER)
            ax.set_xticks([]); ax.set_yticks([])
            ax.text(0.5, 0.5, label.upper(), ha="center", va="center",
                    fontsize=8, color=MUTED, transform=ax.transAxes,
                    fontfamily="monospace")
            ax.set_navigate(False)
            self.tab_axes.append((ax, label.lower()))

        # metric cards row
        self.metric_axes = []
        met_gs = gridspec.GridSpecFromSubplotSpec(
            1, 4, subplot_spec=self.content_gs[1, :], wspace=0.08
        )
        for i in range(4):
            ax = self.fig.add_subplot(met_gs[i])
            ax.set_facecolor(SURFACE)
            for sp in ax.spines.values():
                sp.set_edgecolor(BORDER)
            ax.set_xticks([]); ax.set_yticks([])
            self.metric_axes.append(ax)

        # main chart area (rows 2-9)
        self.chart_gs = gridspec.GridSpecFromSubplotSpec(
            2, 2, subplot_spec=self.content_gs[2:, :],
            hspace=0.5, wspace=0.35
        )

    # ── sidebar ─────────────────────────────────────────────────────────────
    def _draw_sidebar(self):
        ax = self.ax_sidebar
        ax.cla()
        ax.set_facecolor(SURFACE)
        for sp in ax.spines.values():
            sp.set_edgecolor(BORDER)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlim(0, 1)
        ax.set_ylim(0, len(ENGINES) + 1)

        ax.text(0.08, len(ENGINES) + 0.6, "ENGINES",
                color=MUTED, fontsize=7, va="center", fontfamily="monospace")

        for i, eng in enumerate(ENGINES):
            y = len(ENGINES) - i - 0.5
            active = (i == self.selected)
            bg_color = "#1e2a22" if active else SURFACE
            rect = FancyBboxPatch((0.02, y - 0.38), 0.96, 0.76,
                                  boxstyle="round,pad=0.01",
                                  facecolor=bg_color, edgecolor=ACCENT if active else "none",
                                  linewidth=0.8, transform=ax.transAxes,
                                  clip_on=False)
            # convert axes coords manually
            ax.add_patch(FancyBboxPatch(
                (0.03, (y - 0.38) / (len(ENGINES) + 1)),
                0.94, 0.72 / (len(ENGINES) + 1),
                boxstyle="round,pad=0.005",
                facecolor=bg_color, edgecolor=ACCENT if active else "none",
                linewidth=0.8, transform=ax.transAxes
            ))
            color = ACCENT if active else TEXT
            ax.text(0.08, y, eng.name, color=color, fontsize=8,
                    va="center", fontfamily="monospace")
            env_color = ACCENT if active else MUTED
            ax.text(0.92, y, eng.env, color=env_color, fontsize=6,
                    va="center", ha="right", fontfamily="monospace")

        self.fig.canvas.mpl_connect("button_press_event", self._on_click)

    # ── tab bar ─────────────────────────────────────────────────────────────
    def _draw_tabs(self):
        for ax, label in self.tab_axes:
            active = (label == self.tab)
            ax.set_facecolor("#1e2a22" if active else SURFACE2)
            for sp in ax.spines.values():
                sp.set_edgecolor(ACCENT if active else BORDER)
                sp.set_linewidth(1.2 if active else 0.5)
            ax.texts[0].set_color(ACCENT if active else MUTED)

    # ── metric cards ────────────────────────────────────────────────────────
    def _draw_metrics(self):
        eng = ENGINES[self.selected]
        labels = ["Thrust", "Spec. impulse", "Exit Mach", "Exit pressure"]
        values = [
            f"{eng.thrust/1000:.1f} kN",
            f"{eng.isp:.0f} s",
            f"{eng.me:.2f}",
            f"{eng.pe/1000:.1f} kPa",
        ]
        subs = [
            f"{eng.thrust:.3e} N",
            f"ε = {eng.epsilon:.2f}",
            f"γ = {eng.gamma}",
            eng.pressure_label()[0],
        ]
        sub_colors = [MUTED, MUTED, MUTED, eng.pressure_label()[1]]

        for ax, lbl, val, sub, sc in zip(self.metric_axes, labels, values, subs, sub_colors):
            ax.cla()
            ax.set_facecolor(SURFACE)
            for sp in ax.spines.values():
                sp.set_edgecolor(BORDER)
            ax.set_xticks([]); ax.set_yticks([])
            ax.text(0.5, 0.82, lbl.upper(), ha="center", va="center",
                    color=MUTED, fontsize=7, transform=ax.transAxes, fontfamily="monospace")
            ax.text(0.5, 0.50, val, ha="center", va="center",
                    color=TEXT, fontsize=13, fontweight="bold",
                    transform=ax.transAxes, fontfamily="monospace")
            ax.text(0.5, 0.18, sub, ha="center", va="center",
                    color=sc, fontsize=7, transform=ax.transAxes, fontfamily="monospace")

    # ── chart helpers ────────────────────────────────────────────────────────
    def _clear_charts(self):
        for ax in self.fig.axes:
            if ax not in [self.ax_sidebar] + self.metric_axes + [a for a, _ in self.tab_axes]:
                self.fig.delaxes(ax)

    def _new_ax(self, row, col, rowspan=1, colspan=1):
        if rowspan > 1 or colspan > 1:
            sub = gridspec.GridSpecFromSubplotSpec(
                2, 2, subplot_spec=self.chart_gs[:, :]
            )
            ax = self.fig.add_subplot(self.chart_gs[:, :])
        else:
            ax = self.fig.add_subplot(self.chart_gs[row, col])
        ax.set_facecolor(SURFACE)
        for sp in ax.spines.values():
            sp.set_edgecolor(BORDER)
        ax.tick_params(colors=MUTED, labelsize=8)
        ax.grid(True, alpha=0.3)
        return ax

    def _style_ax(self, ax, title="", xlabel="", ylabel=""):
        ax.set_facecolor(SURFACE)
        for sp in ax.spines.values():
            sp.set_edgecolor(BORDER)
        ax.tick_params(colors=MUTED, labelsize=8)
        ax.grid(True, alpha=0.3)
        if title:  ax.set_title(title, color=MUTED, fontsize=8, pad=6, fontfamily="monospace")
        if xlabel: ax.set_xlabel(xlabel, color=MUTED, fontsize=8)
        if ylabel: ax.set_ylabel(ylabel, color=MUTED, fontsize=8, labelpad=8)

    # ── tabs ────────────────────────────────────────────────────────────────
    def _render(self):
        self._clear_charts()
        self._draw_metrics()
        self._draw_tabs()
        if   self.tab == "overview":    self._tab_overview()
        elif self.tab == "contour":     self._tab_contour()
        elif self.tab == "comparison":  self._tab_comparison()
        elif self.tab == "scores":      self._tab_scores()
        self.fig.canvas.draw_idle()

    # overview ---------------------------------------------------------------
    def _tab_overview(self):
        eng   = ENGINES[self.selected]
        color = PALETTE[self.selected % len(PALETTE)]

        # mach profile (top left)
        ax1 = self.fig.add_subplot(self.chart_gs[0, 0])
        self._style_ax(ax1, "Mach number along nozzle", "x (m)", "Mach")
        ax1.plot(eng.x, eng.mach_profile, color=color, linewidth=2)
        ax1.fill_between(eng.x, eng.mach_profile, alpha=0.15, color=color)
        ax1.axhline(1.0, color=MUTED, linewidth=0.8, linestyle="--", alpha=0.6)

        # pressure drop (top right)
        ax2 = self.fig.add_subplot(self.chart_gs[0, 1])
        self._style_ax(ax2, "Static pressure along nozzle", "x (m)", "Pressure (kPa)")
        g  = eng.gamma
        pe_profile = eng.pc * (1 + (g-1)/2 * eng.mach_profile**2) ** (-g/(g-1))
        ax2.plot(eng.x, pe_profile / 1000, color=color, linewidth=2)
        ax2.fill_between(eng.x, pe_profile / 1000, alpha=0.15, color=color)
        if eng.pa > 0:
            ax2.axhline(eng.pa / 1000, color="#ffd166", linewidth=1,
                        linestyle="--", alpha=0.8, label=f"pa = {eng.pa/1000:.1f} kPa")
            ax2.legend(fontsize=7, facecolor=SURFACE2, edgecolor=BORDER, labelcolor=TEXT)

        # info table (bottom, spans both cols)
        ax3 = self.fig.add_subplot(self.chart_gs[1, :])
        ax3.set_facecolor(SURFACE)
        for sp in ax3.spines.values(): sp.set_edgecolor(BORDER)
        ax3.set_xticks([]); ax3.set_yticks([])
        ax3.set_title("Engine parameters", color=MUTED, fontsize=8, pad=6, fontfamily="monospace")

        rows = [
            ["Chamber pressure",  f"{eng.pc/1e6:.2f} MPa",
             "Mass flow rate",    f"{eng.mdot:.2f} kg/s"],
            ["Chamber temp",      f"{eng.tc:.0f} K",
             "Exit velocity",     f"{eng.ve:.0f} m/s"],
            ["Specific heat γ",   f"{eng.gamma}",
             "Exit temperature",  f"{eng.te:.0f} K"],
            ["Gas constant R",    f"{eng.R} J/kg·K",
             "Half-angle",        f"{eng.half_angle_deg:.1f}°"],
            ["Environment",       "Vacuum" if eng.pa == 0 else f"{eng.pa/1000:.1f} kPa",
             "Div. efficiency",   f"{eng.eta_div*100:.2f}%"],
        ]
        col_x = [0.02, 0.25, 0.52, 0.75]
        for r_idx, row in enumerate(rows):
            y = 0.85 - r_idx * 0.18
            for c_idx, cell in enumerate(row):
                color_c = MUTED if c_idx % 2 == 0 else TEXT
                ax3.text(col_x[c_idx], y, cell, transform=ax3.transAxes,
                         color=color_c, fontsize=8, va="top", fontfamily="monospace")

    # contour ----------------------------------------------------------------
    def _tab_contour(self):
        eng   = ENGINES[self.selected]
        color = PALETTE[self.selected % len(PALETTE)]

        ax = self.fig.add_subplot(self.chart_gs[:, :])
        self._style_ax(ax, f"Nozzle contour — {eng.name}", "Axial position (m)", "Radius (m)")

        ax.plot(eng.x,  eng.r_profile, color=color, linewidth=2.5, label="upper wall")
        ax.plot(eng.x, -eng.r_profile, color=color, linewidth=2.5, linestyle="--", label="lower wall")
        ax.fill_between(eng.x,  eng.r_profile, -eng.r_profile, alpha=0.08, color=color)
        ax.scatter(eng.x,  eng.r_profile, color=color, s=30, zorder=5)
        ax.scatter(eng.x, -eng.r_profile, color=color, s=30, zorder=5)
        ax.axhline(0, color=MUTED, linewidth=0.6, linestyle=":")
        ax.axvline(0, color=MUTED, linewidth=0.6, linestyle=":", label="throat")

        throat_r = eng.r_profile[0]
        ax.annotate(f"throat r={throat_r:.4f} m",
                    xy=(0, throat_r), xytext=(eng.L*0.15, throat_r * 1.4),
                    color=MUTED, fontsize=7,
                    arrowprops=dict(arrowstyle="->", color=MUTED, lw=0.8))
        exit_r = eng.r_profile[-1]
        ax.annotate(f"exit r={exit_r:.4f} m  ε={eng.epsilon:.2f}",
                    xy=(eng.L, exit_r), xytext=(eng.L*0.7, exit_r * 1.15),
                    color=MUTED, fontsize=7,
                    arrowprops=dict(arrowstyle="->", color=MUTED, lw=0.8))

        ax.legend(fontsize=7, facecolor=SURFACE2, edgecolor=BORDER, labelcolor=TEXT)

    # comparison -------------------------------------------------------------
    def _tab_comparison(self):
        names  = [e.name for e in ENGINES]
        colors = [PALETTE[i % len(PALETTE)] for i in range(len(ENGINES))]
        hi     = [1.0 if i == self.selected else 0.3 for i in range(len(ENGINES))]
        bar_colors = [(r, g, b, a) for (r, g, b), a in
                      zip([plt.matplotlib.colors.to_rgb(c) for c in colors], hi)]
        y = np.arange(len(ENGINES))

        datasets = [
            ("Thrust (kN)",       [e.thrust/1000 for e in ENGINES], self.chart_gs[0, 0]),
            ("Specific impulse (s)", [e.isp for e in ENGINES],      self.chart_gs[0, 1]),
            ("Exit Mach",         [e.me  for e in ENGINES],         self.chart_gs[1, 0]),
            ("Expansion ratio ε", [e.epsilon for e in ENGINES],     self.chart_gs[1, 1]),
        ]
        for title, data, gs_slot in datasets:
            ax = self.fig.add_subplot(gs_slot)
            self._style_ax(ax, title)
            bars = ax.barh(y, data, color=bar_colors, height=0.6)
            ax.set_yticks(y)
            ax.set_yticklabels(names, fontsize=6.5)
            ax.tick_params(axis="y", length=0)

    # scores -----------------------------------------------------------------
    def _tab_scores(self):
        names  = [e.name for e in ENGINES]
        scores = [e.norm_score for e in ENGINES]
        colors = [PALETTE[i % len(PALETTE)] for i in range(len(ENGINES))]
        hi     = [1.0 if i == self.selected else 0.3 for i in range(len(ENGINES))]
        bar_colors = [(r, g, b, a) for (r, g, b), a in
                      zip([plt.matplotlib.colors.to_rgb(c) for c in colors], hi)]
        y = np.arange(len(ENGINES))

        ax = self.fig.add_subplot(self.chart_gs[:, :])
        self._style_ax(ax, "Normalised optimality score")
        ax.barh(y, scores, color=bar_colors, height=0.6)
        ax.set_yticks(y)
        ax.set_yticklabels(names, fontsize=8)
        ax.set_xlim(0, 1.05)
        ax.axvline(1.0, color=ACCENT, linewidth=0.8, linestyle="--", alpha=0.6)
        ax.tick_params(axis="y", length=0)
        for i, (s, name) in enumerate(zip(scores, names)):
            ax.text(s + 0.01, i, f"{s:.3f}", va="center", color=TEXT, fontsize=7)

    # ── click handler ────────────────────────────────────────────────────────
    def _on_click(self, event):
        # sidebar clicks → select engine
        if event.inaxes == self.ax_sidebar:
            ax = self.ax_sidebar
            n  = len(ENGINES)
            if event.ydata is not None:
                idx = int(n - event.ydata - 0.5)
                if 0 <= idx < n:
                    self.selected = idx
                    self._draw_sidebar()
                    self._render()
            return

        # tab clicks
        for ax, label in self.tab_axes:
            if event.inaxes == ax:
                self.tab = label
                self._render()
                return

    def show(self):
        plt.show()


# ── entry point ─────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("Solving nozzles...")
    dash = Dashboard()
    dash.show()