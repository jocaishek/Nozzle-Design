"""
Interactive nozzle dashboard with optimization — matplotlib
Optimizes all engines at startup, then shows comparison data instantly.
Run with: python nozzle_dashboard.py
"""

import numpy as np
from scipy.optimize import fsolve, minimize
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyBboxPatch

# ── theme ───────────────────────────────────────────────────────────────────
BG       = "#0d0f12"
SURFACE  = "#13161b"
SURFACE2 = "#1a1e25"
BORDER   = "#2a2e38"
TEXT     = "#e8eaf0"
MUTED    = "#6b7280"
ACCENT   = "#4fffb0"
RED      = "#ff6b6b"
YELLOW   = "#ffd166"
BLUE     = "#00d4ff"

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


# ── physics ──────────────────────────────────────────────────────────────────
class Nozzle:
    def __init__(self, name, env, pc, tc, gamma, R, r_profile, pa, L=0.5):
        self.name      = name
        self.env       = env
        self.pc        = pc
        self.tc        = tc
        self.gamma     = gamma
        self.R         = R
        self.r_profile = np.array(r_profile, dtype=float)
        self.pa        = pa
        self.L         = L
        self._solve()

    def _area_mach_eq(self, M, epsilon):
        g = self.gamma
        return (1/M) * ((2/(g+1)) * (1 + (g-1)/2 * M**2)) ** ((g+1)/(2*(g-1))) - epsilon

    def _solve(self):
        A            = np.pi * self.r_profile**2
        self.at      = A[0]
        self.ae      = A[-1]
        self.epsilon = self.ae / self.at
        self.x       = np.linspace(0, self.L, len(self.r_profile))

        self.me   = fsolve(self._area_mach_eq, 2.5, args=(self.epsilon,))[0]
        g         = self.gamma
        self.pe   = self.pc * (1 + (g-1)/2 * self.me**2) ** (-g/(g-1))
        self.te   = self.tc / (1 + (g-1)/2 * self.me**2)
        self.ve   = self.me * np.sqrt(g * self.R * self.te)
        self.mdot = (self.at * self.pc * np.sqrt(g / (self.R * self.tc))
                     * (2/(g+1)) ** ((g+1)/(2*(g-1))))

        half_angle          = np.arctan((self.r_profile[-1] - self.r_profile[0]) / self.L)
        self.eta_div        = np.cos(half_angle)
        self.half_angle_deg = np.degrees(half_angle)
        self.thrust         = self.eta_div * self.mdot * self.ve + (self.pe - self.pa) * self.ae
        self.isp            = self.thrust / (self.mdot * 9.81)
        self.p_ratio        = self.pe / self.pa if self.pa > 0 else None

        raw_score = self.thrust / 1e6
        if self.pa > 0:
            raw_score /= (1 + abs(self.pe - self.pa) / self.pa)
        self.raw_score = raw_score

        self.mach_profile = np.ones(len(self.r_profile))
        for i in range(1, len(self.r_profile)):
            eps_i = (np.pi * self.r_profile[i]**2) / self.at
            self.mach_profile[i] = fsolve(
                lambda M, e=eps_i: (1/M)*((2/(g+1))*(1+(g-1)/2*M**2))**((g+1)/(2*(g-1)))-e,
                2.0
            )[0]

    def pressure_label(self):
        if self.p_ratio is None:
            return "vacuum", ACCENT
        if self.p_ratio < 0.9:
            return f"overexpanded  pe/pa={self.p_ratio:.2f}", RED
        if self.p_ratio > 1.1:
            return f"underexpanded  pe/pa={self.p_ratio:.2f}", YELLOW
        return f"matched  pe/pa={self.p_ratio:.2f}", ACCENT


# ── optimizer ────────────────────────────────────────────────────────────────
def optimize_engine(eng):
    r_throat = eng.r_profile[0]
    r_max    = r_throat * 50
    r_init   = eng.r_profile.copy()

    def objective(r):
        if abs(r[0] - r_throat) > 1e-6:
            return 1e9
        if np.any(np.diff(r) < 0):
            return 1e9
        if np.any(r < r_throat) or np.any(r > r_max):
            return 1e9
        try:
            tmp = Nozzle("", "", eng.pc, eng.tc, eng.gamma, eng.R, r, eng.pa, eng.L)
            return -tmp.thrust
        except Exception:
            return 1e9

    res = minimize(objective, r_init, method="Nelder-Mead",
                   options={"maxiter": 800, "xatol": 1e-7, "fatol": 1e-7})
    return Nozzle(eng.name, eng.env, eng.pc, eng.tc, eng.gamma,
                  eng.R, res.x, eng.pa, eng.L)


# ── engine definitions ───────────────────────────────────────────────────────
def ls(a, b, n=10):
    return np.linspace(a, b, n)

ENGINES = [
    Nozzle("Generic",       "SL",  3e6,    3500, 1.40, 287, ls(0.056,0.110), 101325),
    Nozzle("NASA RS-25",    "SL",  2.1e7,  3550, 1.22, 360, ls(0.056,0.150), 101325),
    Nozzle("Merlin 1D SL",  "SL",  9.7e6,  3400, 1.22, 310, ls(0.056,0.080), 101325),
    Nozzle("Merlin 1D Vac", "VAC", 9.7e6,  3400, 1.22, 310, ls(0.056,0.200), 0),
    Nozzle("Saturn F-1",    "SL",  7e6,    3300, 1.23, 300, ls(0.056,0.080), 101325),
    Nozzle("Raptor Vac",    "VAC", 3e7,    3500, 1.20, 370, ls(0.056,0.220), 0),
    Nozzle("RL-10B-2",      "VAC", 4.4e6,  3400, 1.22, 444, ls(0.056,0.280), 0),
    Nozzle("Vulcain 2",     "SL",  11.7e6, 3500, 1.22, 360, ls(0.056,0.170), 101325),
    Nozzle("RD-180",        "SL",  25.7e6, 3676, 1.20, 390, ls(0.056,0.130), 101325),
    Nozzle("RD-170",        "SL",  24.5e6, 3680, 1.20, 390, ls(0.056,0.130), 101325),
    Nozzle("Vikas",         "SL",  5.85e6, 3200, 1.23, 330, ls(0.056,0.100), 101325),
    Nozzle("BE-3",          "SL",  8.7e6,  3200, 1.26, 460, ls(0.056,0.120), 101325),
    Nozzle("Aestus",        "VAC", 1.1e6,  3050, 1.25, 300, ls(0.056,0.140), 0),
    Nozzle("Rutherford",    "SL",  8.3e6,  3300, 1.22, 320, ls(0.056,0.085), 101325),
]

print(f"Optimizing {len(ENGINES)} engines at startup...")
OPTIMIZED = []
for i, eng in enumerate(ENGINES):
    print(f"  [{i+1}/{len(ENGINES)}] {eng.name}...", end=" ", flush=True)
    opt = optimize_engine(eng)
    pct = (opt.thrust - eng.thrust) / eng.thrust * 100
    OPTIMIZED.append(opt)
    print(f"+{pct:.2f}%")

max_score = max(e.raw_score for e in ENGINES)
for e in ENGINES:
    e.norm_score = e.raw_score / max_score

max_score_opt = max(e.raw_score for e in OPTIMIZED)
for e in OPTIMIZED:
    e.norm_score = e.raw_score / max_score_opt

IMPROVEMENTS = [(o.thrust - e.thrust) / e.thrust * 100
                for e, o in zip(ENGINES, OPTIMIZED)]

print("Done. Launching dashboard...\n")


# ── dashboard ────────────────────────────────────────────────────────────────
class Dashboard:
    def __init__(self):
        self.selected = 0
        self.tab      = "overview"

        self.fig = plt.figure(figsize=(17, 9.5), facecolor=BG)
        self.fig.canvas.manager.set_window_title("Nozzle Optimization Dashboard")

        # store clickable bounding boxes for sidebar and tabs
        self._sidebar_rects = []   # list of (y_min, y_max, index)
        self._tab_rects     = []   # list of (ax, label)

        self._build_layout()
        self._draw_sidebar()
        self._draw_tabs()
        self._render()
        self.fig.canvas.mpl_connect("button_press_event", self._on_click)

    # ── layout ───────────────────────────────────────────────────────────────
    def _build_layout(self):
        outer = gridspec.GridSpec(
            1, 2, figure=self.fig,
            left=0.01, right=0.99, top=0.96, bottom=0.04,
            wspace=0.10, width_ratios=[0.17, 0.83]
        )
        self.ax_sidebar = self.fig.add_subplot(outer[0])
        self.ax_sidebar.set_facecolor(SURFACE)
        for sp in self.ax_sidebar.spines.values():
            sp.set_edgecolor(BORDER)
        self.ax_sidebar.set_xticks([]); self.ax_sidebar.set_yticks([])

        self.content_gs = gridspec.GridSpecFromSubplotSpec(
            10, 2, subplot_spec=outer[1], hspace=0.55, wspace=0.40
        )

        # tab buttons
        self.tab_axes = []
        tab_gs = gridspec.GridSpecFromSubplotSpec(
            1, 5, subplot_spec=self.content_gs[0, :], wspace=0.05
        )
        for i, label in enumerate(["Overview", "Contour", "Comparison", "Scores", "Optimization"]):
            ax = self.fig.add_subplot(tab_gs[i])
            ax.set_facecolor(SURFACE2)
            for sp in ax.spines.values():
                sp.set_edgecolor(BORDER)
            ax.set_xticks([]); ax.set_yticks([])
            ax.text(0.5, 0.5, label.upper(), ha="center", va="center",
                    fontsize=7.5, color=MUTED, transform=ax.transAxes,
                    fontfamily="monospace")
            ax.set_navigate(False)
            self.tab_axes.append((ax, label.lower()))
            self._tab_rects.append((ax, label.lower()))

        # metric cards
        self.metric_axes = []
        met_gs = gridspec.GridSpecFromSubplotSpec(
            1, 6, subplot_spec=self.content_gs[1, :], wspace=0.06
        )
        for i in range(6):
            ax = self.fig.add_subplot(met_gs[i])
            ax.set_facecolor(SURFACE)
            for sp in ax.spines.values():
                sp.set_edgecolor(BORDER)
            ax.set_xticks([]); ax.set_yticks([])
            self.metric_axes.append(ax)

        self.chart_gs = gridspec.GridSpecFromSubplotSpec(
            2, 2, subplot_spec=self.content_gs[2:, :],
            hspace=0.50, wspace=0.35
        )

    # ── sidebar ──────────────────────────────────────────────────────────────
    def _draw_sidebar(self):
        ax = self.ax_sidebar
        ax.cla()
        ax.set_facecolor(SURFACE)
        for sp in ax.spines.values():
            sp.set_edgecolor(BORDER)
        ax.set_xticks([]); ax.set_yticks([])
        n = len(ENGINES)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, n + 1)
        ax.text(0.08, n + 0.65, "ENGINES",
                color=MUTED, fontsize=7, va="center", fontfamily="monospace")

        self._sidebar_rects = []
        for i, eng in enumerate(ENGINES):
            y      = n - i - 0.5
            active = (i == self.selected)
            pct    = IMPROVEMENTS[i]
            y_lo   = (y - 0.45) / (n + 1)
            y_hi   = (y + 0.45) / (n + 1)
            self._sidebar_rects.append((y_lo, y_hi, i))

            ax.add_patch(FancyBboxPatch(
                (0.03, y_lo), 0.94, y_hi - y_lo,
                boxstyle="round,pad=0.005",
                facecolor="#1e2a22" if active else SURFACE,
                edgecolor=ACCENT if active else "none",
                linewidth=0.8, transform=ax.transAxes
            ))
            ax.text(0.08, y + 0.13, eng.name,
                    color=ACCENT if active else TEXT,
                    fontsize=7.5, va="center", fontfamily="monospace")
            ax.text(0.08, y - 0.17, f"+{pct:.2f}%",
                    color=ACCENT, fontsize=6.5, va="center", fontfamily="monospace")
            ax.text(0.93, y + 0.13, eng.env,
                    color=ACCENT if active else MUTED,
                    fontsize=6, va="center", ha="right", fontfamily="monospace")

    # ── tab bar ──────────────────────────────────────────────────────────────
    def _draw_tabs(self):
        for ax, label in self.tab_axes:
            active = (label == self.tab)
            ax.set_facecolor("#1e2a22" if active else SURFACE2)
            for sp in ax.spines.values():
                sp.set_edgecolor(ACCENT if active else BORDER)
                sp.set_linewidth(1.4 if active else 0.5)
            ax.texts[0].set_color(ACCENT if active else MUTED)

    # ── metric cards ─────────────────────────────────────────────────────────
    def _draw_metrics(self):
        eng = ENGINES[self.selected]
        opt = OPTIMIZED[self.selected]
        pct = IMPROVEMENTS[self.selected]

        cards = [
            ("THRUST",        f"{eng.thrust/1000:.1f} kN",  f"{eng.thrust:.3e} N",      MUTED,   False),
            ("OPT THRUST",    f"{opt.thrust/1000:.1f} kN",  f"+{pct:.2f}% better",      ACCENT,  True),
            ("SPEC IMPULSE",  f"{eng.isp:.0f} s",           f"ε = {eng.epsilon:.2f}",   MUTED,   False),
            ("OPT IMPULSE",   f"{opt.isp:.0f} s",           f"ε = {opt.epsilon:.2f}",   ACCENT,  True),
            ("EXIT MACH",     f"{eng.me:.2f}",              f"γ = {eng.gamma}",         MUTED,   False),
            ("EXIT PRESSURE", f"{eng.pe/1000:.1f} kPa",     eng.pressure_label()[0],    eng.pressure_label()[1], False),
        ]
        for ax, (lbl, val, sub, sc, is_opt) in zip(self.metric_axes, cards):
            ax.cla()
            ax.set_facecolor("#0d1a12" if is_opt else SURFACE)
            for sp in ax.spines.values():
                sp.set_edgecolor(ACCENT if is_opt else BORDER)
                sp.set_linewidth(1.0 if is_opt else 0.5)
            ax.set_xticks([]); ax.set_yticks([])
            ax.text(0.5, 0.84, lbl, ha="center", va="center",
                    color=ACCENT if is_opt else MUTED,
                    fontsize=6.5, transform=ax.transAxes, fontfamily="monospace")
            ax.text(0.5, 0.52, val, ha="center", va="center",
                    color=ACCENT if is_opt else TEXT,
                    fontsize=12, fontweight="bold",
                    transform=ax.transAxes, fontfamily="monospace")
            ax.text(0.5, 0.18, sub, ha="center", va="center",
                    color=sc, fontsize=6.5,
                    transform=ax.transAxes, fontfamily="monospace")

    # ── helpers ──────────────────────────────────────────────────────────────
    def _clear_charts(self):
        keep = ([self.ax_sidebar] + self.metric_axes + [a for a, _ in self.tab_axes])
        for ax in self.fig.axes:
            if ax not in keep:
                self.fig.delaxes(ax)

    def _style_ax(self, ax, title="", xlabel="", ylabel="", left_pad=0.18):
        ax.set_facecolor(SURFACE)
        for sp in ax.spines.values():
            sp.set_edgecolor(BORDER)
        ax.tick_params(colors=MUTED, labelsize=8)
        ax.grid(True, alpha=0.25)
        # extra left margin so y-tick labels don't get clipped
        ax.yaxis.set_tick_params(pad=4)
        if title:  ax.set_title(title,  color=MUTED, fontsize=8, pad=6, fontfamily="monospace")
        if xlabel: ax.set_xlabel(xlabel, color=MUTED, fontsize=8)
        if ylabel: ax.set_ylabel(ylabel, color=MUTED, fontsize=8, labelpad=6)

    def _engine_colors(self, alpha_selected=1.0, alpha_rest=0.28):
        """Return RGBA color list — one per engine, dim unselected."""
        result = []
        for i, col in enumerate(PALETTE[:len(ENGINES)]):
            r, g, b = plt.matplotlib.colors.to_rgb(col)
            a = alpha_selected if i == self.selected else alpha_rest
            result.append((r, g, b, a))
        return result

    # ── render ───────────────────────────────────────────────────────────────
    def _render(self):
        self._clear_charts()
        self._draw_metrics()
        self._draw_tabs()
        {
            "overview":     self._tab_overview,
            "contour":      self._tab_contour,
            "comparison":   self._tab_comparison,
            "scores":       self._tab_scores,
            "optimization": self._tab_optimization,
        }[self.tab]()
        self.fig.canvas.draw_idle()

    # ── overview ─────────────────────────────────────────────────────────────
    def _tab_overview(self):
        eng = ENGINES[self.selected]
        opt = OPTIMIZED[self.selected]

        ax1 = self.fig.add_subplot(self.chart_gs[0, 0])
        self._style_ax(ax1, "Mach along nozzle — original vs optimized", "x (m)", "Mach")
        ax1.plot(eng.x, eng.mach_profile, color=YELLOW, lw=1.8, ls="--", label="Original")
        ax1.plot(opt.x, opt.mach_profile, color=ACCENT,  lw=2.2, label="Optimized")
        ax1.axhline(1.0, color=MUTED, lw=0.8, ls=":", alpha=0.5)
        ax1.legend(fontsize=7, facecolor=SURFACE2, edgecolor=BORDER, labelcolor=TEXT)

        ax2 = self.fig.add_subplot(self.chart_gs[0, 1])
        self._style_ax(ax2, "Static pressure along nozzle — original vs optimized", "x (m)", "Pressure (kPa)")
        g     = eng.gamma
        p_o   = eng.pc * (1 + (g-1)/2 * eng.mach_profile**2) ** (-g/(g-1))
        g2    = opt.gamma
        p_p   = opt.pc * (1 + (g2-1)/2 * opt.mach_profile**2) ** (-g2/(g2-1))
        ax2.plot(eng.x, p_o/1000, color=YELLOW, lw=1.8, ls="--", label="Original")
        ax2.plot(opt.x, p_p/1000, color=ACCENT,  lw=2.2,          label="Optimized")
        if eng.pa > 0:
            ax2.axhline(eng.pa/1000, color=BLUE, lw=1, ls="--", alpha=0.7,
                        label=f"Pa = {eng.pa/1000:.1f} kPa")
        ax2.legend(fontsize=7, facecolor=SURFACE2, edgecolor=BORDER, labelcolor=TEXT)

        ax3 = self.fig.add_subplot(self.chart_gs[1, :])
        ax3.set_facecolor(SURFACE); ax3.set_xticks([]); ax3.set_yticks([])
        for sp in ax3.spines.values(): sp.set_edgecolor(BORDER)
        ax3.set_title("Engine parameters", color=MUTED, fontsize=8, pad=6, fontfamily="monospace")
        rows = [
            ["Chamber pressure", f"{eng.pc/1e6:.2f} MPa", "Mass flow rate", f"{eng.mdot:.2f} kg/s"],
            ["Chamber temp",     f"{eng.tc:.0f} K",        "Exit velocity",  f"{eng.ve:.0f} m/s"],
            ["Specific heat γ",  f"{eng.gamma}",           "Exit temp",      f"{eng.te:.0f} K"],
            ["Gas constant R",   f"{eng.R} J/kg·K",        "Half-angle",     f"{eng.half_angle_deg:.1f}°"],
            ["Environment",      "Vacuum" if eng.pa==0 else f"{eng.pa/1000:.1f} kPa",
             "Div. efficiency",  f"{eng.eta_div*100:.2f}%"],
        ]
        for r_idx, row in enumerate(rows):
            yy = 0.85 - r_idx * 0.17
            for c_idx, (cx, cell) in enumerate(zip([0.02,0.26,0.52,0.76], row)):
                ax3.text(cx, yy, cell, transform=ax3.transAxes,
                         color=MUTED if c_idx%2==0 else TEXT,
                         fontsize=8, va="top", fontfamily="monospace")

    # ── contour ──────────────────────────────────────────────────────────────
    def _tab_contour(self):
        eng = ENGINES[self.selected]
        opt = OPTIMIZED[self.selected]
        pct = IMPROVEMENTS[self.selected]

        ax = self.fig.add_subplot(self.chart_gs[:, :])
        self._style_ax(ax, f"Nozzle contour — {eng.name}  (original vs optimized)",
                       "Axial position (m)", "Radius (m)")

        ax.plot(eng.x,  eng.r_profile, color=YELLOW, lw=2.0, ls="--", label="Original")
        ax.plot(eng.x, -eng.r_profile, color=YELLOW, lw=2.0, ls="--")
        ax.plot(opt.x,  opt.r_profile, color=ACCENT,  lw=2.5, label="Optimized")
        ax.plot(opt.x, -opt.r_profile, color=ACCENT,  lw=2.5)
        ax.fill_between(eng.x,  eng.r_profile, -eng.r_profile, alpha=0.05, color=YELLOW)
        ax.fill_between(opt.x,  opt.r_profile, -opt.r_profile, alpha=0.09, color=ACCENT)
        ax.scatter(opt.x,  opt.r_profile, color=ACCENT, s=25, zorder=5)
        ax.scatter(opt.x, -opt.r_profile, color=ACCENT, s=25, zorder=5)
        ax.axhline(0, color=MUTED, lw=0.5, ls=":")
        ax.axvline(0, color=MUTED, lw=0.5, ls=":", alpha=0.5)
        ax.legend(fontsize=8, facecolor=SURFACE2, edgecolor=BORDER, labelcolor=TEXT)
        ax.text(0.98, 0.97, f"Thrust improvement: +{pct:.2f}%",
                ha="right", va="top", transform=ax.transAxes,
                color=ACCENT, fontsize=9, fontfamily="monospace")

    # ── comparison ───────────────────────────────────────────────────────────
    def _tab_comparison(self):
        colors = self._engine_colors()
        # darken each color for the "original" bar (same hue, less alpha)
        colors_orig = [(r, g, b, a*0.55) for r, g, b, a in colors]
        names  = [e.name for e in ENGINES]
        y      = np.arange(len(ENGINES))
        h      = 0.36

        datasets = [
            ("Thrust (kN)",         [e.thrust/1000    for e in ENGINES], [o.thrust/1000    for o in OPTIMIZED], self.chart_gs[0,0]),
            ("Specific impulse (s)",[e.isp            for e in ENGINES], [o.isp            for o in OPTIMIZED], self.chart_gs[0,1]),
            ("Exit Mach",           [e.me             for e in ENGINES], [o.me             for o in OPTIMIZED], self.chart_gs[1,0]),
            ("Expansion ratio ε",   [e.epsilon        for e in ENGINES], [o.epsilon        for o in OPTIMIZED], self.chart_gs[1,1]),
        ]

        for title, orig_data, opt_data, gs_slot in datasets:
            ax = self.fig.add_subplot(gs_slot)
            self._style_ax(ax, title)
            ax.barh(y - h/2, orig_data, h, color=colors_orig, label="Original")
            ax.barh(y + h/2, opt_data,  h, color=colors,      label="Optimized")
            ax.set_yticks(y)
            ax.set_yticklabels(names, fontsize=6.5)
            ax.tick_params(axis="y", length=0, pad=4)
            # push left edge so labels don't clip
            ax.margins(y=0.02)
            if gs_slot == self.chart_gs[0, 0]:
                ax.legend(fontsize=6.5, facecolor=SURFACE2, edgecolor=BORDER, labelcolor=TEXT)

    # ── scores ───────────────────────────────────────────────────────────────
    def _tab_scores(self):
        colors      = self._engine_colors()
        colors_orig = [(r, g, b, a*0.55) for r, g, b, a in colors]
        names       = [e.name for e in ENGINES]
        scores_orig = [e.norm_score for e in ENGINES]
        scores_opt  = [o.norm_score for o in OPTIMIZED]
        y           = np.arange(len(ENGINES))
        h           = 0.36

        ax = self.fig.add_subplot(self.chart_gs[:, :])
        self._style_ax(ax, "Normalised optimality score — original vs optimized")
        ax.barh(y - h/2, scores_orig, h, color=colors_orig, label="Original")
        ax.barh(y + h/2, scores_opt,  h, color=colors,      label="Optimized")
        ax.set_yticks(y)
        ax.set_yticklabels(names, fontsize=8)
        ax.tick_params(axis="y", length=0, pad=4)
        ax.set_xlim(0, 1.18)
        ax.axvline(1.0, color=ACCENT, lw=0.8, ls="--", alpha=0.5)
        ax.margins(y=0.02)
        for i, (so, sp) in enumerate(zip(scores_orig, scores_opt)):
            ax.text(max(so, sp) + 0.02, i, f"{so:.3f} → {sp:.3f}",
                    va="center", color=TEXT, fontsize=6.5)
        ax.legend(fontsize=8, facecolor=SURFACE2, edgecolor=BORDER, labelcolor=TEXT)

    # ── optimization ─────────────────────────────────────────────────────────
    def _tab_optimization(self):
        eng   = ENGINES[self.selected]
        opt   = OPTIMIZED[self.selected]
        pct   = IMPROVEMENTS[self.selected]
        names = [e.name for e in ENGINES]
        y     = np.arange(len(ENGINES))
        colors = self._engine_colors()

        # top-left: thrust vs Pa sweep
        ax1 = self.fig.add_subplot(self.chart_gs[0, 0])
        self._style_ax(ax1, "Thrust vs ambient pressure", "Pa (kPa)", "Thrust (kN)")
        pa_sweep = np.linspace(0, 101325, 150)

        def sweep(r_profile):
            out = []
            for pa in pa_sweep:
                tmp = Nozzle("","", eng.pc, eng.tc, eng.gamma, eng.R, r_profile, pa, eng.L)
                out.append(tmp.thrust / 1000)
            return out

        T_o = sweep(eng.r_profile)
        T_p = sweep(opt.r_profile)
        ax1.plot(pa_sweep/1000, T_o, color=YELLOW, lw=1.8, ls="--", label="Original")
        ax1.plot(pa_sweep/1000, T_p, color=ACCENT,  lw=2.2,          label="Optimized")
        ax1.fill_between(pa_sweep/1000, T_o, T_p,
                         where=[p>=o for o,p in zip(T_o,T_p)],
                         alpha=0.12, color=ACCENT)
        if eng.pa > 0:
            ax1.axvline(eng.pa/1000, color=MUTED, lw=0.8, ls=":", alpha=0.6)
        ax1.legend(fontsize=7, facecolor=SURFACE2, edgecolor=BORDER, labelcolor=TEXT)

        # top-right: % improvement bar for all engines
        ax2 = self.fig.add_subplot(self.chart_gs[0, 1])
        self._style_ax(ax2, "Thrust improvement by engine (%)", "Improvement (%)", "")
        ax2.barh(y, IMPROVEMENTS, color=colors, height=0.55)
        ax2.set_yticks(y)
        ax2.set_yticklabels(names, fontsize=6.5)
        ax2.tick_params(axis="y", length=0, pad=4)
        ax2.axvline(0, color=BORDER, lw=0.8)
        ax2.margins(y=0.02)
        for i, val in enumerate(IMPROVEMENTS):
            col = colors[i]
            ax2.text(val + 0.05, i, f"+{val:.2f}%", va="center",
                     color=col if colors[i][3] > 0.5 else MUTED,
                     fontsize=6.5, fontfamily="monospace")

        # bottom-left: contour
        ax3 = self.fig.add_subplot(self.chart_gs[1, 0])
        self._style_ax(ax3, "Contour comparison", "x (m)", "r (m)")
        ax3.plot(eng.x,  eng.r_profile, color=YELLOW, lw=1.8, ls="--", label="Original")
        ax3.plot(opt.x,  opt.r_profile, color=ACCENT,  lw=2.2,          label="Optimized")
        ax3.plot(eng.x, -eng.r_profile, color=YELLOW, lw=1.8, ls="--")
        ax3.plot(opt.x, -opt.r_profile, color=ACCENT,  lw=2.2)
        ax3.fill_between(eng.x,  eng.r_profile, -eng.r_profile, alpha=0.05, color=YELLOW)
        ax3.fill_between(opt.x,  opt.r_profile, -opt.r_profile, alpha=0.09, color=ACCENT)
        ax3.axhline(0, color=MUTED, lw=0.4, ls=":")
        ax3.legend(fontsize=7, facecolor=SURFACE2, edgecolor=BORDER, labelcolor=TEXT)

        # bottom-right: delta table
        ax4 = self.fig.add_subplot(self.chart_gs[1, 1])
        ax4.set_facecolor(SURFACE); ax4.set_xticks([]); ax4.set_yticks([])
        for sp in ax4.spines.values(): sp.set_edgecolor(BORDER)
        ax4.set_title("Optimization delta", color=MUTED, fontsize=8, pad=6, fontfamily="monospace")

        delta_rows = [
            ("Thrust",        f"{eng.thrust/1000:.2f} kN", f"{opt.thrust/1000:.2f} kN", f"+{pct:.2f}%"),
            ("Isp",           f"{eng.isp:.1f} s",          f"{opt.isp:.1f} s",          f"{opt.isp-eng.isp:+.1f} s"),
            ("Exit Mach",     f"{eng.me:.3f}",             f"{opt.me:.3f}",             f"{opt.me-eng.me:+.3f}"),
            ("Epsilon",       f"{eng.epsilon:.2f}",        f"{opt.epsilon:.2f}",        f"{opt.epsilon-eng.epsilon:+.2f}"),
            ("Half-angle",    f"{eng.half_angle_deg:.1f}°",f"{opt.half_angle_deg:.1f}°",f"{opt.half_angle_deg-eng.half_angle_deg:+.1f}°"),
            ("Exit pressure", f"{eng.pe/1000:.1f} kPa",   f"{opt.pe/1000:.1f} kPa",   f"{(opt.pe-eng.pe)/1000:+.1f} kPa"),
        ]
        col_x   = [0.02, 0.30, 0.57, 0.80]
        headers = ["Metric", "Original", "Optimized", "Delta"]
        for cx, h_lbl in zip(col_x, headers):
            ax4.text(cx, 0.93, h_lbl, transform=ax4.transAxes,
                     color=MUTED, fontsize=7.5, va="top",
                     fontfamily="monospace", fontweight="bold")
        for r_idx, (metric, orig, opt_val, delta) in enumerate(delta_rows):
            yy = 0.80 - r_idx * 0.125
            for cx, val, col in zip(col_x, [metric, orig, opt_val, delta],
                                    [MUTED,  TEXT,   ACCENT, ACCENT]):
                ax4.text(cx, yy, val, transform=ax4.transAxes,
                         color=col, fontsize=7.5, va="top", fontfamily="monospace")

    # ── click handler — pixel-space hit testing ───────────────────────────────
    def _on_click(self, event):
        if event.inaxes is None:
            return

        # tab clicks
        for ax, label in self._tab_rects:
            if event.inaxes is ax:
                if self.tab != label:
                    self.tab = label
                    self._render()
                return

        # sidebar clicks — use display coords of the sidebar axes
        if event.inaxes is self.ax_sidebar:
            ax  = self.ax_sidebar
            # convert click to axes coordinates (0-1)
            inv = ax.transAxes.inverted()
            ax_x, ax_y = inv.transform((event.x, event.y))
            for (y_lo, y_hi, idx) in self._sidebar_rects:
                if y_lo <= ax_y <= y_hi:
                    if self.selected != idx:
                        self.selected = idx
                        self._draw_sidebar()
                        self._render()
                    return

    def show(self):
        plt.show()


# ── entry point ──────────────────────────────────────────────────────────────
if __name__ == "__main__":
    dash = Dashboard()
    dash.show()