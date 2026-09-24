#!/usr/bin/env python3
"""
pvt_explorer.py - interactive p-v-T surface explorer for a pure substance
=========================================================================

Built for teaching the vapor dome and phase change: the three classic 2-D
property diagrams (p-v, T-v, p-T) are slices and projections of ONE 3-D
equilibrium surface, and this script lets students see that directly.

Layout
------
* Left:  the 3-D p-v-T surface, colored by phase (compressed liquid,
         liquid-vapor mixture, superheated vapor / gas, supercritical fluid),
         with the saturation lines (vapor dome) and the critical point.
* Right: the p-v, T-v and p-T diagrams.

Controls
--------
* "Hold constant": pick p, v or T.
* Slider: the value held constant.  A translucent plane cuts the 3-D surface
  and the red curve is the intersection (isobar, isochore or isotherm).
  The SAME red curve is drawn in all three 2-D diagrams.  The panel outlined
  in red is the plane the slice actually lies in; the other two panels show
  its projections.  Example: an isotherm is a curve in p-v, a horizontal
  line in T-v, and a vertical line in p-T, where the entire two-phase
  segment collapses to a single point on the saturation curve.
* "Wall projections": shadows of the dome and the red curve on the back
  walls of the 3-D box.  The 2-D diagrams are exactly these shadows.
* "Show surface": hide the surface to see the dome and slice clearly.
* Click-drag the 3-D plot to rotate it.

Equations of state
------------------
* Default: real water from CoolProp (IAPWS-95).   pip install CoolProp
  Any CoolProp fluid works:     python pvt_explorer.py --fluid R134a
                                python pvt_explorer.py --fluid CarbonDioxide
* Generic van der Waals fluid in reduced units (Maxwell equal-area rule
  for the dome):                python pvt_explorer.py --eos vdw
  This is used automatically if CoolProp is not installed.

The solid phase is not modeled: the surface starts just above the triple
point, so solid, solid-liquid and solid-vapor regions are not shown.

Other options
-------------
  --save PREFIX   render one picture per slice mode to PREFIX_p.png,
                  PREFIX_v.png, PREFIX_T.png and exit (for lecture slides)

Requirements: numpy, scipy, matplotlib (CoolProp optional).
"""

import argparse
import sys

import numpy as np
from scipy.optimize import brentq


# =============================================================================
#  Equations of state
# =============================================================================
class Fluid:
    """Common machinery. Subclasses provide Tc, pc, vc, Tmin, Tmax, pmax and
    _sat_point(T) -> (psat, vf, vg), _p_single(T, v), v_PT(p, T)."""

    # display conversions (overridden for real fluids)
    T_label, p_label, v_label = "T", "p", "v"
    T_unit, p_unit, v_unit = "", "", ""

    def T_disp(self, T):
        return T

    def T_si(self, Td):
        return Td

    def p_disp(self, p):
        return p

    def p_si(self, pd):
        return pd

    # ---- saturation table ---------------------------------------------------
    def _build_sat_table(self, n=260):
        # points clustered near the critical point, where the dome closes
        s = np.geomspace(1e-6, 1.0, n)
        Ts = np.sort(self.Tc - s * (self.Tc - self.Tmin))
        rows = []
        for T in Ts:
            try:
                ps, vf, vg = self._sat_point(T)
                if np.isfinite([ps, vf, vg]).all() and vg > vf:
                    rows.append((T, ps, vf, vg))
            except Exception:
                pass
        rows.append((self.Tc, self.pc, self.vc, self.vc))
        self.sat_T, self.sat_p, self.sat_vf, self.sat_vg = np.array(rows).T

    def psat(self, T):
        return np.exp(np.interp(T, self.sat_T, np.log(self.sat_p)))

    def vf(self, T):
        return np.exp(np.interp(T, self.sat_T, np.log(self.sat_vf)))

    def vg(self, T):
        return np.exp(np.interp(T, self.sat_T, np.log(self.sat_vg)))

    def Tsat(self, p):
        return np.interp(np.log(p), np.log(self.sat_p), self.sat_T)

    # ---- p(T, v) including the two-phase region ----------------------------
    def p_Tv(self, T, v):
        T, v = np.broadcast_arrays(np.asarray(T, float), np.asarray(v, float))
        out = np.full(T.shape, np.nan)
        two = (T < self.Tc) & (v >= self.vf(T)) & (v <= self.vg(T))
        if two.any():
            out[two] = self.psat(T[two])
        one = ~two
        if one.any():
            out[one] = self._p_single(T[one], v[one])
        return float(out) if out.ndim == 0 else out


class VanDerWaals(Fluid):
    """van der Waals fluid in reduced variables: p_r = 8T_r/(3v_r-1) - 3/v_r^2"""

    name = "van der Waals fluid (reduced units)"
    T_label, p_label, v_label = "$T_r$", "$p_r$", "$v_r$"

    def __init__(self):
        self.Tc = self.pc = self.vc = 1.0
        self.Tmin, self.Tmax, self.pmax = 0.55, 1.6, 4.0
        self._build_sat_table()
        self.vmax = 1.5 * self.vg(self.Tmin)
        self.defaults = dict(p=0.4, v=3.0, T=0.85)

    @staticmethod
    def _pr(T, v):
        return 8 * T / (3 * v - 1) - 3 / v**2

    @staticmethod
    def _roots(p, T):
        # 3p v^3 - (p + 8T) v^2 + 9 v - 3 = 0
        r = np.roots([3 * p, -(p + 8 * T), 9.0, -3.0])
        r = r[np.abs(r.imag) < 1e-9 * np.abs(r.real).max()].real
        return np.sort(r[r > 1 / 3])

    def _sat_point(self, T):
        # spinodals: 4T v^3 - 9v^2 + 6v - 1 = 0
        r = np.roots([4 * T, -9.0, 6.0, -1.0])
        r = np.sort(r[np.abs(r.imag) < 1e-12].real)
        r = r[r > 1 / 3]
        p_lo = max(self._pr(T, r[0]), 0.0)
        p_hi = self._pr(T, r[-1])
        eps = 1e-7 * (p_hi - p_lo)

        def area(p):  # (integral of p dv) - p (vg - vl); zero at saturation
            rr = self._roots(p, T)
            vl, vv = rr[0], rr[-1]
            return (8 * T / 3 * np.log((3 * vv - 1) / (3 * vl - 1))
                    + 3 * (1 / vv - 1 / vl) - p * (vv - vl))

        ps = brentq(area, p_lo + eps, p_hi - eps, xtol=1e-14)
        rr = self._roots(ps, T)
        return ps, rr[0], rr[-1]

    def _p_single(self, T, v):
        return self._pr(T, v)

    def v_PT(self, p, T):
        r = self._roots(p, T)
        if len(r) == 1 or T >= self.Tc:
            return r[-1] if T >= self.Tc else r[0]
        return r[0] if p > self.psat(T) else r[-1]


class CoolPropFluid(Fluid):
    """Real fluid from CoolProp's Helmholtz-energy equations of state."""

    T_label, p_label, v_label = "$T$", "$p$", "$v$"
    T_unit, p_unit, v_unit = "°C", "kPa", "m³/kg"

    def __init__(self, fluid="Water"):
        import CoolProp.CoolProp as CP
        self.CP = CP
        self.AS = CP.AbstractState("HEOS", fluid)
        self.name = f"{fluid} (CoolProp)"
        self.Tc = self.AS.T_critical()
        self.pc = self.AS.p_critical()
        self.vc = 1.0 / self.AS.rhomass_critical()
        self.Tmin = max(self.AS.Ttriple(), self.AS.Tmin()) + 1.0
        self.Tmax = 1.4 * self.Tc
        self.pmax = 4.5 * self.pc
        # Some fluids (e.g. CO2) would be solid at (Tmin, pmax); since the solid
        # isn't modeled, nudge Tmin up until the liquid exists there.
        while not np.isfinite(self.v_PT(self.pmax, self.Tmin)) and self.Tmin < 0.9 * self.Tc:
            self.Tmin += 1.0
        self._build_sat_table()
        self.vmax = 1.5 * self.vg(self.Tmin)
        if fluid.lower() == "water":   # textbook-friendly starting values
            self.defaults = dict(p=1.0e6, v=0.05, T=473.15)
        else:
            self.defaults = dict(p=0.045 * self.pc, v=5 * self.vc,
                                 T=self.Tmin + 0.55 * (self.Tc - self.Tmin))

    def T_disp(self, T):
        return np.asarray(T) - 273.15

    def T_si(self, Td):
        return np.asarray(Td) + 273.15

    def p_disp(self, p):
        return np.asarray(p) / 1e3

    def p_si(self, pd):
        return np.asarray(pd) * 1e3

    def _sat_point(self, T):
        self.AS.update(self.CP.QT_INPUTS, 0.0, T)
        ps, vf = self.AS.p(), 1.0 / self.AS.rhomass()
        self.AS.update(self.CP.QT_INPUTS, 1.0, T)
        return ps, vf, 1.0 / self.AS.rhomass()

    def _p_single(self, T, v):
        out = np.empty(T.size)
        for i, (Ti, vi) in enumerate(zip(T.ravel(), v.ravel())):
            try:
                self.AS.update(self.CP.DmassT_INPUTS, 1.0 / vi, Ti)
                out[i] = self.AS.p()
            except Exception:
                out[i] = np.nan
        return out.reshape(T.shape)

    def v_PT(self, p, T):
        for dT in (0.0, 1e-4, -1e-4):
            try:
                self.AS.update(self.CP.PT_INPUTS, p, T + dT)
                return 1.0 / self.AS.rhomass()
            except Exception:
                pass
        return np.nan


# =============================================================================
#  Slices through the surface (SI units in, SI units out)
# =============================================================================
def isotherm(fl, T, n=300):
    v = np.geomspace(fl.v_PT(fl.pmax, T), fl.vmax, n)
    if T < fl.Tc:
        v = np.sort(np.concatenate([v, [fl.vf(T), fl.vg(T)]]))
    return v, np.full_like(v, T), fl.p_Tv(T, v)


def isochore(fl, v, n=300):
    T = np.linspace(fl.Tmin, fl.Tmax, n)
    Tx = dome_crossing_T(fl, v)
    if Tx is not None:
        T = np.sort(np.append(T, Tx))
    p = fl.p_Tv(T, v)
    p[p > fl.pmax * 1.0001] = np.nan
    return np.full_like(T, v), T, p


def isobar(fl, p, n=220):
    v1 = fl.v_PT(p, fl.Tmin)
    v2 = min(fl.v_PT(p, fl.Tmax), fl.vmax)
    v = np.geomspace(v1, v2, n)
    below = p < fl.pc
    if below:
        Ts = float(fl.Tsat(p))
        vf, vg = fl.vf(Ts), fl.vg(Ts)
        v = np.sort(np.concatenate([v, [x for x in (vf, vg) if v1 < x < v2]]))
    T = np.full_like(v, np.nan)
    for i, vi in enumerate(v):
        if below and vf <= vi <= vg:
            T[i] = Ts
            continue
        fa = fl.p_Tv(fl.Tmin, vi) - p
        fb = fl.p_Tv(fl.Tmax, vi) - p
        if abs(fa) < 1e-6 * p:
            T[i] = fl.Tmin
        elif abs(fb) < 1e-6 * p:
            T[i] = fl.Tmax
        elif np.isfinite(fa) and np.isfinite(fb) and fa * fb < 0:
            try:
                T[i] = brentq(lambda t: fl.p_Tv(t, vi) - p,
                              fl.Tmin, fl.Tmax, xtol=1e-7 * fl.Tc)
            except Exception:
                pass
    return v, T, np.full_like(v, p)


def dome_crossing_T(fl, v):
    """Temperature at which a constant-v (rigid tank) heating line meets the
    saturation curve, or None if it never does in the plotted range."""
    arr = fl.sat_vf if v < fl.vc else fl.sat_vg
    d = np.log(arr) - np.log(v)
    idx = np.where(np.sign(d[:-1]) != np.sign(d[1:]))[0]
    if len(idx) == 0:
        return None
    k = idx[-1]
    f = d[k] / (d[k] - d[k + 1])
    return fl.sat_T[k] + f * (fl.sat_T[k + 1] - fl.sat_T[k])


# =============================================================================
#  Look and feel
# =============================================================================
BG = "#F4F5F7"          # figure background
CARD = "#FFFFFF"        # panel background
EDGE = "#D9DCE2"        # panel borders
INK = "#1F2430"         # main text / dome
MUTED = "#6B7280"       # secondary text
FAINT = "#C9CED6"       # background curves
ACCENT = "#D6336C"      # the selected slice
ACCENT_TINT = "#FCF1F5"

PHASE_COLORS = {
    "Compressed liquid": "#3A78C9",
    "Liquid + vapor": "#5DB28C",
    "Superheated vapor / gas": "#EE9A55",
    "Supercritical fluid": "#A3AABA",
}
PHASE_TEXT = {  # darker versions for labels
    "liquid": "#2A5C9E", "mix": "#3C8A66", "vapor": "#C06A25", "super": "#6B7285",
}


def apply_style():
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Inter", "Helvetica Neue", "Helvetica", "Arial",
                            "Segoe UI", "DejaVu Sans"],
        "mathtext.fontset": "dejavusans",
        "font.size": 9.5,
        "text.color": INK,
        "axes.edgecolor": EDGE,
        "axes.labelcolor": INK,
        "axes.labelsize": 9.5,
        "axes.titlesize": 10.5,
        "axes.titleweight": "bold",
        "axes.facecolor": CARD,
        "axes.linewidth": 0.9,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "xtick.minor.size": 1.5,
        "ytick.minor.size": 1.5,
        "grid.color": "#E6E8EC",
        "grid.linewidth": 0.7,
        "figure.facecolor": BG,
        "savefig.facecolor": BG,
        "toolbar": "toolbar2",
    })
    import logging
    logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)


# =============================================================================
#  The 3-D surface
# =============================================================================
def build_surface(fl, nT=61, nL=22, n2=22, nV=40):
    """Grid whose columns follow the phase boundaries, so the dome edges are
    crisp: each isotherm row = [liquid | two-phase | vapor] segments."""
    from matplotlib.colors import to_rgba
    Ts = np.unique(np.append(np.linspace(fl.Tmin, fl.Tmax, nT), fl.Tc))
    colors = {k: np.array(to_rgba(c)) for k, c in PHASE_COLORS.items()}
    V, P, C = [], [], []
    for T in Ts:
        vlow = fl.v_PT(fl.pmax, T)
        if T < fl.Tc:
            a, b = fl.vf(T), fl.vg(T)
        else:
            a = b = max(fl.vc, vlow * 1.01)
        a = max(a, vlow * 1.0001)
        b = min(max(b, a), fl.vmax / 1.001)
        v = np.concatenate([np.geomspace(vlow, a, nL)[:-1],
                            np.geomspace(a, b, n2),
                            np.geomspace(b, fl.vmax, nV)[1:]])
        p = fl.p_Tv(T, v)
        if T < fl.Tc:
            p[nL - 1:nL - 1 + n2] = fl.psat(T)
        c = np.empty((len(v), 4))
        for j in range(len(v)):
            pm = np.nanmean(p[j:j + 2])
            if T < fl.Tc:
                if j < nL - 1:
                    key = "Compressed liquid"
                elif j < nL - 1 + n2 - 1:
                    key = "Liquid + vapor"
                else:
                    key = "Superheated vapor / gas"
            else:
                # soft blend across p = p_c, so the boundary isn't a staircase
                t = np.clip(np.log(pm / fl.pc) / 0.25 + 0.5, 0, 1)
                c[j] = (t * colors["Supercritical fluid"]
                        + (1 - t) * colors["Superheated vapor / gas"])
                continue
            c[j] = colors[key]
        V.append(v)
        P.append(p)
        C.append(c)
    V, P, C = np.array(V), np.array(P), np.array(C)
    return V, np.repeat(Ts[:, None], V.shape[1], axis=1), P, C


def shade(X, Y, Z, C, light=(0.35, -0.5, 0.8), ambient=0.58):
    """Soft Lambert shading computed in normalized (box) coordinates, so the
    very different axis scales don't distort the lighting."""
    def norm(A):
        return (A - np.nanmin(A)) / (np.nanmax(A) - np.nanmin(A))
    x, y, z = norm(X), norm(Y), norm(Z)
    du = np.stack([np.gradient(a, axis=1) for a in (x, y, z)], -1)
    dv = np.stack([np.gradient(a, axis=0) for a in (x, y, z)], -1)
    n = np.cross(du, dv)
    mag = np.linalg.norm(n, axis=-1, keepdims=True)
    with np.errstate(invalid="ignore", divide="ignore"):
        n = n / mag
    L = np.asarray(light) / np.linalg.norm(light)
    lam = np.abs(np.nan_to_num(n @ L, nan=0.8))
    k = ambient + (1 - ambient) * lam
    out = C.copy()
    out[..., :3] = np.clip(out[..., :3] * k[..., None] + 0.04, 0, 1)
    return out


# =============================================================================
#  The interactive figure
# =============================================================================
MODES = ["p  (isobar)", "v  (isochore)", "T  (isotherm)"]


class Explorer:
    def __init__(self, fl):
        import textwrap
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D
        from matplotlib.patches import FancyBboxPatch
        from matplotlib.ticker import FuncFormatter
        from matplotlib.widgets import CheckButtons, RadioButtons, Slider

        apply_style()
        self.plt, self.fl, self.wrap = plt, fl, textwrap
        self.mode = "T"
        self.dyn, self.wall_dyn = [], []

        print("Building surface ...", flush=True)
        V, T, P, C = build_surface(fl)
        self.X = np.log10(V)
        self.Y = fl.T_disp(T)
        self.Z = np.log10(fl.p_disp(P))
        dT, dp = fl.T_disp, fl.p_disp
        self.dome_v = np.concatenate([fl.sat_vf, fl.sat_vg[::-1]])
        self.dome_T = np.concatenate([fl.sat_T, fl.sat_T[::-1]])
        self.dome_p = np.concatenate([fl.sat_p, fl.sat_p[::-1]])

        # ---------- figure & header ---------------------------------------
        self.fig = fig = plt.figure(figsize=(15.5, 9.2))
        try:
            fig.canvas.manager.set_window_title("p-v-T Surface Explorer")
        except Exception:
            pass
        fig.text(0.02, 0.968, "p–v–T Surface Explorer", fontsize=17,
                 weight="bold", color=INK)
        fig.text(0.02, 0.940, f"{fl.name}   ·   drag the 3-D plot to rotate   ·   "
                 "hold one property constant to slice the surface",
                 fontsize=9.5, color=MUTED)
        handles = [Line2D([], [], marker="s", ls="", ms=10, mfc=c, mec="none",
                          label=k) for k, c in PHASE_COLORS.items()]
        handles += [Line2D([], [], color=INK, lw=2, label="Saturation (vapor dome)"),
                    Line2D([], [], color=ACCENT, lw=2.5, label="Selected slice")]
        fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.015, 0.93),
                   ncol=3, frameon=False, fontsize=9, handletextpad=0.4,
                   columnspacing=1.4)

        # ---------- 3-D axes ------------------------------------------------
        ax = self.ax3 = fig.add_axes([-0.03, 0.175, 0.62, 0.725], projection="3d")
        ax.set_facecolor("none")
        ax.computed_zorder = False
        for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
            axis.set_pane_color((1, 1, 1, 0.55))
            axis.pane.set_edgecolor(EDGE)
            axis.line.set_color(EDGE)

        xl = (np.nanmin(self.X) - 0.35, np.nanmax(self.X) + 0.05)
        yl = (np.nanmin(self.Y), np.nanmax(self.Y))
        zl = (np.nanmin(self.Z) - 0.2, np.nanmax(self.Z) + 0.05)
        self.lims = xl, yl, zl
        xw, yw, zw = xl[0], yl[1], zl[0]

        # wall shadows of the dome (drawn first = behind everything)
        lv, lT, lp = np.log10(self.dome_v), dT(self.dome_T), np.log10(dp(self.dome_p))
        kw = dict(color="#8C93A1", lw=1.2, zorder=0)
        self.wall_static = [
            ax.plot(lv, np.full_like(lv, yw), lp, **kw)[0],
            ax.plot(lv, lT, np.full_like(lv, zw), **kw)[0],
            ax.plot(np.full_like(lT, xw), lT, lp, **kw)[0],
        ]

        FC = shade(self.X, self.Y, self.Z, C)
        surf = ax.plot_surface(self.X, self.Y, self.Z, facecolors=FC,
                               rstride=1, cstride=1, linewidth=0.3, shade=False,
                               antialiased=True, zorder=1)
        try:   # edges the same color as faces hides hairline seams between cells
            surf.set_edgecolor(surf._facecolor3d)
        except Exception:
            surf.set_linewidth(0)
        self.surf = [surf]
        # a few isotherms drawn on the surface as subtle contour lines
        for i in range(0, self.X.shape[0], 6):
            self.surf.append(ax.plot(self.X[i], self.Y[i], self.Z[i], color="white",
                                     lw=0.6, alpha=0.55, zorder=2)[0])
        ax.plot(lv, lT, lp, color=INK, lw=2.2, zorder=4)
        self.cp3 = (np.log10(fl.vc), float(dT(fl.Tc)), np.log10(float(dp(fl.pc))))
        ax.scatter(*self.cp3, s=36, color="white", edgecolor=INK, linewidth=1.5,
                   depthshade=False, zorder=6)

        ax.set_xlim(xl); ax.set_ylim(yl); ax.set_zlim(zl)
        pow10 = FuncFormatter(lambda x, _: f"$10^{{{int(round(x))}}}$")
        ax.set_xticks(np.arange(np.ceil(xl[0]), np.floor(xl[1]) + 1))
        ax.set_zticks(np.arange(np.ceil(zl[0]), np.floor(zl[1]) + 1))
        ax.xaxis.set_major_formatter(pow10)
        ax.zaxis.set_major_formatter(pow10)
        ax.tick_params(labelsize=8.5, colors=MUTED, pad=1)
        ax.set_xlabel(self._lab("v"), labelpad=7)
        ax.set_ylabel(self._lab("T"), labelpad=7)
        ax.set_zlabel(self._lab("p"), labelpad=5)
        try:
            ax.set_box_aspect((1.25, 1.0, 0.85), zoom=1.08)
        except TypeError:
            ax.set_box_aspect((1.25, 1.0, 0.85))
        ax.view_init(elev=22, azim=-58)

        # ---------- 2-D panels ----------------------------------------------
        w, h, x0 = 0.37, 0.195, 0.608
        self.ax_pv = fig.add_axes([x0, 0.745, w, h])
        self.ax_Tv = fig.add_axes([x0, 0.475, w, h])
        self.ax_pT = fig.add_axes([x0, 0.205, w, h])
        self._setup_2d()

        # ---------- bottom cards -------------------------------------------
        def card(x, y, cw, ch):
            fig.add_artist(FancyBboxPatch((x, y), cw, ch, transform=fig.transFigure,
                                          boxstyle="round,pad=0,rounding_size=0.008",
                                          fc=CARD, ec=EDGE, lw=0.9, zorder=-1))
        card(0.015, 0.018, 0.56, 0.15)
        card(0.608, 0.018, 0.37, 0.15)

        def caps(x, y, s):
            fig.text(x, y, s, fontsize=8, weight="bold", color=MUTED)

        caps(0.03, 0.143, "HOLD CONSTANT")
        caps(0.165, 0.143, "DISPLAY")
        caps(0.622, 0.143, "WHAT'S HAPPENING")

        rax = fig.add_axes([0.025, 0.028, 0.13, 0.11])
        cax = fig.add_axes([0.16, 0.058, 0.12, 0.08])
        for a in (rax, cax):
            a.set_facecolor("none")
            for sp in a.spines.values():
                sp.set_visible(False)
        try:
            self.radio = RadioButtons(
                rax, MODES, active=2, activecolor=ACCENT,
                label_props={"fontsize": [10] * 3, "color": [INK] * 3},
                radio_props={"s": [70] * 3, "edgecolor": ["#9AA0AB"] * 3})
        except TypeError:   # matplotlib < 3.7
            self.radio = RadioButtons(rax, MODES, active=2, activecolor=ACCENT)
        try:
            self.checks = CheckButtons(
                cax, ["Wall projections", "Show surface"], [True, True],
                label_props={"fontsize": [10] * 2, "color": [INK] * 2},
                frame_props={"s": [70] * 2, "edgecolor": ["#9AA0AB"] * 2},
                check_props={"color": [ACCENT] * 2, "s": [45] * 2})
        except TypeError:
            self.checks = CheckButtons(cax, ["Wall projections", "Show surface"], [True, True])

        self.slider_title = fig.text(0.31, 0.143, "", fontsize=8, weight="bold", color=MUTED)
        sax = fig.add_axes([0.31, 0.085, 0.2, 0.03])
        sax.set_facecolor("none")
        self.slider = Slider(sax, "", 0, 1, valinit=0.5, color=ACCENT,
                             track_color="#E4E6EA",
                             handle_style={"facecolor": "white", "edgecolor": ACCENT,
                                           "size": 13})
        self.slider.vline.set_visible(False)
        self.slider.valtext.set_visible(False)
        self.value_text = fig.text(0.31, 0.028, "", fontsize=15, weight="bold", color=ACCENT)
        self.range_text = fig.text(0.51, 0.064, "", fontsize=8, color=MUTED, ha="right")

        self.status = fig.text(0.622, 0.128, "", fontsize=9, va="top", color=INK,
                               linespacing=1.45)

        self.radio.on_clicked(self._on_mode)
        self.checks.on_clicked(self._on_check)
        self.slider.on_changed(self._update)
        self._on_mode(MODES[2])

    # ------------------------------------------------------------------
    def _lab(self, q):
        fl = self.fl
        name = {"p": fl.p_label, "v": fl.v_label, "T": fl.T_label}[q]
        unit = {"p": fl.p_unit, "v": fl.v_unit, "T": fl.T_unit}[q]
        return f"{name}  [{unit}]" if unit else name

    def _setup_2d(self):
        fl, dT, dp = self.fl, self.fl.T_disp, self.fl.p_disp
        g = dict(color=FAINT, lw=0.8, zorder=1)
        dome = dict(color=INK, lw=2, zorder=3)

        Tbg = np.unique(np.concatenate([
            np.linspace(fl.Tmin, fl.Tc, 5)[1:-1], [fl.Tc],
            np.linspace(fl.Tc, fl.Tmax, 4)[1:]]))
        pbg = np.unique(np.append(np.geomspace(fl.psat(fl.Tmin) * 3, fl.pmax / 1.5, 6), fl.pc))
        vbg = np.unique(np.append(np.geomspace(fl.v_PT(fl.pmax, fl.Tmin) * 1.02,
                                               fl.vmax / 2, 7), fl.vc))

        a = self.ax_pv
        for T in Tbg:
            v, _, p = isotherm(fl, T)
            a.plot(v, dp(p), **g)
        a.plot(self.dome_v, dp(self.dome_p), **dome)
        a.set(xscale="log", yscale="log", xlabel=self._lab("v"), ylabel=self._lab("p"))

        a = self.ax_Tv
        for p in pbg:
            v, T, _ = isobar(fl, p)
            a.plot(v, dT(T), **g)
        a.plot(self.dome_v, dT(self.dome_T), **dome)
        a.set(xscale="log", xlabel=self._lab("v"), ylabel=self._lab("T"))

        a = self.ax_pT
        for v in vbg:
            _, T, p = isochore(fl, v)
            a.plot(dT(T), dp(p), **g)
        a.plot(dT(fl.sat_T), dp(fl.sat_p), **dome)
        a.set(yscale="log", xlabel=self._lab("T"), ylabel=self._lab("p"))

        self.panel_names = {self.ax_pv: ("p–v diagram", "gray lines: isotherms"),
                            self.ax_Tv: ("T–v diagram", "gray lines: isobars"),
                            self.ax_pT: ("p–T diagram", "gray lines: isochores")}
        self.role_titles = {}
        for a, (name, sub) in self.panel_names.items():
            a.set_title(f"{name}", loc="left", pad=6)
            a.text(0.5, 1.0, sub, transform=a.transAxes, ha="center", va="bottom",
                   fontsize=8, color=MUTED)
            self.role_titles[a] = a.set_title("", loc="right", pad=6, fontsize=8.5)

        # region labels in phase colors. Positions are computed from the
        # saturation data (not fixed screen fractions) so every label is
        # guaranteed to sit inside the region it names, for any fluid.
        lk = dict(fontsize=9, style="italic", weight="bold", ha="center",
                  va="center", zorder=4)
        xl, yl, zl = self.lims
        p_floor = fl.p_si(10 ** zl[0])            # bottom of the p axes (SI)

        def T_on_vg(v):   # temperature where the sat. vapor line has volume v
            return np.interp(np.log(v), np.log(fl.sat_vg[::-1]), fl.sat_T[::-1])

        # liquid + vapor: middle of the dome
        Tm = fl.Tmin + 0.6 * (fl.Tc - fl.Tmin)
        vm = np.sqrt(fl.vf(Tm) * fl.vg(Tm))
        self.ax_pv.text(vm, dp(np.sqrt(fl.psat(fl.Tmin) * fl.psat(Tm))), "liquid + vapor",
                        color=PHASE_TEXT["mix"], **lk)
        self.ax_Tv.text(vm, dT(fl.Tmin + 0.3 * (Tm - fl.Tmin)), "liquid + vapor",
                        color=PHASE_TEXT["mix"], **lk)

        # vapor: right of the saturated-vapor line, below the critical isotherm
        v_lab = 10 ** (xl[0] + 0.8 * (xl[1] - xl[0]))
        if v_lab < fl.vg(fl.Tmin):
            T_star = T_on_vg(v_lab)
        else:                                     # dome doesn't reach this far
            T_star = fl.Tmin
        # In p-v the band between the dome and the critical isotherm is thin at
        # large v, so put the label above the hottest isotherm instead: that is
        # still vapor/gas (T > T_c but p << p_c) and there is open space there.
        p_hot = fl.p_Tv(fl.Tmax, v_lab) * 10 ** (0.12 * (zl[1] - zl[0]))
        self.ax_pv.text(v_lab, dp(min(p_hot, 0.5 * fl.pc)), "vapor",
                        color=PHASE_TEXT["vapor"], **lk)
        self.ax_Tv.text(v_lab, dT(T_star + 0.45 * (fl.Tc - T_star)), "vapor",
                        color=PHASE_TEXT["vapor"], **lk)

        # p-T: liquid above the vaporization curve, vapor below it,
        # supercritical above and to the right of the critical point
        TL = fl.Tmin + 0.25 * (fl.Tc - fl.Tmin)
        self.ax_pT.text(dT(TL), dp(np.sqrt(fl.psat(TL) * fl.pmax)), "liquid",
                        color=PHASE_TEXT["liquid"], **lk)
        TV = fl.Tmin + 0.6 * (fl.Tc - fl.Tmin)
        self.ax_pT.text(dT(TV), dp(np.sqrt(p_floor * fl.psat(TV))), "vapor / gas",
                        color=PHASE_TEXT["vapor"], **lk)
        self.ax_pT.text(dT(fl.Tc + 0.55 * (fl.Tmax - fl.Tc)), dp(np.sqrt(fl.pc * fl.pmax)),
                        "supercritical", color=PHASE_TEXT["super"], **lk)

        for a, x, y in ((self.ax_pv, fl.vc, dp(fl.pc)), (self.ax_Tv, fl.vc, dT(fl.Tc)),
                        (self.ax_pT, dT(fl.Tc), dp(fl.pc))):
            a.plot(x, y, "o", ms=6, mfc="white", mec=INK, mew=1.5, zorder=6)
            a.annotate("critical point", (x, y), xytext=(7, 5),
                       textcoords="offset points", fontsize=8, color=INK, zorder=6)
            a.grid(True, which="major")
            a.set_axisbelow(True)
            a.autoscale(False)

        xl, yl, zl = self.lims
        self.ax_pv.set_xlim(10 ** xl[0], 10 ** xl[1]); self.ax_pv.set_ylim(10 ** zl[0], 10 ** zl[1])
        self.ax_Tv.set_xlim(10 ** xl[0], 10 ** xl[1]); self.ax_Tv.set_ylim(*yl)
        self.ax_pT.set_xlim(*yl); self.ax_pT.set_ylim(10 ** zl[0], 10 ** zl[1])

    # ------------------------------------------------------------------
    def _on_mode(self, label):
        fl = self.fl
        self.mode = label[0]
        d = fl.defaults
        if self.mode == "p":
            lo, hi = np.log10(fl.p_disp(fl.psat(fl.Tmin) * 1.5)), np.log10(fl.p_disp(fl.pmax / 1.2))
            val, title = np.log10(fl.p_disp(d["p"])), "PRESSURE"
            rng = (f"{10**lo:.3g} – {10**hi:.3g} {fl.p_unit}".strip(), "log scale")
        elif self.mode == "v":
            lo = np.log10(fl.v_PT(fl.pmax, fl.Tmin) * 1.02)
            hi = np.log10(fl.vmax / 1.5)
            val, title = np.log10(d["v"]), "SPECIFIC VOLUME"
            rng = (f"{10**lo:.3g} – {10**hi:.3g} {fl.v_unit}".strip(), "log scale")
        else:
            lo, hi = float(fl.T_disp(fl.Tmin)), float(fl.T_disp(fl.Tmax))
            val, title = float(fl.T_disp(d["T"])), "TEMPERATURE"
            rng = (f"{lo:.3g} – {hi:.3g} {fl.T_unit}".strip(), "")
        self.slider_title.set_text(title)
        self.range_text.set_text(rng[0] + (f"  ({rng[1]})" if rng[1] else ""))
        s = self.slider
        s.valmin, s.valmax = lo, hi
        s.ax.set_xlim(lo, hi)
        if hasattr(s.poly, "set_x"):
            s.poly.set_x(lo)
        else:
            s.poly.set_visible(False)
        s.set_val(np.clip(val, lo, hi))

    def _on_check(self, label):
        on = self.checks.get_status()
        for a in self.wall_static + self.wall_dyn:
            a.set_visible(on[0])
        for a in self.surf:
            a.set_visible(on[1])
        self.fig.canvas.draw_idle()

    # ------------------------------------------------------------------
    def _update(self, val):
        fl, dT, dp = self.fl, self.fl.T_disp, self.fl.p_disp
        ax = self.ax3
        for a in self.dyn + self.wall_dyn:
            a.remove()
        self.dyn, self.wall_dyn = [], []
        xl, yl, zl = self.lims

        if self.mode == "p":
            p = float(fl.p_si(10 ** val))
            v, T, P = isobar(fl, p)
            txt = f"{fl.p_label} = {dp(p):.4g} {fl.p_unit}"
            msg = self._describe_p(p)
            plane_axis, plane_val, home = "z", np.log10(dp(p)), self.ax_Tv
        elif self.mode == "v":
            vv = 10 ** val
            v, T, P = isochore(fl, vv)
            txt = f"{fl.v_label} = {vv:.4g} {fl.v_unit}"
            msg = self._describe_v(vv)
            plane_axis, plane_val, home = "x", val, self.ax_pT
        else:
            Tk = float(fl.T_si(val))
            v, T, P = isotherm(fl, Tk)
            txt = f"{fl.T_label} = {val:.4g} {fl.T_unit}"
            msg = self._describe_T(Tk)
            plane_axis, plane_val, home = "y", val, self.ax_pv
        self.value_text.set_text(txt)
        paras = [self.wrap.fill(s, 82) for s in msg.split("\n")]
        self.status.set_text("\n".join(paras))

        # --- 3-D: cutting plane with outline + intersection curve
        c = np.array([0.0, 1.0, 1.0, 0.0, 0.0])
        r = np.array([0.0, 0.0, 1.0, 1.0, 0.0])
        A, B = np.meshgrid([0.0, 1.0], [0.0, 1.0])
        sx, sy, sz = (lambda t: xl[0] + t * np.ptp(xl)), (lambda t: yl[0] + t * np.ptp(yl)), \
                     (lambda t: zl[0] + t * np.ptp(zl))
        if plane_axis == "x":
            surf = (np.full_like(A, plane_val), sy(A), sz(B))
            edge = (np.full_like(c, plane_val), sy(c), sz(r))
        elif plane_axis == "y":
            surf = (sx(A), np.full_like(A, plane_val), sz(B))
            edge = (sx(c), np.full_like(c, plane_val), sz(r))
        else:
            surf = (sx(A), sy(B), np.full_like(A, plane_val))
            edge = (sx(c), sy(r), np.full_like(c, plane_val))
        self.dyn.append(ax.plot_surface(*surf, color=ACCENT, alpha=0.09, shade=False,
                                        linewidth=0, zorder=3))
        self.dyn.append(ax.plot(*edge, color=ACCENT, lw=0.8, alpha=0.45, zorder=3)[0])

        lv, lT, lp = np.log10(v), dT(T), np.log10(dp(P))
        self.dyn.append(ax.plot(lv, lT, lp, color="white", lw=5.5, alpha=0.9,
                                solid_capstyle="round", zorder=7)[0])
        self.dyn.append(ax.plot(lv, lT, lp, color=ACCENT, lw=3, solid_capstyle="round",
                                zorder=8)[0])

        kw = dict(color=ACCENT, lw=1.3, ls=(0, (4, 2)), alpha=0.8, zorder=0.5)
        self.wall_dyn = [
            ax.plot(lv, np.full_like(lv, yl[1]), lp, **kw)[0],
            ax.plot(lv, lT, np.full_like(lv, zl[0]), **kw)[0],
            ax.plot(np.full_like(lv, xl[0]), lT, lp, **kw)[0],
        ]
        for a in self.wall_dyn:
            a.set_visible(self.checks.get_status()[0])

        # --- 2-D: same curve in every plane; tag the slice plane
        ck = dict(color=ACCENT, lw=2.6, zorder=5, solid_capstyle="round")
        self.dyn.append(self.ax_pv.plot(v, dp(P), **ck)[0])
        self.dyn.append(self.ax_Tv.plot(v, dT(T), **ck)[0])
        self.dyn.append(self.ax_pT.plot(dT(T), dp(P), **ck)[0])
        for a in (self.ax_pv, self.ax_Tv, self.ax_pT):
            is_home = a is home
            a.set_facecolor(ACCENT_TINT if is_home else CARD)
            for sp in a.spines.values():
                sp.set_color(ACCENT if is_home else EDGE)
                sp.set_linewidth(1.6 if is_home else 0.9)
            t = self.role_titles[a]
            t.set_text("● SLICE PLANE" if is_home else "projection")
            t.set_color(ACCENT if is_home else MUTED)
            t.set_fontweight("bold" if is_home else "normal")
        self.fig.canvas.draw_idle()

    # ---- plain-language explanations for the status card ----------------
    def _f(self, q, x):
        fl = self.fl
        if q == "T":
            return f"{float(fl.T_disp(x)):.4g} {fl.T_unit}".strip()
        if q == "p":
            return f"{float(fl.p_disp(x)):.4g} {fl.p_unit}".strip()
        return f"{x:.4g} {fl.v_unit}".strip()

    def _describe_p(self, p):
        fl = self.fl
        if abs(p / fl.pc - 1) < 0.02:
            return ("Critical isobar: the boiling plateau has shrunk to a single "
                    "inflection point at the critical point.")
        if p > fl.pc:
            return (f"Supercritical isobar (p > p_c = {self._f('p', fl.pc)}): heating turns "
                    "the liquid-like fluid into a gas-like fluid continuously. There is "
                    "no boiling and no plateau in the T–v diagram.")
        Ts = float(fl.Tsat(p))
        return (f"Constant-pressure heating below p_c: compressed liquid → saturated liquid "
                f"(v_f = {self._f('v', fl.vf(Ts))}) → boils at constant T_sat = {self._f('T', Ts)} "
                f"→ saturated vapor (v_g = {self._f('v', fl.vg(Ts))}) → superheated vapor.\n"
                "The flat segment in T–v is the boiling process; in p–T the whole thing "
                "collapses to a single point on the vaporization curve.")

    def _describe_T(self, T):
        fl = self.fl
        if abs(T - fl.Tc) < 0.005 * fl.Tc:
            return ("Critical isotherm: it has a horizontal inflection point at the "
                    "critical point (∂p/∂v = ∂²p/∂v² = 0).")
        if T > fl.Tc:
            return (f"Isotherm above T_c = {self._f('T', fl.Tc)}: no phase change at any "
                    "pressure. At large v it approaches ideal-gas behavior, pv ≈ RT.")
        ps = fl.psat(T)
        return (f"Isotherm below T_c: compressed liquid (steep, nearly incompressible) → "
                f"saturated liquid → liquid-vapor mixture at constant p_sat = {self._f('p', ps)} "
                f"(v from {self._f('v', fl.vf(T))} to {self._f('v', fl.vg(T))}) → superheated vapor.\n"
                "Inside the dome, constant T means constant p, so in p–T the mixture is one point.")

    def _describe_v(self, v):
        fl = self.fl
        if abs(v / fl.vc - 1) < 0.03:
            return ("Rigid tank filled at the critical specific volume: on heating, the "
                    "liquid-vapor interface stays put and vanishes at the critical point.")
        Tx = dome_crossing_T(fl, v)
        if Tx is None:
            return "This isochore does not enter the vapor dome in the plotted range."
        px = fl.psat(Tx)
        if v < fl.vc:
            return (f"Rigid tank, v < v_c: heating a liquid-vapor mixture raises p along the "
                    f"saturation curve until the tank is full of saturated liquid at "
                    f"T = {self._f('T', Tx)}, p = {self._f('p', px)}.\n"
                    "After that, p shoots up (compressed liquid). Note the kink in p–T where "
                    "the isochore leaves the vaporization curve.")
        return (f"Rigid tank, v > v_c: heating a liquid-vapor mixture evaporates the liquid "
                f"until it is saturated vapor at T = {self._f('T', Tx)}, p = {self._f('p', px)}.\n"
                "After that it is superheated vapor and p rises gently (roughly p ∝ T, like an "
                "ideal gas).")


# =============================================================================
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--eos", choices=["coolprop", "vdw"], default="coolprop")
    ap.add_argument("--fluid", default="Water", help="CoolProp fluid name (default Water)")
    ap.add_argument("--save", metavar="PREFIX", help="save PNGs of the three modes and exit")
    args = ap.parse_args()

    if args.save:
        import matplotlib
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fl = None
    if args.eos == "coolprop":
        try:
            fl = CoolPropFluid(args.fluid)
        except ImportError:
            print("CoolProp not found (pip install CoolProp); using a van der Waals fluid instead.")
    if fl is None:
        fl = VanDerWaals()

    ex = Explorer(fl)
    if args.save:
        for i, lab in enumerate(MODES):
            ex.radio.set_active(i)
            ex.fig.savefig(f"{args.save}_{lab[0]}.png", dpi=120)
            print("saved", f"{args.save}_{lab[0]}.png")
        return
    plt.show()


if __name__ == "__main__":
    main()