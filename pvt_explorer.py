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
    T_label, p_label, v_label = "T_r", "p_r", "v_r"

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

    T_label, p_label, v_label = "T", "p", "v"
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
#  The 3-D surface
# =============================================================================
PHASE_COLORS = {
    "Compressed liquid": "#3d7fd1",
    "Liquid + vapor": "#6cbf6a",
    "Superheated vapor / gas": "#f08c3a",
    "Supercritical fluid": "#a9a9b8",
}


def build_surface(fl, nT=55, nL=20, n2=20, nV=36):
    """Grid whose columns follow the phase boundaries, so the dome edges are
    crisp: each isotherm row = [liquid | two-phase | vapor] segments."""
    Ts = np.unique(np.append(np.linspace(fl.Tmin, fl.Tmax, nT), fl.Tc))
    from matplotlib.colors import to_rgba
    colors = {k: np.array(to_rgba(c, 0.92)) for k, c in PHASE_COLORS.items()}
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
        # color of the cell to the right of each node
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
                key = "Supercritical fluid" if pm > fl.pc else "Superheated vapor / gas"
            c[j] = colors[key]
        V.append(v)
        P.append(p)
        C.append(c)
    V, P, C = np.array(V), np.array(P), np.array(C)
    return V, np.repeat(Ts[:, None], V.shape[1], axis=1), P, C


# =============================================================================
#  The interactive figure
# =============================================================================
class Explorer:
    def __init__(self, fl):
        import matplotlib.pyplot as plt
        from matplotlib.patches import Patch
        from matplotlib.ticker import FuncFormatter
        from matplotlib.widgets import CheckButtons, RadioButtons, Slider

        self.plt, self.fl = plt, fl
        self.mode = "T"
        self.dyn = []           # artists that change with the slider
        self.wall_dyn = []      # red-curve wall projections

        # ---------- data --------------------------------------------------
        print("Building surface ...", flush=True)
        V, T, P, C = build_surface(fl)
        self.X = np.log10(V)
        self.Y = fl.T_disp(T)
        self.Z = np.log10(fl.p_disp(P))
        dT, dp = fl.T_disp, fl.p_disp

        # dome (sat. liquid line up to critical point, sat. vapor line down)
        self.dome_v = np.concatenate([fl.sat_vf, fl.sat_vg[::-1]])
        self.dome_T = np.concatenate([fl.sat_T, fl.sat_T[::-1]])
        self.dome_p = np.concatenate([fl.sat_p, fl.sat_p[::-1]])

        # ---------- figure -----------------------------------------------
        self.fig = fig = plt.figure(figsize=(15.5, 9))
        fig.canvas.manager.set_window_title("p-v-T surface explorer")
        fig.suptitle(f"p–v–T surface: {fl.name}", fontsize=14, weight="bold",
                     x=0.27, y=0.985)

        ax = self.ax3 = fig.add_axes([0.0, 0.2, 0.55, 0.77], projection="3d")
        self.surf = ax.plot_surface(self.X, self.Y, self.Z, facecolors=C,
                                    rstride=1, cstride=1, linewidth=0.15,
                                    edgecolor=(0, 0, 0, 0.25), shade=False,
                                    antialiased=False)
        ax.plot(np.log10(self.dome_v), dT(self.dome_T),
                np.log10(dp(self.dome_p)), "k-", lw=2.2, zorder=10)
        self.cp3 = (np.log10(fl.vc), float(dT(fl.Tc)), np.log10(float(dp(fl.pc))))
        ax.scatter(*self.cp3, color="k", s=40, zorder=11, depthshade=False)

        xl = (np.nanmin(self.X) - 0.35, np.nanmax(self.X) + 0.05)
        yl = (np.nanmin(self.Y), np.nanmax(self.Y))
        zl = (np.nanmin(self.Z) - 0.2, np.nanmax(self.Z) + 0.05)
        self.lims = xl, yl, zl
        ax.set_xlim(xl); ax.set_ylim(yl); ax.set_zlim(zl)
        pow10 = FuncFormatter(lambda x, _: f"$10^{{{int(round(x))}}}$")
        ax.set_xticks(np.arange(np.ceil(xl[0]), np.floor(xl[1]) + 1))
        ax.set_zticks(np.arange(np.ceil(zl[0]), np.floor(zl[1]) + 1))
        ax.xaxis.set_major_formatter(pow10)
        ax.zaxis.set_major_formatter(pow10)
        ax.set_xlabel(self._lab("v"), labelpad=8)
        ax.set_ylabel(self._lab("T"), labelpad=8)
        ax.set_zlabel(self._lab("p"), labelpad=6)
        ax.view_init(elev=22, azim=-58)
        ax.legend(handles=[Patch(color=c, label=k) for k, c in PHASE_COLORS.items()],
                  loc="upper left", fontsize=8.5, framealpha=0.9)

        # static wall projections of the dome (far walls for default view)
        xw, yw, zw = xl[0], yl[1], zl[0]
        lv, lT, lp = np.log10(self.dome_v), dT(self.dome_T), np.log10(dp(self.dome_p))
        kw = dict(color="0.35", lw=1.3)
        self.wall_static = [
            ax.plot(lv, np.full_like(lv, yw), lp, **kw)[0],        # p-v wall
            ax.plot(lv, lT, np.full_like(lv, zw), **kw)[0],        # T-v wall
            ax.plot(np.full_like(lT, xw), lT, lp, **kw)[0],        # p-T wall
        ]

        # ---------- 2-D panels -------------------------------------------
        w, h, x0 = 0.37, 0.215, 0.61
        self.ax_pv = fig.add_axes([x0, 0.745, w, h])
        self.ax_Tv = fig.add_axes([x0, 0.445, w, h])
        self.ax_pT = fig.add_axes([x0, 0.145, w, h])
        self._setup_2d()

        # ---------- widgets ----------------------------------------------
        fig.text(0.02, 0.165, "Hold constant:", fontsize=10, weight="bold")
        rax = fig.add_axes([0.02, 0.035, 0.1, 0.12])
        self.radio = RadioButtons(rax, ["p (isobar)", "v (isochore)", "T (isotherm)"],
                                  active=2)
        cax = fig.add_axes([0.13, 0.06, 0.13, 0.09])
        self.checks = CheckButtons(cax, ["Wall projections", "Show surface"],
                                   [True, True])
        sax = fig.add_axes([0.31, 0.12, 0.22, 0.03])
        self.slider = Slider(sax, "", 0, 1, valinit=0.5, color="tab:red")
        self.slider.vline.set_visible(False)
        self.slider_title = fig.text(0.31, 0.16, "", fontsize=10, weight="bold")
        self.status = fig.text(0.29, 0.005, "", fontsize=9, va="bottom",
                               wrap=True, linespacing=1.35)

        self.radio.on_clicked(self._on_mode)
        self.checks.on_clicked(self._on_check)
        self.slider.on_changed(self._update)
        self._on_mode("T (isotherm)")

    # ------------------------------------------------------------------
    def _lab(self, q):
        fl = self.fl
        name = {"p": fl.p_label, "v": fl.v_label, "T": fl.T_label}[q]
        unit = {"p": fl.p_unit, "v": fl.v_unit, "T": fl.T_unit}[q]
        return f"{name} [{unit}]" if unit else name

    def _setup_2d(self):
        fl, dT, dp = self.fl, self.fl.T_disp, self.fl.p_disp
        g = dict(color="0.72", lw=0.9)
        dome = dict(color="k", lw=2)

        # background families of curves
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
        a.plot(self.dome_v, dp(self.dome_p), **dome, label="Saturation (vapor dome)")
        a.set(xscale="log", yscale="log", xlabel=self._lab("v"), ylabel=self._lab("p"))
        a.set_title("p–v plane   (gray: isotherms)", fontsize=10)

        a = self.ax_Tv
        for p in pbg:
            v, T, _ = isobar(fl, p)
            a.plot(v, dT(T), **g)
        a.plot(self.dome_v, dT(self.dome_T), **dome)
        a.set(xscale="log", xlabel=self._lab("v"), ylabel=self._lab("T"))
        a.set_title("T–v plane   (gray: isobars)", fontsize=10)

        a = self.ax_pT
        for v in vbg:
            _, T, p = isochore(fl, v)
            a.plot(dT(T), dp(p), **g)
        a.plot(dT(fl.sat_T), dp(fl.sat_p), **dome)
        a.set(yscale="log", xlabel=self._lab("T"), ylabel=self._lab("p"))
        a.set_title("p–T plane   (gray: isochores)", fontsize=10)
        tk = dict(transform=a.transAxes, fontsize=9, color="0.25", style="italic")
        a.text(0.08, 0.82, "Liquid", **tk)
        a.text(0.62, 0.15, "Vapor / gas", **tk)
        a.text(0.72, 0.9, "Supercritical", **tk)
        a.text(0.33, 0.47, "vaporization\ncurve", **tk)

        # region labels and critical points
        Tm = fl.Tmin + 0.6 * (fl.Tc - fl.Tmin)
        vm = np.sqrt(fl.vf(Tm) * fl.vg(Tm))
        self.ax_pv.text(vm, dp(fl.psat(Tm)) * 0.55, "L + V", ha="center",
                        fontsize=9, style="italic", color="0.25")
        self.ax_Tv.text(vm, dT(Tm - 0.12 * (fl.Tc - fl.Tmin)), "L + V",
                        ha="center", fontsize=9, style="italic", color="0.25")
        for a, x, y in ((self.ax_pv, fl.vc, dp(fl.pc)), (self.ax_Tv, fl.vc, dT(fl.Tc)),
                        (self.ax_pT, dT(fl.Tc), dp(fl.pc))):
            a.plot(x, y, "ko", ms=5)
            a.annotate("critical point", (x, y), xytext=(6, 6),
                       textcoords="offset points", fontsize=8)
            a.grid(True, which="major", alpha=0.3)
            a.tick_params(labelsize=8)
            a.autoscale(False)

        # sensible limits
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
            val, title = np.log10(fl.p_disp(d["p"])), "Pressure (log scale)"
        elif self.mode == "v":
            lo = np.log10(fl.v_PT(fl.pmax, fl.Tmin) * 1.02)
            hi = np.log10(fl.vmax / 1.5)
            val, title = np.log10(d["v"]), "Specific volume (log scale)"
        else:
            lo, hi = float(fl.T_disp(fl.Tmin)), float(fl.T_disp(fl.Tmax))
            val, title = float(fl.T_disp(d["T"])), "Temperature"
        self.slider_title.set_text(title)
        s = self.slider
        s.valmin, s.valmax = lo, hi
        s.ax.set_xlim(lo, hi)
        if hasattr(s.poly, "set_x"):      # keep the filled bar anchored at the left end
            s.poly.set_x(lo)
        else:
            s.poly.set_visible(False)
        s.set_val(np.clip(val, lo, hi))

    def _on_check(self, label):
        on = self.checks.get_status()
        for a in self.wall_static + self.wall_dyn:
            a.set_visible(on[0])
        self.surf.set_visible(on[1])
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
            txt = f"p = {dp(p):.4g} {fl.p_unit}"
            self.status.set_text(self._describe_p(p))
            plane_axis, plane_val, home = "z", np.log10(dp(p)), self.ax_Tv
        elif self.mode == "v":
            vv = 10 ** val
            v, T, P = isochore(fl, vv)
            txt = f"v = {vv:.4g} {fl.v_unit}"
            self.status.set_text(self._describe_v(vv))
            plane_axis, plane_val, home = "x", val, self.ax_pT
        else:
            Tk = float(fl.T_si(val))
            v, T, P = isotherm(fl, Tk)
            txt = f"T = {val:.4g} {fl.T_unit}"
            self.status.set_text(self._describe_T(Tk))
            plane_axis, plane_val, home = "y", val, self.ax_pv
        self.slider.valtext.set_text(txt)

        # --- 3-D: cutting plane + intersection curve
        g = np.linspace(0, 1, 2)
        A, B = np.meshgrid(g, g)
        if plane_axis == "x":
            PX, PY, PZ = np.full_like(A, plane_val), yl[0] + A * np.ptp(yl), zl[0] + B * np.ptp(zl)
        elif plane_axis == "y":
            PX, PY, PZ = xl[0] + A * np.ptp(xl), np.full_like(A, plane_val), zl[0] + B * np.ptp(zl)
        else:
            PX, PY, PZ = xl[0] + A * np.ptp(xl), yl[0] + B * np.ptp(yl), np.full_like(A, plane_val)
        self.dyn.append(ax.plot_surface(PX, PY, PZ, color="tab:red", alpha=0.13,
                                        shade=False, zorder=1))
        lv, lT, lp = np.log10(v), dT(T), np.log10(dp(P))
        self.dyn.append(ax.plot(lv, lT, lp, color="red", lw=3, zorder=20)[0])

        kw = dict(color="red", lw=1.6, ls="--")
        self.wall_dyn = [
            ax.plot(lv, np.full_like(lv, yl[1]), lp, **kw)[0],
            ax.plot(lv, lT, np.full_like(lv, zl[0]), **kw)[0],
            ax.plot(np.full_like(lv, xl[0]), lT, lp, **kw)[0],
        ]
        for a in self.wall_dyn:
            a.set_visible(self.checks.get_status()[0])

        # --- 2-D: same curve in every plane
        ck = dict(color="red", lw=2.6, zorder=5)
        self.dyn.append(self.ax_pv.plot(v, dp(P), **ck)[0])
        self.dyn.append(self.ax_Tv.plot(v, dT(T), **ck)[0])
        self.dyn.append(self.ax_pT.plot(dT(T), dp(P), **ck)[0])
        for a in (self.ax_pv, self.ax_Tv, self.ax_pT):
            for sp in a.spines.values():
                sp.set_color("red" if a is home else "black")
                sp.set_linewidth(2.2 if a is home else 0.8)
        self.fig.canvas.draw_idle()

    # ---- plain-language explanations for the status line ----------------
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
                    "the liquid-like fluid into a gas-like fluid continuously; "
                    "no boiling, no plateau in T–v.")
        Ts = float(fl.Tsat(p))
        return (f"Isobar below p_c: compressed liquid → saturated liquid (v_f = {self._f('v', fl.vf(Ts))}) "
                f"→ boils at constant T_sat = {self._f('T', Ts)} (flat line in T–v) → saturated vapor "
                f"(v_g = {self._f('v', fl.vg(Ts))}) → superheated vapor.\n"
                "In p–T the whole boiling process is a single point on the vaporization curve.")

    def _describe_T(self, T):
        fl = self.fl
        if abs(T - fl.Tc) < 0.005 * fl.Tc:
            return ("Critical isotherm: horizontal inflection at the critical point "
                    "(∂p/∂v = ∂²p/∂v² = 0).")
        if T > fl.Tc:
            return (f"Isotherm above T_c = {self._f('T', fl.Tc)}: no phase change at any pressure; "
                    "at large v it approaches ideal-gas behaviour, pv ≈ RT.")
        ps = fl.psat(T)
        return (f"Isotherm below T_c: compressed liquid (steep, nearly incompressible) → saturated "
                f"liquid → liquid-vapor mixture at constant p_sat = {self._f('p', ps)} "
                f"(v from {self._f('v', fl.vf(T))} to {self._f('v', fl.vg(T))}) → superheated vapor.\n"
                "Constant T and constant p go together inside the dome: in p–T the mixture is one point.")

    def _describe_v(self, v):
        fl = self.fl
        if abs(v / fl.vc - 1) < 0.03:
            return ("Rigid tank filled at the critical specific volume: on heating, the meniscus "
                    "stays put and vanishes at the critical point.")
        Tx = dome_crossing_T(fl, v)
        if Tx is None:
            return "This isochore does not enter the vapor dome in the plotted range."
        px = fl.psat(Tx)
        if v < fl.vc:
            return (f"Rigid tank, v < v_c: heating a liquid-vapor mixture raises p along the saturation "
                    f"curve until the tank is full of saturated liquid at T = {self._f('T', Tx)}, "
                    f"p = {self._f('p', px)};\nafter that p shoots up (compressed liquid). "
                    "Note the kink in p–T where the isochore leaves the vaporization curve.")
        return (f"Rigid tank, v > v_c: heating a liquid-vapor mixture evaporates the liquid; it becomes "
                f"saturated vapor at T = {self._f('T', Tx)}, p = {self._f('p', px)};\n"
                "after that it is superheated vapor and p rises gently (roughly p ∝ T, like an ideal gas).")


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
        for lab in ("p (isobar)", "v (isochore)", "T (isotherm)"):
            ex.radio.set_active(["p (isobar)", "v (isochore)", "T (isotherm)"].index(lab))
            ex.fig.savefig(f"{args.save}_{lab[0]}.png", dpi=110)
            print("saved", f"{args.save}_{lab[0]}.png")
        return
    plt.show()


if __name__ == "__main__":
    main()
