"""
Drug Dosing Schedule Model — Streamlit App
Kyle Cheng | HDR Framework (Evolutionary Game Theory)
"""

import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

# ─────────────────────────────────────────────
# Page Config
# ─────────────────────────────────────────────
st.set_page_config(
    page_title="Drug Dosing Schedule Model",
    layout="wide",
)

# ─────────────────────────────────────────────
# Custom CSS
# ─────────────────────────────────────────────
st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;600&family=IBM+Plex+Sans:wght@300;400;600&display=swap');

html, body, [class*="css"] {
    font-family: 'IBM Plex Sans', sans-serif;
}

h1, h2, h3 {
    font-family: 'IBM Plex Mono', monospace !important;
    letter-spacing: -0.02em;
}

.block-container {
    padding-top: 2rem;
}

.stSidebar {
    background-color: #0d0d0d !important;
}
.stSidebar * {
    color: #e0e0e0 !important;
}

.metric-card {
    background: #f7f7f2;
    border-left: 4px solid #1a1a2e;
    padding: 0.75rem 1rem;
    border-radius: 4px;
    margin-bottom: 0.5rem;
    font-family: 'IBM Plex Mono', monospace;
    font-size: 0.85rem;
}

.metric-card.positive { border-left-color: #e53935; }
.metric-card.negative { border-left-color: #1976d2; }
.metric-card.neutral  { border-left-color: #43a047; }

.elim-badge {
    display: inline-block;
    padding: 0.3rem 0.8rem;
    border-radius: 20px;
    font-family: 'IBM Plex Mono', monospace;
    font-size: 0.8rem;
    font-weight: 600;
}
.elim-yes { background: #e8f5e9; color: #2e7d32; border: 1px solid #81c784; }
.elim-no  { background: #fff3e0; color: #e65100; border: 1px solid #ffb74d; }

.section-header {
    font-family: 'IBM Plex Mono', monospace;
    font-size: 0.7rem;
    letter-spacing: 0.12em;
    text-transform: uppercase;
    color: #888;
    margin-bottom: 0.25rem;
    margin-top: 1.5rem;
}
</style>
""", unsafe_allow_html=True)

# ─────────────────────────────────────────────
# Header
# ─────────────────────────────────────────────
st.markdown("# Drug Dosing Schedule Model")
st.markdown("**Based on Hawk Dove Retaliator Game** — Sensitive · Resistant · Drug-Tolerant Persisters")
st.divider()

# ─────────────────────────────────────────────
# Sidebar Inputs
# ─────────────────────────────────────────────
with st.sidebar:
    st.markdown("## ⚙️ Parameters")

    st.markdown('<div class="section-header">Synergy</div>', unsafe_allow_html=True)
    syn_score = st.number_input("Synergy Score", value=0.5, step=0.05, format="%.3f",
                                help="Drug synergy score from ML model")

    st.markdown('<div class="section-header">Sensitive Cells (S / Dove)</div>', unsafe_allow_html=True)
    s0       = st.number_input("Initial Count s₀", value=700.0, min_value=1.0, step=50.0)
    lambda_s = st.number_input("Growth rate under Drug 1 (λ_S)", value=-0.3, step=0.05, format="%.3f")
    mu_s     = st.number_input("Growth rate under Drug 2 (μ_S)", value=0.2,  step=0.05, format="%.3f")

    st.markdown('<div class="section-header">Resistant Cells (R / Hawk)</div>', unsafe_allow_html=True)
    r0       = st.number_input("Initial Count r₀", value=200.0, min_value=0.0, step=50.0)
    lambda_r = st.number_input("Growth rate under Drug 1 (λ_R)", value=0.1,  step=0.05, format="%.3f")
    mu_r     = st.number_input("Growth rate under Drug 2 (μ_R)", value=0.05, step=0.05, format="%.3f")

    st.markdown('<div class="section-header">DTP Cells (P / Retaliator)</div>', unsafe_allow_html=True)
    p0_count = st.number_input("Initial Count p₀", value=100.0, min_value=0.0, step=50.0)
    lambda_p = st.number_input("Growth rate under Drug 1 (λ_P)", value=-0.01, step=0.005, format="%.4f",
                               help="Near 0 — dormancy under drug")
    mu_p     = st.number_input("Growth rate under Drug 2 (μ_P)", value=-0.01, step=0.005, format="%.4f",
                               help="Near 0 — dormancy under drug")

    st.markdown('<div class="section-header">HDR Interaction Matrix</div>', unsafe_allow_html=True)
    with st.expander("Expand to edit payoffs"):
        a_SS = st.number_input("S vs S (Dove-Dove)",          value=0.0,  step=0.05, format="%.2f")
        a_SR = st.number_input("S vs R (Dove-Hawk suppressed)",value=-0.3, step=0.05, format="%.2f")
        a_SP = st.number_input("S vs P (Dove-Ret, neutral)",   value=0.0,  step=0.05, format="%.2f")
        a_RS = st.number_input("R vs S (Hawk exploits Dove)",  value=0.3,  step=0.05, format="%.2f")
        a_RR = st.number_input("R vs R (Hawk-Hawk cost)",      value=-0.4, step=0.05, format="%.2f")
        a_RP = st.number_input("R vs P (Ret retaliates)",      value=-0.2, step=0.05, format="%.2f")
        a_PS = st.number_input("P vs S (Ret-Dove neutral)",    value=0.0,  step=0.05, format="%.2f")
        a_PR = st.number_input("P vs R (dormancy survival)",   value=-0.1, step=0.05, format="%.2f")
        a_PP = st.number_input("P vs P (cooperative)",         value=0.0,  step=0.05, format="%.2f")

    st.divider()
    run_btn = st.button("▶ Run Simulation", type="primary", use_container_width=True)

# ─────────────────────────────────────────────
# Model Core 
# ─────────────────────────────────────────────
SYNERGY_SCALE = 0.015

A = np.array([
    [a_SS, a_SR, a_SP],
    [a_RS, a_RR, a_RP],
    [a_PS, a_PR, a_PP],
])

def base_fitness(X, lam, mu, syn):
    penalty = SYNERGY_SCALE * syn * X * (1 - X)
    return X * lam + (1 - X) * mu - penalty

def fitness_all(p_S, p_R, X, syn):
    p_P = max(0.0, 1.0 - p_S - p_R)
    freqs = np.array([p_S, p_R, p_P])
    W_S = base_fitness(X, lambda_s, mu_s, syn) + A[0] @ freqs
    W_R = base_fitness(X, lambda_r, mu_r, syn) + A[1] @ freqs
    W_P = base_fitness(X, lambda_p, mu_p, syn) + A[2] @ freqs
    return W_S, W_R, W_P

def w_bar(p_S, p_R, W_S, W_R, W_P):
    p_P = max(0.0, 1.0 - p_S - p_R)
    return p_S * W_S + p_R * W_R + p_P * W_P

def delta_p(p_S, p_R, W_S, W_R, W_P):
    p_P = max(0.0, 1.0 - p_S - p_R)
    wb = w_bar(p_S, p_R, W_S, W_R, W_P)
    dp_S = p_S * (W_S - wb)
    dp_R = p_R * (W_R - wb)
    dp_P = p_P * (W_P - wb)
    return dp_S, dp_R, dp_P

def simulate(p_S0, p_R0, X, N0, dt=0.1, t_max=500.0, K=1e6):
    p_S, p_R = p_S0, p_R0
    N = N0
    t = 0.0
    history = [(t, p_S, p_R, max(0, 1-p_S-p_R), N)]

    while t < t_max:
        W_S, W_R, W_P = fitness_all(p_S, p_R, X, syn_score)
        wb = w_bar(p_S, p_R, W_S, W_R, W_P)
        N = N * np.exp(wb * dt * (1 - N / K))
        if N < 1e-9:
            return float('inf'), history
        dp_S, dp_R, dp_P = delta_p(p_S, p_R, W_S, W_R, W_P)
        p_S = np.clip(p_S + dp_S * dt, 0.0, 1.0)
        p_R = np.clip(p_R + dp_R * dt, 0.0, 1.0)
        p_P_new = max(0.0, 1.0 - p_S - p_R)
        norm = p_S + p_R + p_P_new
        if norm > 0:
            p_S /= norm
            p_R /= norm
        t += dt
        if N > N0:
            return t, history
        history.append((t, p_S, p_R, max(0, 1 - p_S - p_R), N))

    return float('inf'), history

def optimize(p_S0, p_R0, N0):
    best_X, best_time = 0.0, 0.0
    elim_candidates = []
    results = []

    for X in np.linspace(0, 1, 101):
        W_S, W_R, W_P = fitness_all(p_S0, p_R0, X, syn_score)
        wb_init = w_bar(p_S0, p_R0, W_S, W_R, W_P)
        W_S_pure, _, _ = fitness_all(1, 0, X, syn_score)
        _, W_R_pure, _ = fitness_all(0, 1, X, syn_score)
        _, _, W_P_pure = fitness_all(0, 0, X, syn_score)
        can_elim = (W_S_pure < 0) and (W_R_pure < 0) and (W_P_pure < 0)

        if can_elim:
            elim_candidates.append((X, wb_init))
            results.append((X, float('inf'), W_S, W_R, W_P, wb_init))
            continue

        t_pod, _ = simulate(p_S0, p_R0, X, N0)
        results.append((X, t_pod, W_S, W_R, W_P, wb_init))
        if t_pod > best_time:
            best_time = t_pod
            best_X = X

    elim_possible = len(elim_candidates) > 0
    if elim_possible:
        best_X = min(elim_candidates, key=lambda x: x[1])[0]
    return elim_possible, best_X, best_time, results

# ─────────────────────────────────────────────
# Plotting Functions
# ─────────────────────────────────────────────
def ternary_to_cart(p_S, p_R, p_P):
    x = 0.5 * p_S + p_P
    y = (np.sqrt(3) / 2) * p_S
    return x, y

def make_ternary_fig(opt_X, N0, p_S0, p_R0):
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.patch.set_facecolor('#fafaf8')
    fig.suptitle("HDR Replicator Dynamics — Sensitive / Resistant / DTP",
                 fontsize=13, fontweight='bold', fontfamily='monospace', color='#1a1a2e')

    # ── Left: Ternary plot ──
    ax = axes[0]
    ax.set_facecolor('#fafaf8')
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title(f"Population Trajectories on Simplex\n(Drug 1 fraction X = {opt_X:.2f})",
                 fontsize=11, fontfamily='monospace')

    corners = np.array([[0, 0], [1, 0], [0.5, np.sqrt(3) / 2]])
    triangle = plt.Polygon(corners, fill=True, facecolor='#f0f0ea', edgecolor='#1a1a2e', linewidth=1.8)
    ax.add_patch(triangle)

    offset = 0.07
    ax.text(0.5,  np.sqrt(3)/2 + offset, "S\n(Sensitive/Dove)",
            ha='center', va='bottom', fontsize=9, color='#1565c0', fontweight='bold')
    ax.text(-offset - 0.04, -offset, "R\n(Resistant/Hawk)",
            ha='center', va='top', fontsize=9, color='#c62828', fontweight='bold')
    ax.text(1 + offset + 0.04, -offset, "P\n(DTP/Retaliator)",
            ha='center', va='top', fontsize=9, color='#2e7d32', fontweight='bold')

    for v in [0.25, 0.5, 0.75]:
        for i in range(3):
            p1 = corners[i] * v + corners[(i+1) % 3] * (1 - v)
            p2 = corners[i] * v + corners[(i+2) % 3] * (1 - v)
            ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color='#aaaaaa', linewidth=0.4, alpha=0.6)

    init_points = [
        (p_S0, p_R0), (0.8, 0.1), (0.1, 0.8),
        (0.1, 0.1), (0.5, 0.4), (0.4, 0.1),
        (0.1, 0.4), (0.6, 0.3), (0.3, 0.6),
    ]
    for idx, (ps, pr) in enumerate(init_points):
        if ps + pr > 0.98:
            continue
        _, hist = simulate(ps, pr, opt_X, N0, dt=0.2, t_max=300)
        xs = [ternary_to_cart(h[1], h[2], h[3])[0] for h in hist]
        ys = [ternary_to_cart(h[1], h[2], h[3])[1] for h in hist]
        is_actual = (ps == p_S0 and pr == p_R0)
        color = '#1a1a2e' if is_actual else f'C{idx}'
        lw = 2.5 if is_actual else 0.9
        alpha = 1.0 if is_actual else 0.5
        ax.plot(xs, ys, color=color, lw=lw, alpha=alpha)
        ax.plot(xs[0], ys[0], 'o', color=color, markersize=5 if is_actual else 3)
        ax.plot(xs[-1], ys[-1], 's', color=color, markersize=5 if is_actual else 3)

    ax.set_xlim(-0.18, 1.18)
    ax.set_ylim(-0.18, np.sqrt(3)/2 + 0.18)

    legend_elements = [
        Line2D([0], [0], color='#1a1a2e', lw=2.5, label='Initial condition'),
        Line2D([0], [0], marker='o', color='gray', lw=0, markersize=6, label='Start'),
        Line2D([0], [0], marker='s', color='gray', lw=0, markersize=6, label='End'),
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=8, framealpha=0.8)

    # ── Right: Population fractions ──
    ax2 = axes[1]
    ax2.set_facecolor('#fafaf8')
    ax2.set_title(f"Population Fractions Over Time\n(Optimal X = {opt_X:.2f})",
                  fontsize=11, fontfamily='monospace')

    _, hist = simulate(p_S0, p_R0, opt_X, N0, dt=0.1, t_max=500)
    times = [h[0] for h in hist]
    ps_t  = [h[1] for h in hist]
    pr_t  = [h[2] for h in hist]
    pp_t  = [h[3] for h in hist]
    N_t   = [h[4] for h in hist]

    ax2.plot(times, ps_t, color='#1565c0', lw=2.2, label='Sensitive (S/Dove)')
    ax2.plot(times, pr_t, color='#c62828', lw=2.2, label='Resistant (R/Hawk)')
    ax2.plot(times, pp_t, color='#2e7d32', lw=2.2, label='DTP (P/Retaliator)')
    ax2.fill_between(times, ps_t, alpha=0.07, color='#1565c0')
    ax2.fill_between(times, pr_t, alpha=0.07, color='#c62828')
    ax2.fill_between(times, pp_t, alpha=0.07, color='#2e7d32')

    ax2.set_xlabel("Time (generations)", fontsize=10)
    ax2.set_ylabel("Population Fraction", fontsize=10)
    ax2.set_ylim(0, 1.05)
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.25, linestyle='--')

    final_ps, final_pr, final_pp = ps_t[-1], pr_t[-1], pp_t[-1]
    for val, color in [(final_ps, '#1565c0'), (final_pr, '#c62828'), (final_pp, '#2e7d32')]:
        ax2.axhline(y=val, color=color, lw=0.7, ls='--', alpha=0.5)
    ax2.text(times[-1]*0.97, final_ps+0.02, f'{final_ps:.2f}', color='#1565c0', fontsize=8, ha='right')
    ax2.text(times[-1]*0.97, final_pr+0.02, f'{final_pr:.2f}', color='#c62828', fontsize=8, ha='right')
    ax2.text(times[-1]*0.97, final_pp+0.02, f'{final_pp:.2f}', color='#2e7d32', fontsize=8, ha='right')

    ax3 = ax2.twinx()
    ax3.plot(times, N_t, color='#555', lw=1.2, ls='--', alpha=0.45, label='Total N')
    ax3.axhline(y=N0, color='#555', lw=0.7, ls=':', alpha=0.4)
    ax3.set_ylabel("Total Cell Count", fontsize=9, color='gray')
    ax3.tick_params(axis='y', labelcolor='gray')
    ax3.legend(loc='lower right', fontsize=8)

    plt.tight_layout()
    return fig

def make_sweep_fig(results, opt_X):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    fig.patch.set_facecolor('#fafaf8')
    ax1.set_facecolor('#fafaf8')
    ax2.set_facecolor('#fafaf8')
    fig.suptitle("Drug Schedule Optimization Sweep", fontsize=13, fontweight='bold',
                 fontfamily='monospace', color='#1a1a2e')

    Xs       = [r[0] for r in results]
    t_pods   = [r[1] if r[1] != float('inf') else None for r in results]
    W_S_vals = [r[2] for r in results]
    W_R_vals = [r[3] for r in results]
    W_P_vals = [r[4] for r in results]

    t_finite = [(x, t) for x, t in zip(Xs, t_pods) if t is not None]
    t_elim   = [x for x, t in zip(Xs, t_pods) if t is None]

    if t_finite:
        ax1.scatter([t[0] for t in t_finite], [t[1] for t in t_finite],
                    color='steelblue', s=22, alpha=0.8, label='Finite POD', zorder=3)
    if t_elim:
        ax1.scatter(t_elim,
                    [max([t[1] for t in t_finite], default=500)] * len(t_elim),
                    color='gold', s=50, marker='*', label='Elimination possible', zorder=4)

    ax1.axvline(x=opt_X, color='#c62828', lw=1.8, ls='--', label=f'Optimal X = {opt_X:.2f}')
    ax1.set_ylabel("Time to Progression", fontsize=10)
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.25, linestyle='--')

    ax2.plot(Xs, W_S_vals, color='#1565c0', lw=2.2, label='W_S (Sensitive)')
    ax2.plot(Xs, W_R_vals, color='#c62828', lw=2.2, label='W_R (Resistant)')
    ax2.plot(Xs, W_P_vals, color='#2e7d32', lw=2.2, label='W_P (DTP/Persister)')
    ax2.axhline(y=0, color='black', lw=0.9, ls=':')
    ax2.axvline(x=opt_X, color='#c62828', lw=1.8, ls='--')
    ax2.fill_between(Xs, W_S_vals, 0, alpha=0.05, color='#1565c0')
    ax2.fill_between(Xs, W_R_vals, 0, alpha=0.05, color='#c62828')
    ax2.set_xlabel("Fraction of time on Drug 1 (X)", fontsize=10)
    ax2.set_ylabel("Strategy Fitness W_i(X)", fontsize=10)
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.25, linestyle='--')

    plt.tight_layout()
    return fig

# ─────────────────────────────────────────────
# Main Results Area
# ─────────────────────────────────────────────
if not run_btn:
    st.info("👈 Configure parameters in the sidebar, then click **▶ Run Simulation** to begin.")

    # Show equations summary
    with st.expander("Mathematical Equations"):
        st.markdown("""
| # | Equation | Description |
|---|----------|-------------|
| 1 | `base_i(X) = X·λ_i + (1-X)·μ_i − a_syn·X·(1-X)` |  Base Fitness (Under drugs) |
| 2 | `W_i = base_i(X) + Σ p_j · A[i,j]` | Interactions Between Populations and Game Payoff |
| 3 | ``w̄ = Σ p_i · W_i` | Average Population Fitness |
| 4 | `Δp_i = p_i · (W_i − w̄)` | Replicator Equation |
| 5 | `p_i(t+Δt) = p_i(t) + Δt·Δp_i` | Euler's Method to Update Cell Fractions  |
| 6 | `N(t+Δt) = N·exp(w̄·Δt·(1−N/K)` | Ricker Model for Updating Population Size |
""")
else:
    # ── Run ──
    N0   = s0 + r0 + p0_count
    p_S0 = s0 / N0
    p_R0 = r0 / N0
    p_P0 = max(0.0, 1.0 - p_S0 - p_R0)

    with st.spinner("Iterating between drug schedules..."):
        elim_poss, opt_X, max_t, all_results = optimize(p_S0, p_R0, N0)

    W_S_opt, W_R_opt, W_P_opt = fitness_all(p_S0, p_R0, opt_X, syn_score)
    wb_opt = w_bar(p_S0, p_R0, W_S_opt, W_R_opt, W_P_opt)
    dp_S, dp_R, dp_P = delta_p(p_S0, p_R0, W_S_opt, W_R_opt, W_P_opt)

    # ── Summary Metrics ──
    st.markdown("## Results")

    col1, col2, col3, col4 = st.columns(4)
    col1.metric("Optimal Drug 1 Fraction (X)", f"{opt_X:.2f}")
    col2.metric("Drug 2 Fraction (1-X)", f"{1-opt_X:.2f}")
    col3.metric("Avg Fitness w̄", f"{wb_opt:+.5f}")
    col4.metric("Time to POD", "∞ (Elim.)" if elim_poss else f"{round(max_t, 1)}")

    st.divider()

    # ── Elimination badge ──
    if elim_poss:
        st.markdown(
            '<span class="elim-badge elim-yes">Elimination Possible</span>',
            unsafe_allow_html=True)
    else:
        st.markdown(
            f'<span class="elim-badge elim-no">Elimination Is Not Possible — Max time to POD: {round(max_t, 2)}</span>',
            unsafe_allow_html=True)

    st.markdown("")

    # ── Fitness & Replicator Dynamics Table ──
    col_left, col_right = st.columns(2)

    with col_left:
        st.markdown("### Fitness at Optimal X")
        st.markdown(f"""
<div class="metric-card positive">W_S (Sensitive / Dove):  &nbsp;<strong>{W_S_opt:+.5f}</strong></div>
<div class="metric-card positive">W_R (Resistant / Hawk):  &nbsp;<strong>{W_R_opt:+.5f}</strong></div>
<div class="metric-card neutral">W_P (DTP / Retaliator):  &nbsp;<strong>{W_P_opt:+.5f}</strong></div>
<div class="metric-card negative">w̄  (avg population):    &nbsp;<strong>{wb_opt:+.5f}</strong></div>
""", unsafe_allow_html=True)

    with col_right:
        st.markdown("### Replicator Dynamics at t=0")

        def trend(val, pos_label, neg_label):
            if val > 0:
                return f"▲ {pos_label}"
            elif val < 0:
                return f"▼ {neg_label}"
            return "◆ Equilibrium"

        st.markdown(f"""
<div class="metric-card positive">
  Δp_S = <strong>{dp_S:+.6f}</strong> — {trend(dp_S, "Sensitive RECOVERING", "Sensitive SHRINKING")}
</div>
<div class="metric-card positive">
  Δp_R = <strong>{dp_R:+.6f}</strong> — {trend(dp_R, "Resistance SPREADING", "Resistance SHRINKING")}
</div>
<div class="metric-card neutral">
  Δp_P = <strong>{dp_P:+.6f}</strong> — {trend(dp_P, "DTP EXPANDING", "DTP SHRINKING")}
</div>
""", unsafe_allow_html=True)

    st.divider()

    # ── Interpretation ──
    st.markdown("### Interpretation")
    if W_R_opt < W_S_opt and W_R_opt < W_P_opt:
        st.success("**Resistant (Hawk) has the lowest fitness** — the drug schedule suppresses Hawks. DTPs (Retaliators) act as a buffer against resistance re-emergence.")
    elif W_P_opt > W_R_opt:
        st.info("**DTPs have fitness advantage over Resistant** at this schedule — consistent with Retaliator suppressing Hawk in game theory.")
    else:
        st.warning("**Resistant cells remain competitive.** Consider adjusting interaction parameters or increasing drug doses.")

    st.divider()

    # ── Plots ──
    st.markdown("### 📈 Visuals")

    tab1, tab2 = st.tabs(["Ternary Simplex & Population Dynamics", "Drug Schedule Sweep"])

    with tab1:
        with st.spinner("Rendering ternary plot..."):
            fig1 = make_ternary_fig(opt_X, N0, p_S0, p_R0)
        st.pyplot(fig1, use_container_width=True)

    with tab2:
        with st.spinner("Rendering drug sweep..."):
            fig2 = make_sweep_fig(all_results, opt_X)
        st.pyplot(fig2, use_container_width=True)

    # ── Initial Conditions Summary ──
    with st.expander(" Initial Conditions"):
        st.markdown(f"""
| Parameter | Value |
|-----------|-------|
| p_S (Sensitive) | `{p_S0:.4f}` |
| p_R (Resistant) | `{p_R0:.4f}` |
| p_P (DTP) | `{p_P0:.4f}` |
| Total N₀ | `{N0:.0f}` |
| Synergy Score | `{syn_score}` |
""")
