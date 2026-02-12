import streamlit as st
import numpy as np

# ============================================================
# 1. APP TITLE
# ============================================================
st.title("Synergy-Based Optimal Drug Scheduling Model")

# ============================================================
# 2. USER INPUTS
# ============================================================

st.header("Synergy Score")
syn_score = st.number_input("Enter synergy score:", value=1.0)

st.header("Sensitive Cell Measurements (at dosing concentration)")
s0 = st.number_input("Initial Sensitive Cell Count (s0)", value=1e6)
lambda_s = st.number_input("Net Growth Rate under Drug 1 (λ_S)", value=-0.05)
mu_s = st.number_input("Net Growth Rate under Drug 2 (μ_S)", value=-0.02)

st.header("Resistant Cell Measurements (at dosing concentration)")
r0 = st.number_input("Initial Resistant Cell Count (r0)", value=1e4)
lambda_r = st.number_input("Net Growth Rate under Drug 1 (λ_R)", value=0.03)
mu_r = st.number_input("Net Growth Rate under Drug 2 (μ_R)", value=0.01)

# Synergy scaling constant
SYNERGY_SCALE = 0.005

# ============================================================
# 3. CORE MODEL FUNCTION
# ============================================================

def calculate_plan():
    best_time = 0.0
    best_X = 0.0
    elimination_possible = False

    for X in np.linspace(0, 1, 101):
        penalty = SYNERGY_SCALE * syn_score * X * (1 - X)

        aS = (X * lambda_s + (1 - X) * mu_s) - penalty
        aR = (X * lambda_r + (1 - X) * mu_r) - penalty

        # Elimination condition
        if aS < 0 and aR < 0:
            elimination_possible = True
            return True, X, float("inf")

        # Time to progression (POD)
        if aR > aS and aS < 0 and aR > 0:
            try:
                t_pod = (1 / (aR - aS)) * np.log(-s0 * aS / (r0 * aR))
                if t_pod > best_time:
                    best_time = t_pod
                    best_X = X
            except:
                pass

        elif aR <= 0 and aS <= 0:
            if 0 > aR > -1e-9:
                t_pod = 1000.0
                if t_pod > best_time:
                    best_time = t_pod
                    best_X = X

    return elimination_possible, best_X, best_time


# ============================================================
# 4. RUN BUTTON
# ============================================================

if st.button("Run Optimization"):
    elim_poss, opt_x, max_t = calculate_plan()

    st.subheader("Results")

    st.write(f"**Elimination Possible:** {elim_poss}")
    st.write(f"**Optimal Fraction Drug 1:** {round(opt_x, 2)}")
    st.write(f"**Optimal Fraction Drug 2:** {round(1 - opt_x, 2)}")

    if max_t == float("inf"):
        st.write("**Max Time to Progression:** ∞ (Elimination achieved)")
    else:
        st.write(f"**Max Time to Progression:** {round(max_t, 5)}")
