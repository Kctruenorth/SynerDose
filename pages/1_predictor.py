import streamlit as st
import joblib
import numpy as np
import pandas as pd


st.set_page_config(
    page_title="Drug Synergy Predictor",
    page_icon="🧬",
    layout="centered"
)

st.markdown("""
<style>
    @import url('https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;600&family=IBM+Plex+Sans:wght@300;400;600&display=swap');

    html, body, [class*="css"] {
        font-family: 'IBM Plex Sans', sans-serif;
    }

    .stApp {
        background-color: #0a0e1a;
        color: #e0e6f0;
    }

    h1, h2, h3 {
        font-family: 'IBM Plex Mono', monospace !important;
        color: #00d4ff !important;
    }

    .header-block {
        border-left: 3px solid #00d4ff;
        padding: 0.4rem 1rem;
        margin-bottom: 1.5rem;
        background: rgba(0, 212, 255, 0.05);
    }

    .result-card {
        padding: 1.5rem 2rem;
        border-radius: 4px;
        margin-top: 1.5rem;
        font-family: 'IBM Plex Mono', monospace;
    }

    .synergistic {
        background: rgba(0, 255, 150, 0.08);
        border: 1px solid #00ff96;
        color: #00ff96;
    }

    .not-synergistic {
        background: rgba(255, 80, 80, 0.08);
        border: 1px solid #ff5050;
        color: #ff5050;
    }

    .score-value {
        font-size: 2.5rem;
        font-weight: 600;
        letter-spacing: -1px;
    }

    .score-label {
        font-size: 0.75rem;
        opacity: 0.6;
        text-transform: uppercase;
        letter-spacing: 2px;
        margin-bottom: 0.25rem;
    }

    .warning-box {
        background: rgba(255, 200, 0, 0.08);
        border: 1px solid #ffc800;
        color: #ffc800;
        padding: 0.75rem 1rem;
        border-radius: 4px;
        font-size: 0.85rem;
        margin-top: 0.5rem;
        font-family: 'IBM Plex Mono', monospace;
    }

    .section-label {
        font-family: 'IBM Plex Mono', monospace;
        font-size: 0.7rem;
        text-transform: uppercase;
        letter-spacing: 3px;
        color: #00d4ff;
        opacity: 0.7;
        margin-bottom: 0.5rem;
        margin-top: 1.5rem;
    }

    div[data-testid="stSelectbox"] label,
    div[data-testid="stTextInput"] label,
    div[data-testid="stNumberInput"] label {
        font-family: 'IBM Plex Mono', monospace !important;
        font-size: 0.78rem !important;
        color: #8899bb !important;
        text-transform: uppercase;
        letter-spacing: 1px;
    }

    div[data-testid="stSelectbox"] > div > div,
    div[data-testid="stTextInput"] > div > div > input,
    div[data-testid="stNumberInput"] > div > div > input {
        background-color: #111827 !important;
        border: 1px solid #1e2d45 !important;
        color: #e0e6f0 !important;
        border-radius: 3px !important;
        font-family: 'IBM Plex Mono', monospace !important;
    }

    div[data-testid="stSelectbox"] > div > div:focus-within,
    div[data-testid="stTextInput"] > div > div > input:focus,
    div[data-testid="stNumberInput"] > div > div > input:focus {
        border-color: #00d4ff !important;
        box-shadow: 0 0 0 1px #00d4ff22 !important;
    }

    .stButton > button {
        background: #00d4ff !important;
        color: #0a0e1a !important;
        border: none !important;
        border-radius: 3px !important;
        font-family: 'IBM Plex Mono', monospace !important;
        font-weight: 600 !important;
        font-size: 0.85rem !important;
        letter-spacing: 2px !important;
        text-transform: uppercase !important;
        padding: 0.6rem 2rem !important;
        width: 100% !important;
        margin-top: 1rem !important;
        transition: all 0.2s ease !important;
    }

    .stButton > button:hover {
        background: #00b8d9 !important;
        transform: translateY(-1px) !important;
    }

    hr {
        border-color: #1e2d45 !important;
        margin: 1.5rem 0 !important;
    }
</style>
""", unsafe_allow_html=True)


@st.cache_resource
def load_models():
    model          = joblib.load("xgboost_synergy_model.pkl")
    encoder        = joblib.load("target_encoder.pkl")
    final_features = joblib.load("final_features.pkl")
    genetics_df    = joblib.load("genetics_df.pkl")
    return model, encoder, final_features, genetics_df

def add_interactions(df):
    df = df.copy()
    df["IC50_A_log"] = np.log1p(df["IC50_A"])
    df["IC50_B_log"] = np.log1p(df["IC50_B"])
    return df


# ── Header ───────────────────────────────────────────────────────────────────
st.markdown("""
<div class="header-block">
    <h1 style="margin:0; font-size:1.6rem;">DRUG SYNERGY PREDICTOR</h1>
    <p style="margin:0.3rem 0 0; color:#8899bb; font-size:0.85rem; font-family:'IBM Plex Mono',monospace;">
        XGBoost · Clinical + Genomic Features · Loewe Synergy Score
    </p>
</div>
""", unsafe_allow_html=True)

try:
    model, encoder, final_features, genetics_df = load_models()
except Exception as e:
    st.error(f"Failed to load model files: {e}")
    st.stop()

cell_lines = sorted(genetics_df["cell_line_name"].unique().tolist())

# ── Cell Line & Compounds ─────────────────────────────────────────────────────
st.markdown('<p class="section-label">Cell Line & Compounds</p>', unsafe_allow_html=True)

cell_line  = st.selectbox("Cell Line", options=cell_lines, index=cell_lines.index("HT-29") if "HT-29" in cell_lines else 0)
col1, col2 = st.columns(2)
with col1:
    compound_a = st.text_input("Compound A", value="MAP2K_1", placeholder="e.g. MAP2K_1")
with col2:
    compound_b = st.text_input("Compound B", value="TKI", placeholder="e.g. TKI")

# ── Drug A Parameters ─────────────────────────────────────────────────────────
st.markdown('<p class="section-label">Drug A — Pharmacological Parameters</p>', unsafe_allow_html=True)
a1, a2, a3, a4 = st.columns(4)
with a1:
    max_conc_a = st.number_input("Max Conc A (µM)", min_value=0.0, value=1.0, step=0.1, help="Typical: 1 – 3 µM")
with a2:
    ic50_a = st.number_input("IC50 A (µM)", min_value=0.0, value=0.33, step=0.01, help="Typical: 0.06 – 1 µM")
with a3:
    h_a = st.number_input("Hill Coeff A", min_value=0.0, value=1.3, step=0.1, help="Typical: 0.57 – 2.9")
with a4:
    einf_a = st.number_input("Einf A (%)", min_value=0.0, max_value=100.0, value=38.0, step=1.0, help="Typical: 10.4 – 79.2%")

# ── Drug B Parameters ─────────────────────────────────────────────────────────
st.markdown('<p class="section-label">Drug B — Pharmacological Parameters</p>', unsafe_allow_html=True)
b1, b2, b3, b4 = st.columns(4)
with b1:
    max_conc_b = st.number_input("Max Conc B (µM)", min_value=0.0, value=1.0, step=0.1, help="Typical: 1 – 3 µM")
with b2:
    ic50_b = st.number_input("IC50 B (µM)", min_value=0.0, value=0.41, step=0.01, help="Typical: 0.12 – 1 µM")
with b3:
    h_b = st.number_input("Hill Coeff B", min_value=0.0, value=1.2, step=0.1, help="Typical: 0.49 – 2.89")
with b4:
    einf_b = st.number_input("Einf B (%)", min_value=0.0, max_value=100.0, value=45.0, step=1.0, help="Typical: 11 – 80.4%")

# ── Predict ───────────────────────────────────────────────────────────────────
st.markdown("---")
if st.button("RUN PREDICTION"):
    if not compound_a or not compound_b:
        st.error("Please enter both compound names.")
    else:
        sample = pd.DataFrame([{
            "CELL_LINE":  cell_line,
            "COMPOUND_A": compound_a,
            "COMPOUND_B": compound_b,
            "MAX_CONC_A": max_conc_a,
            "MAX_CONC_B": max_conc_b,
            "IC50_A":     ic50_a,
            "H_A":        h_a,
            "Einf_A":     einf_a,
            "IC50_B":     ic50_b,
            "H_B":        h_b,
            "Einf_B":     einf_b
        }])

        sample = add_interactions(sample)
        cell_genetics = genetics_df[genetics_df["cell_line_name"] == cell_line]

        if cell_genetics.empty:
            st.error(f"Cell line '{cell_line}' not found in genetics data.")
        else:
            sample     = pd.merge(sample, cell_genetics, left_on="CELL_LINE", right_on="cell_line_name", how="left")
            sample     = sample[final_features]
            sample_enc = encoder.transform(sample)
            pred       = model.predict(sample_enc)[0]

            is_synergistic = pred >= 20
            card_class     = "synergistic" if is_synergistic else "not-synergistic"
            label          = "SYNERGISTIC" if is_synergistic else "NOT SYNERGISTIC"
            emoji          = "✦" if is_synergistic else "✕"

            st.markdown(f"""
            <div class="result-card {card_class}">
                <div class="score-label">Predicted Loewe Synergy Score</div>
                <div class="score-value">{pred:.2f}</div>
                <div style="margin-top:0.75rem; font-size:1rem; letter-spacing:3px;">
                    {emoji} &nbsp; {label}
                </div>
                <div style="margin-top:0.5rem; font-size:0.75rem; opacity:0.6;">
                    {cell_line} · {compound_a} + {compound_b} · threshold ≥ 20
                </div>
            </div>
            """, unsafe_allow_html=True)
