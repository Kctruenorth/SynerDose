import streamlit as st
import pandas as pd
import time
from google.oauth2.service_account import Credentials
from googleapiclient.discovery import build

st.set_page_config(
    page_title="Contribute — Drug Synergy",
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

    .info-box {
        background: rgba(0, 212, 255, 0.05);
        border: 1px solid #1e2d45;
        padding: 1rem 1.25rem;
        border-radius: 4px;
        font-size: 0.85rem;
        color: #8899bb;
        font-family: 'IBM Plex Mono', monospace;
        margin-bottom: 1.5rem;
        line-height: 1.7;
    }

    div[data-testid="stSelectbox"] label,
    div[data-testid="stTextInput"] label,
    div[data-testid="stNumberInput"] label,
    div[data-testid="stSlider"] label {
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

# ── Google Sheets Setup ───────────────────────────────────────────────────────
SCOPE = ['https://www.googleapis.com/auth/spreadsheets']

def get_sheet_service():
    credentials_info = dict(st.secrets["google_credentials"])
    credentials = Credentials.from_service_account_info(credentials_info, scopes=SCOPE)
    service = build("sheets", "v4", credentials=credentials)
    return service.spreadsheets()

SHEET_ID = "1e0jkPmkRL3VuBg6Ei_RWtzRPsI00RMS7cMZyg5rqWOs"  

def append_to_sheet(data: dict):
    sheet = get_sheet_service()
    values = [[
        data["CELL_LINE"],
        data["COMPOUND_A"],
        data["COMPOUND_B"],
        data["MAX_CONC_A"],
        data["MAX_CONC_B"],
        data["IC50_A"],
        data["H_A"],
        data["Einf_A"],
        data["IC50_B"],
        data["H_B"],
        data["Einf_B"],
        data["SYNERGY_SCORE"],
    ]]
    body = {"values": values}
    sheet.values().append(
        spreadsheetId=SHEET_ID,
        range="Sheet1!A:L",
        valueInputOption="RAW",
        body=body
    ).execute()

# ── Header ────────────────────────────────────────────────────────────────────
st.markdown("""
<div class="header-block">
    <h1 style="margin:0; font-size:1.6rem;">CONTRIBUTE DATA</h1>
    <p style="margin:0.3rem 0 0; color:#8899bb; font-size:0.85rem; font-family:'IBM Plex Mono',monospace;">
        Help improve the model by submitting real experimental results
    </p>
</div>
""", unsafe_allow_html=True)

st.markdown("""
<div class="info-box">
    By submitting your experimental data, you help expand the training dataset.<br>
    All fields are required. The synergy score should be the experimentally observed Loewe score.
</div>
""", unsafe_allow_html=True)

# ── Cell Line & Compounds ─────────────────────────────────────────────────────
st.markdown('<p class="section-label">Cell Line & Compounds</p>', unsafe_allow_html=True)

cell_line  = st.text_input("Cell Line", placeholder="e.g. HT-29").strip().upper()
col1, col2 = st.columns(2)
with col1:
    compound_a = st.text_input("Compound A", placeholder="e.g. MAP2K_1").strip()
with col2:
    compound_b = st.text_input("Compound B", placeholder="e.g. TKI").strip()

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

# ── Observed Synergy Score ────────────────────────────────────────────────────
st.markdown('<p class="section-label">Observed Outcome</p>', unsafe_allow_html=True)
synergy_score = st.number_input(
    "Observed Synergy Score (Loewe)",
    value=0.0,
    step=0.1,
    help="The experimentally measured synergy score. ≥ 20 = synergistic."
)

# ── Submit ────────────────────────────────────────────────────────────────────
st.markdown("---")
if st.button("SUBMIT DATA"):
    if not cell_line or not compound_a or not compound_b:
        st.error("Please fill in Cell Line, Compound A, and Compound B.")
    else:
        with st.spinner("Submitting..."):
            time.sleep(1)
            try:
                append_to_sheet({
                    "CELL_LINE":     cell_line,
                    "COMPOUND_A":    compound_a,
                    "COMPOUND_B":    compound_b,
                    "MAX_CONC_A":    max_conc_a,
                    "MAX_CONC_B":    max_conc_b,
                    "IC50_A":        ic50_a,
                    "H_A":           h_a,
                    "Einf_A":        einf_a,
                    "IC50_B":        ic50_b,
                    "H_B":           h_b,
                    "Einf_B":        einf_b,
                    "SYNERGY_SCORE": synergy_score,
                })
                st.success("Data submitted successfully — thank you for contributing!")
            except Exception as e:
                st.error(f"Submission failed: {e}")
