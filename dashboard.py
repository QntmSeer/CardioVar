"""
dashboard.py — CardioVar v1.2 Production Dashboard
Scientifically approvable Plotly visualizations + clean UX.
"""
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

import json, io, time, os, base64, collections, math
import streamlit as st
import pandas as pd
import numpy as np
import requests
import plotly.graph_objects as go

from plots import (
    plot_variant_impact_profile,
    plot_background_kde,
    plot_pathogenicity_radar,
    plot_evidence_stack,
    plot_gnomad_context,
    plot_clinvar_lollipop,
    plot_tissue_heatmap,
    plot_enformer_tracks,
)

# ── Page config ───────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="CardioVar",
    page_icon="❤️",
    layout="wide",
    initial_sidebar_state="expanded",
)

API_URL = os.getenv("API_URL", "http://localhost:8000")

COLORS = {
    "primary":    "#F46036",
    "secondary":  "#5B85AA",
    "accent":     "#9C7A97",
    "background": "#F4F5F7",
    "card_bg":    "#FFFFFF",
    "text_dark":  "#171123",
    "text_light": "#888DA7",
    "success":    "#2A9D8F",
    "warning":    "#E9C46A",
    "danger":     "#C44027",
}

# ── Global CSS ────────────────────────────────────────────────────────────────
st.markdown(f"""
<style>
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700;800&display=swap');

html, body, [class*="css"] {{
    font-family: 'Inter', sans-serif;
    background-color: {COLORS['background']};
    color: {COLORS['text_dark']};
}}

/* ── Sidebar ─────────────────────────────────────────────────── */
[data-testid="stSidebar"] {{
    background: linear-gradient(175deg, #1a1035 0%, #0e1f3d 100%);
    border-right: 1px solid rgba(255,255,255,0.07);
}}
[data-testid="stSidebar"] * {{
    color: #e8eaf0 !important;
}}
[data-testid="stSidebar"] .stRadio label {{
    padding: 6px 12px;
    border-radius: 8px;
    transition: background 0.15s;
    font-size: 14px;
    font-weight: 500;
}}
[data-testid="stSidebar"] .stRadio label:hover {{
    background: rgba(244,96,54,0.18) !important;
}}
[data-testid="stSidebar"] hr {{
    border-color: rgba(255,255,255,0.1) !important;
}}

/* ── Cards ───────────────────────────────────────────────────── */
.cv-card {{
    background: {COLORS['card_bg']};
    border-radius: 14px;
    padding: 24px 28px;
    box-shadow: 0 2px 12px rgba(0,0,0,0.06);
    margin-bottom: 20px;
    border: 1px solid #E8E8EC;
}}

/* ── Metric tiles ─────────────────────────────────────────────── */
.metric-tile {{
    background: {COLORS['card_bg']};
    border-radius: 12px;
    padding: 16px 20px;
    border: 1px solid #E8E8EC;
    box-shadow: 0 2px 8px rgba(0,0,0,0.04);
    text-align: center;
}}
.metric-label {{
    font-size: 11px;
    text-transform: uppercase;
    letter-spacing: 1.2px;
    color: {COLORS['text_light']};
    margin-bottom: 6px;
    font-weight: 600;
}}
.metric-value {{
    font-size: 26px;
    font-weight: 700;
    color: {COLORS['text_dark']};
    line-height: 1.1;
}}
.metric-sub {{
    font-size: 11px;
    color: {COLORS['text_light']};
    margin-top: 4px;
}}

/* ── Classification badge ─────────────────────────────────────── */
.badge {{
    display: inline-block;
    padding: 6px 18px;
    border-radius: 20px;
    font-weight: 700;
    font-size: 14px;
    letter-spacing: 0.5px;
    text-transform: uppercase;
}}
.badge-pathogenic {{ background: rgba(196,64,39,0.12); color: {COLORS['danger']}; border: 1.5px solid {COLORS['danger']}; }}
.badge-vus        {{ background: rgba(233,196,106,0.15); color: #9C7A00; border: 1.5px solid {COLORS['warning']}; }}
.badge-benign     {{ background: rgba(42,157,143,0.12); color: {COLORS['success']}; border: 1.5px solid {COLORS['success']}; }}

/* ── Section title ─────────────────────────────────────────────── */
.section-title {{
    font-size: 22px;
    font-weight: 800;
    color: {COLORS['text_dark']};
    margin-bottom: 4px;
}}
.section-sub {{
    font-size: 13px;
    color: {COLORS['text_light']};
    margin-bottom: 20px;
}}

/* ── Buttons ────────────────────────────────────────────────────── */
.stButton>button {{
    background: {COLORS['primary']};
    color: white;
    border-radius: 9px;
    font-weight: 600;
    font-size: 14px;
    border: none;
    padding: 0.5rem 1.2rem;
    transition: all 0.2s ease;
}}
.stButton>button:hover {{
    background: #D94E28;
    box-shadow: 0 4px 12px rgba(244,96,54,0.3);
    transform: translateY(-1px);
}}

/* ── Tabs ────────────────────────────────────────────────────────── */
.stTabs [data-baseweb="tab-list"] {{ gap: 20px; }}
.stTabs [data-baseweb="tab"] {{
    height: 46px;
    background: transparent;
    border-radius: 6px 6px 0 0;
    font-weight: 600;
    font-size: 14px;
    padding: 8px 4px;
}}
.stTabs [aria-selected="true"] {{
    color: {COLORS['primary']} !important;
    border-bottom: 2.5px solid {COLORS['primary']} !important;
}}

/* ── Gene card ────────────────────────────────────────────────────── */
.gene-card {{
    background: white;
    border-radius: 14px;
    padding: 20px;
    border: 1.5px solid #EAEAEA;
    box-shadow: 0 2px 8px rgba(0,0,0,0.05);
    transition: all 0.2s ease;
}}
.gene-card:hover {{
    border-color: {COLORS['primary']};
    box-shadow: 0 4px 16px rgba(244,96,54,0.12);
}}
</style>
""", unsafe_allow_html=True)

# ── Session state defaults ──────────────────────────────────────────────────
if "analysis_results" not in st.session_state:
    st.session_state.analysis_results = None
if "inputs" not in st.session_state:
    st.session_state.inputs = {"chrom": "chr22", "pos": 36191400, "ref": "A", "alt": "C"}
# nav_radio is the SINGLE source of truth for navigation.
# Initialize it once; never pass index= to st.radio again.
if "nav_radio" not in st.session_state:
    st.session_state.nav_radio = "Inference"

def set_page(name: str):
    """Programmatic navigation — write directly to the radio key and rerun."""
    st.session_state.nav_radio = name

# ── Sidebar ─────────────────────────────────────────────────────────────────
PAGES = ["Inference", "Gene Panel", "Benchmarking", "MSA Explorer", "Models", "Report", "Settings"]

with st.sidebar:
    # Logo / branding
    try:
        st.image("logo.png", width=54)
    except Exception:
        st.markdown(f"<h2 style='color:#F46036;margin:0'>&#9829;</h2>", unsafe_allow_html=True)

    st.markdown("""
        <div style="margin: 4px 0 16px 0;">
            <span style="font-size:22px;font-weight:800;color:white;">CARDIO<span style="color:#F46036;">VAR</span></span><br>
            <span style="font-size:11px;color:#888DA7;letter-spacing:1.5px;text-transform:uppercase;">Precision Genomics</span>
        </div>
    """, unsafe_allow_html=True)

    # No index= — key alone controls the radio; avoids the 2-click lag
    selected_page = st.radio(
        "Navigate",
        PAGES,
        key="nav_radio",
        label_visibility="collapsed",
    )

    st.markdown("---")

    # Workspace
    st.markdown("<span style='font-size:12px;font-weight:600;text-transform:uppercase;letter-spacing:1px;color:#888DA7;'>Workspace</span>", unsafe_allow_html=True)

    if st.button("💾 Save Workspace", use_container_width=True):
        state_dump = {
            "inputs": st.session_state.get("inputs", {}),
            "analysis_results": st.session_state.get("analysis_results", {}),
            "page": st.session_state.get("page", "Inference"),
        }
        st.download_button("⬇ Download JSON", json.dumps(state_dump, indent=2),
                           "cardiovar_workspace.json", "application/json",
                           use_container_width=True)

    ws_file = st.file_uploader("📂 Load Workspace", type=["json"], label_visibility="collapsed")
    if ws_file:
        try:
            loaded = json.load(ws_file)
            st.session_state.inputs.update(loaded.get("inputs", {}))
            st.session_state.analysis_results = loaded.get("analysis_results", {})
            st.session_state.page = loaded.get("page", "Inference")
            st.success("Workspace loaded!")
            time.sleep(0.5)
            st.rerun()
        except Exception as e:
            st.error(f"Load failed: {e}")

    st.markdown("---")
    st.markdown(f"<div style='font-size:11px;color:#555;'>v1.2.0 · GRCh38 · Enformer</div>",
                unsafe_allow_html=True)

# ══════════════════════════════════════════════════════════════════════════════
# Page: Inference
# ══════════════════════════════════════════════════════════════════════════════
if selected_page == "Inference":
    st.markdown("<div class='section-title'>Variant Inference</div>", unsafe_allow_html=True)
    st.markdown("<div class='section-sub'>Enter a genomic variant (GRCh38) to predict regulatory impact via Enformer.</div>", unsafe_allow_html=True)

    # Input card
    with st.container():
        st.markdown("<div class='cv-card'>", unsafe_allow_html=True)
        demo_col, _ = st.columns([1, 4])
        with demo_col:
            if st.button("⚡ Load Demo (MYH9)"):
                st.session_state.inputs.update({"chrom": "chr22", "pos": 36191400, "ref": "A", "alt": "C"})
                st.rerun()

        c1, c2, c3, c4, c5 = st.columns([2, 2, 1, 1, 2])
        CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
        with c1:
            chrom = st.selectbox("Chromosome", CHROMS,
                                 index=CHROMS.index(st.session_state.inputs.get("chrom", "chr22")))
        with c2:
            pos   = st.number_input("Position (bp)", min_value=1,
                                    value=int(st.session_state.inputs.get("pos", 36191400)), step=1)
        with c3:
            ref   = st.text_input("Ref", st.session_state.inputs.get("ref", "A")).upper().strip()
        with c4:
            alt   = st.text_input("Alt", st.session_state.inputs.get("alt", "C")).upper().strip()
        with c5:
            st.markdown("<div style='height:28px'></div>", unsafe_allow_html=True)
            run_btn = st.button("▶ Run Analysis", use_container_width=True)
        st.markdown("</div>", unsafe_allow_html=True)

    st.session_state.inputs.update({"chrom": chrom, "pos": pos, "ref": ref, "alt": alt})

    if run_btn:
        if not ref or not alt:
            st.error("Reference and Alternate alleles are required.")
        elif len(ref) > 50 or len(alt) > 50:
            st.error("Allele strings should be ≤50 bp.")
        else:
            payload = {"assembly": "GRCh38", "chrom": chrom, "pos": pos,
                       "ref": ref, "alt": alt, "force_live": False}
            try:
                with st.spinner("Running Enformer inference — this may take 30–60 s on first run…"):
                    resp = requests.post(f"{API_URL}/variant-impact", json=payload, timeout=300)
                    resp.raise_for_status()
                    st.session_state.analysis_results = resp.json()
            except Exception as e:
                st.error(f"Analysis failed: {e}")

    # ── Results ───────────────────────────────────────────────────────────────
    if st.session_state.analysis_results:
        data    = st.session_state.analysis_results
        metrics = data.get("metrics", {})
        var_id  = f"{chrom}:{pos} {ref}→{alt}"

        # ── Evidence stack → classification badge ─────────────────────────
        try:
            _, label, badge_color = plot_evidence_stack(data)
        except Exception:
            label, badge_color = "VUS", COLORS["warning"]

        badge_cls = {
            "LIKELY PATHOGENIC": "badge-pathogenic",
            "VUS": "badge-vus",
            "LIKELY BENIGN": "badge-benign",
        }.get(label, "badge-vus")

        # Header row: variant ID + badge
        hc1, hc2 = st.columns([3, 1])
        with hc1:
            st.markdown(f"<div style='font-size:18px;font-weight:700;color:{COLORS['text_dark']};margin-bottom:4px;'>Results — {var_id}</div>", unsafe_allow_html=True)
            model_used = metrics.get("model_used", "Unknown")
            st.markdown(f"<div style='font-size:12px;color:{COLORS['text_light']};'>Model: {model_used} · {time.strftime('%Y-%m-%d %H:%M')}</div>", unsafe_allow_html=True)
        with hc2:
            st.markdown(f"<div style='text-align:right;margin-top:8px;'><span class='badge {badge_cls}'>{label}</span></div>", unsafe_allow_html=True)

        # Metrics tiles
        conf = metrics.get("confidence", 0.0)
        conf_color = COLORS["success"] if conf > 80 else COLORS["warning"] if conf > 50 else COLORS["danger"]
        zscore_color = COLORS["danger"] if abs(metrics.get("z_score", 0)) > 2 else COLORS["text_dark"]

        m1, m2, m3, m4, m5, m6 = st.columns(6)
        for col, lbl, val, sub in [
            (m1, "Gene",        metrics.get("gene_symbol", "N/A"), "Symbol"),
            (m2, "Max |Δ|",     f"{abs(metrics.get('max_delta', 0)):.4f}", "Enformer"),
            (m3, "Z-Score",     f"{metrics.get('z_score', 0):.2f}", "vs background"),
            (m4, "Percentile",  f"{metrics.get('percentile', 0):.1f}%", "in gene"),
            (m5, "gnomAD Freq", f"{metrics.get('gnomad_freq', 0):.2e}", "population"),
            (m6, "Confidence",  f"{conf:.0f}%", "data quality"),
        ]:
            with col:
                st.markdown(f"""
                <div class='metric-tile'>
                    <div class='metric-label'>{lbl}</div>
                    <div class='metric-value'>{val}</div>
                    <div class='metric-sub'>{sub}</div>
                </div>""", unsafe_allow_html=True)

        st.markdown("<br>", unsafe_allow_html=True)

        # ── Result tabs ───────────────────────────────────────────────────────
        t_impact, t_stats, t_popn, t_tissue, t_gene, t_related = st.tabs([
            "🔬 Impact Profile",
            "📊 Statistical Context",
            "🌍 Population & Clinical",
            "🫀 Tissue Effects",
            "🧬 Gene Context",
            "📋 Related Data",
        ])

        # ── Tab: Impact Profile ────────────────────────────────────────────
        with t_impact:
            st.markdown("<div class='cv-card'>", unsafe_allow_html=True)
            st.markdown("**Plot A — Multi-Track Genomic Browser** | Enformer ΔRNA-seq signal with gene structure and PhyloP conservation")
            try:
                fig_a = plot_variant_impact_profile(data)
                st.plotly_chart(fig_a, use_container_width=True)
            except Exception as e:
                st.error(f"Plot A failed: {e}")
            st.markdown("</div>", unsafe_allow_html=True)

        # ── Tab: Statistical Context ───────────────────────────────────────
        with t_stats:
            col_b, col_c = st.columns(2)
            with col_b:
                st.markdown("<div class='cv-card'>", unsafe_allow_html=True)
                st.markdown("**Plot B — Background Distribution** | Gene-specific KDE with variant Z-score")
                try:
                    fig_b = plot_background_kde(data)
                    st.plotly_chart(fig_b, use_container_width=True)
                except Exception as e:
                    st.error(f"Plot B failed: {e}")
                st.markdown("</div>", unsafe_allow_html=True)
            with col_c:
                st.markdown("<div class='cv-card'>", unsafe_allow_html=True)
                st.markdown("**Plot C — Pathogenicity Radar** | Multi-evidence profile vs known pathogenic median")
                try:
                    fig_c = plot_pathogenicity_radar(data)
                    st.plotly_chart(fig_c, use_container_width=True)
                except Exception as e:
                    st.error(f"Plot C failed: {e}")
                st.markdown("</div>", unsafe_allow_html=True)

            st.markdown("<div class='cv-card'>", unsafe_allow_html=True)
            st.markdown("**Plot D — Evidence Breakdown**")
            try:
                fig_d, _, _ = plot_evidence_stack(data)
                st.plotly_chart(fig_d, use_container_width=True)
            except Exception as e:
                st.error(f"Plot D failed: {e}")
            st.markdown("</div>", unsafe_allow_html=True)

        # ── Tab: Population & Clinical ─────────────────────────────────────
        with t_popn:
            col_e, col_f = st.columns(2)
            with col_e:
                st.markdown("<div class='cv-card'>", unsafe_allow_html=True)
                st.markdown("**Plot E — gnomAD Allele Frequency Context** | Population distribution")
                try:
                    fig_e = plot_gnomad_context(data)
                    st.plotly_chart(fig_e, use_container_width=True)
                except Exception as e:
                    st.error(f"Plot E failed: {e}")
                st.markdown("</div>", unsafe_allow_html=True)
            with col_f:
                st.markdown("<div class='cv-card'>", unsafe_allow_html=True)
                st.markdown("**Plot F — ClinVar Region Map** | Known variants ±50 kb (lollipop)")
                try:
                    fig_f = plot_clinvar_lollipop(data, chrom, pos)
                    st.plotly_chart(fig_f, use_container_width=True)
                except Exception as e:
                    st.error(f"Plot F failed: {e}")
                st.markdown("</div>", unsafe_allow_html=True)

        # ── Tab: Tissue Effects ────────────────────────────────────────────
        with t_tissue:
            col_g, col_h = st.columns([1, 1])
            with col_g:
                st.markdown("<div class='cv-card'>", unsafe_allow_html=True)
                st.markdown("**Plot G — Tissue-Specific Impact** | Sorted by predicted effect magnitude")
                try:
                    te = data.get("tissue_effects", [])
                    fig_g = plot_tissue_heatmap(te)
                    st.plotly_chart(fig_g, use_container_width=True)
                except Exception as e:
                    st.error(f"Plot G failed: {e}")
                st.markdown("</div>", unsafe_allow_html=True)
            with col_h:
                st.markdown("<div class='cv-card'>", unsafe_allow_html=True)
                st.markdown("**Plot H — Enformer Track Heatmap** | Top 20 differentially predicted chromatin/CAGE tracks")
                try:
                    fig_h = plot_enformer_tracks(data)
                    st.plotly_chart(fig_h, use_container_width=True)
                except Exception as e:
                    st.error(f"Plot H failed: {e}")
                st.markdown("</div>", unsafe_allow_html=True)

        # ── Tab: Gene Context ──────────────────────────────────────────────
        with t_gene:
            st.markdown("<div class='cv-card'>", unsafe_allow_html=True)
            g_data = data.get("gene", {})
            if g_data:
                gc1, gc2 = st.columns([2, 1])
                with gc1:
                    st.markdown(f"### {g_data.get('name', metrics.get('gene_symbol', 'N/A'))}")
                    st.markdown(f"**Biotype:** `{g_data.get('biotype', 'N/A')}`")
                    desc = g_data.get("description", "")
                    if desc:
                        st.markdown(f"**Description:** {desc}")
                with gc2:
                    links = g_data.get("links", {})
                    if links:
                        st.markdown("**External Resources**")
                        for name_l, url in links.items():
                            st.markdown(f"[{name_l}]({url})")
                    if st.button("🧬 Open in MSA Explorer", key="btn_to_msa_inf"):
                        st.session_state.page = "MSA Explorer"
                        st.rerun()

                if "expression" in g_data:
                    st.markdown("#### Tissue Expression (GTEx)")
                    expr_df = pd.DataFrame(g_data["expression"])
                    if not expr_df.empty:
                        fig_expr = go.Figure(go.Bar(
                            x=expr_df["tpm"], y=expr_df["tissue"], orientation="h",
                            marker=dict(
                                color=expr_df["tpm"],
                                colorscale=[[0,"#DEE8F5"],[1,COLORS["secondary"]]],
                                showscale=False,
                            ),
                            hovertemplate="%{y}: %{x:.1f} TPM<extra></extra>",
                        ))
                        fig_expr.update_layout(
                            height=max(250, 30*len(expr_df)),
                            xaxis_title="Expression (TPM)",
                            yaxis=dict(autorange="reversed"),
                            margin=dict(l=0, r=0, t=10, b=0),
                            plot_bgcolor="rgba(0,0,0,0)",
                            paper_bgcolor="rgba(0,0,0,0)",
                        )
                        st.plotly_chart(fig_expr, use_container_width=True)
            else:
                st.info("Gene context unavailable (Ensembl API may be offline).")
            st.markdown("</div>", unsafe_allow_html=True)

        # ── Tab: Related Data ──────────────────────────────────────────────
        with t_related:
            st.markdown("<div class='cv-card'>", unsafe_allow_html=True)
            st.markdown("### Known ClinVar / GWAS Associations")
            try:
                payload_rd = {"chrom": chrom, "pos": pos, "ref": ref, "alt": alt}
                r_resp = requests.post(f"{API_URL}/related-data", json=payload_rd, timeout=15)
                if r_resp.status_code == 200:
                    r_data = r_resp.json()
                    if r_data:
                        st.dataframe(pd.DataFrame(r_data), use_container_width=True)
                    else:
                        st.info("No known associations found in ClinVar or GWAS Catalog for this position.")
            except Exception:
                st.warning("Related data API unavailable (API server may be offline).")
            st.markdown("</div>", unsafe_allow_html=True)


# ══════════════════════════════════════════════════════════════════════════════
# Page: Gene Panel
# ══════════════════════════════════════════════════════════════════════════════
GENE_PANEL = [
    {"symbol": "MYH9",  "name": "Myosin Heavy Chain 9",      "disease": "HCM / MYH9-RD",
     "variant": {"chrom": "chr22", "pos": 36191400,   "ref": "A", "alt": "C"}},
    {"symbol": "LMNA",  "name": "Lamin A/C",                  "disease": "DCM / EDMD",
     "variant": {"chrom": "chr1",  "pos": 156104762,  "ref": "G", "alt": "A"}},
    {"symbol": "PCSK9", "name": "Proprotein Convertase 9",    "disease": "Familial HCholesterol",
     "variant": {"chrom": "chr1",  "pos": 55039974,   "ref": "G", "alt": "A"}},
    {"symbol": "ACTN2", "name": "Actinin Alpha 2",            "disease": "HCM / DCM",
     "variant": {"chrom": "chr1",  "pos": 236893000,  "ref": "C", "alt": "T"}},
    {"symbol": "TTN",   "name": "Titin",                      "disease": "DCM / HCM",
     "variant": {"chrom": "chr2",  "pos": 178525989,  "ref": "G", "alt": "A"}},
    {"symbol": "APOB",  "name": "Apolipoprotein B",           "disease": "Familial HCholesterol",
     "variant": {"chrom": "chr2",  "pos": 21009300,   "ref": "C", "alt": "T"}},
    {"symbol": "MYL3",  "name": "Myosin Light Chain 3",       "disease": "HCM",
     "variant": {"chrom": "chr3",  "pos": 46900000,   "ref": "A", "alt": "G"}},
]

if selected_page == "Gene Panel":
    st.markdown("<div class='section-title'>Gene Panel</div>", unsafe_allow_html=True)
    st.markdown("<div class='section-sub'>Pre-loaded pathogenic variants for key cardiovascular disease genes.</div>", unsafe_allow_html=True)

    cols = st.columns(3)
    for i, gene in enumerate(GENE_PANEL):
        v = gene["variant"]
        with cols[i % 3]:
            st.markdown(f"""
            <div class='gene-card'>
                <div style='font-size:20px;font-weight:800;color:{COLORS["primary"]};'>{gene['symbol']}</div>
                <div style='font-size:13px;color:{COLORS["text_dark"]};font-weight:500;margin:4px 0;'>{gene['name']}</div>
                <div style='font-size:11px;color:{COLORS["text_light"]};margin-bottom:10px;'>🫀 {gene['disease']}</div>
                <div style='font-size:11px;font-family:monospace;background:#F4F5F7;padding:6px 10px;border-radius:6px;'>
                    {v['chrom']}:{v['pos']:,} {v['ref']}→{v['alt']}
                </div>
            </div>
            """, unsafe_allow_html=True)
            if st.button(f"Analyze {gene['symbol']}", key=f"gene_btn_{gene['symbol']}",
                         use_container_width=True):
                st.session_state.inputs.update(gene["variant"])
                st.session_state.page = "Inference"
                st.rerun()


# ══════════════════════════════════════════════════════════════════════════════
# Page: Benchmarking
# ══════════════════════════════════════════════════════════════════════════════
if selected_page == "Benchmarking":
    st.markdown("<div class='section-title'>Batch Analysis</div>", unsafe_allow_html=True)
    st.markdown("<div class='section-sub'>Upload a CSV (`chrom,pos,ref,alt`) or VCF file.</div>", unsafe_allow_html=True)

    st.markdown("<div class='cv-card'>", unsafe_allow_html=True)
    bc1, bc2 = st.columns([3, 1])
    with bc1:
        uploaded_file = st.file_uploader("Choose File", type=["csv", "vcf"])
    with bc2:
        st.markdown("<div style='height:28px'></div>", unsafe_allow_html=True)
        if st.button("⚡ Load Demo Batch"):
            demo_csv = "chrom,pos,ref,alt\nchr22,36191400,A,C\nchr1,55039974,G,A\nchr11,47332560,T,C"
            uploaded_file = io.BytesIO(demo_csv.encode("utf-8"))
            uploaded_file.name = "demo_batch.csv"
            st.success("Demo loaded!")

    if uploaded_file:
        if st.button("▶ Process Batch"):
            try:
                if uploaded_file.name.endswith(".vcf"):
                    from vcf_parser import parse_vcf
                    df = parse_vcf(uploaded_file.getvalue().decode("utf-8"))
                    st.info(f"Parsed {len(df)} variants from VCF.")
                else:
                    df = pd.read_csv(uploaded_file)

                required_cols = {"chrom", "pos", "ref", "alt"}
                if not required_cols.issubset(df.columns):
                    st.error(f"Missing columns: {required_cols - set(df.columns)}")
                else:
                    variants = df.to_dict(orient="records")
                    with st.spinner("Starting batch job…"):
                        bres = requests.post(f"{API_URL}/batch-start", json={"variants": variants})
                        bres.raise_for_status()
                        batch_id = bres.json()["batch_id"]

                    pb = st.progress(0)
                    st_txt = st.empty()
                    while True:
                        status = requests.get(f"{API_URL}/batch-status/{batch_id}").json()
                        if status["total"] > 0:
                            pb.progress(status["processed"] / status["total"])
                            st_txt.text(f"Processed {status['processed']}/{status['total']}")
                        if status["status"] in ["completed", "failed"]:
                            break
                        time.sleep(1)

                    if status["status"] == "completed":
                        res_df = pd.DataFrame(status["results"])
                        st.dataframe(res_df, use_container_width=True)
                        st.download_button("⬇ Download CSV", res_df.to_csv(index=False).encode(),
                                           "batch_results.csv", "text/csv")
            except Exception as e:
                st.error(f"Batch error: {e}")
    st.markdown("</div>", unsafe_allow_html=True)


# ══════════════════════════════════════════════════════════════════════════════
# Page: MSA Explorer
# ══════════════════════════════════════════════════════════════════════════════
if selected_page == "MSA Explorer":
    st.markdown("<div class='section-title'>MSA Explorer</div>", unsafe_allow_html=True)
    st.markdown("<div class='section-sub'>Visualize conservation and sequence diversity from FASTA alignments.</div>", unsafe_allow_html=True)

    st.markdown("<div class='cv-card'>", unsafe_allow_html=True)
    mc1, mc2 = st.columns([3, 1])
    with mc1:
        uploaded_msa = st.file_uploader("Upload FASTA Alignment", type=["fasta", "fa", "txt"])
    with mc2:
        st.markdown("<div style='height:28px'></div>", unsafe_allow_html=True)
        if st.button("⚡ Load Demo MSA"):
            demo_fasta = ">seq1\nATGCGTACGTTAG\n>seq2\nATGCGTACGT-AG\n>seq3\nATGCGTACGTTAG\n>seq4\nATGCATACGTTCG"
            uploaded_msa = io.BytesIO(demo_fasta.encode("utf-8"))
            uploaded_msa.name = "demo.fasta"
            st.success("Demo MSA loaded!")

    if uploaded_msa:
        try:
            from msaexplorer.explore import MSA
            import msaexplorer.draw as msa_draw
            import matplotlib.pyplot as plt

            content = uploaded_msa.getvalue().decode("utf-8")
            msa = MSA(content)
            st.success(f"Loaded alignment — {len(msa.alignment)} sequences × {msa.length} positions")

            msa_seq, msa_sim, msa_stat = st.tabs(["Sequence Similarity", "Statistics", "Shannon Entropy"])

            with msa_seq:
                st.markdown("<div class='cv-card'>", unsafe_allow_html=True)
                fig, ax = plt.subplots(figsize=(12, 5))
                msa_draw.similarity_alignment(msa, ax=ax)
                ax.set_title("Sequence Similarity Matrix", fontsize=13, fontweight="bold")
                st.pyplot(fig, use_container_width=True)
                st.markdown("</div>", unsafe_allow_html=True)

            with msa_stat:
                st.markdown("<div class='cv-card'>", unsafe_allow_html=True)
                fig2, ax2 = plt.subplots(figsize=(12, 4))
                msa_draw.stat_plot(msa, stat_type="identity", ax=ax2)
                ax2.set_title("Per-Position Identity", fontsize=13, fontweight="bold")
                st.pyplot(fig2, use_container_width=True)
                st.markdown("</div>", unsafe_allow_html=True)

            with msa_stat:
                st.markdown("<div class='cv-card'>", unsafe_allow_html=True)
                st.markdown("**Sequence Diversity (Shannon Entropy)**")
                st.markdown("")

                entropy_vals = []
                for i in range(msa.length):
                    col_res = [seq.seq[i] for seq in msa.alignment if i < len(seq.seq)]
                    counts  = collections.Counter(col_res)
                    total   = len(col_res)
                    H = -sum((c / total) * math.log2(c / total) for c in counts.values() if c > 0)
                    entropy_vals.append(H)

                max_h = math.log2(4) if entropy_vals else 2  # max for DNA (4 states)
                colors_ent = [
                    f"rgba({int(244*(h/max_h))},{int(96*(1-h/max_h))},{int(54*(1-h/max_h))},0.85)"
                    for h in entropy_vals
                ]

                fig_ent = go.Figure(go.Bar(
                    x=list(range(1, msa.length + 1)),
                    y=entropy_vals,
                    marker_color=colors_ent,
                    name="Shannon Entropy",
                    hovertemplate="Position %{x}: H=%{y:.3f} bits<extra></extra>",
                ))
                fig_ent.add_hline(y=max_h, line_dash="dot", line_color="#888DA7",
                                  annotation_text="Max entropy (4 states)",
                                  annotation_position="top right")
                fig_ent.update_layout(
                    xaxis_title="Alignment Position",
                    yaxis_title="Shannon Entropy (bits)",
                    height=380,
                    paper_bgcolor="rgba(0,0,0,0)",
                    plot_bgcolor="rgba(0,0,0,0)",
                    font=dict(family="Inter, sans-serif"),
                )
                st.plotly_chart(fig_ent, use_container_width=True)
                st.markdown("</div>", unsafe_allow_html=True)

        except Exception as e:
            st.error(f"MSA error: {e}")
    st.markdown("</div>", unsafe_allow_html=True)


# ══════════════════════════════════════════════════════════════════════════════
# Page: Report
# ══════════════════════════════════════════════════════════════════════════════
if selected_page == "Report":
    st.markdown("<div class='section-title'>Variant Report Generator</div>", unsafe_allow_html=True)
    st.markdown("<div class='cv-card'>", unsafe_allow_html=True)

    if st.session_state.analysis_results:
        data    = st.session_state.analysis_results
        metrics = data.get("metrics", {})
        gene    = data.get("gene", {})

        st.success(f"Ready to generate report — {metrics.get('gene_symbol')} {chrom}:{pos} {ref}→{alt}")

        if st.button("📄 Generate HTML Report"):
            try:
                _, label, _ = plot_evidence_stack(data)
            except Exception:
                label = "VUS"

            html = f"""<!DOCTYPE html>
<html><head><meta charset="utf-8">
<title>CardioVar Report — {metrics.get('gene_symbol')}</title>
<style>
  body{{font-family:Inter,sans-serif;padding:40px;max-width:860px;margin:auto;color:#171123;}}
  h1{{color:#F46036;border-bottom:3px solid #F46036;padding-bottom:10px;}}
  h2{{color:#5B85AA;margin-top:28px;}}
  .badge{{display:inline-block;padding:5px 16px;border-radius:16px;font-weight:700;font-size:13px;background:#fee;color:#c44027;border:1.5px solid #c44027;}}
  table{{width:100%;border-collapse:collapse;margin-top:14px;}}
  th,td{{border:1px solid #ddd;padding:8px 12px;text-align:left;font-size:13px;}}
  th{{background:#f7f7f7;font-weight:600;}}
  .footer{{margin-top:40px;font-size:11px;color:#888;border-top:1px solid #eee;padding-top:12px;}}
</style></head><body>
<h1>CardioVar Variant Report</h1>
<p>Generated: {time.strftime("%Y-%m-%d %H:%M")} UTC | Model: {metrics.get("model_used","N/A")}</p>
<span class='badge'>{label}</span>
<h2>Variant Summary</h2>
<table>
<tr><th>Field</th><th>Value</th></tr>
<tr><td>Variant ID</td><td>{chrom}:{pos} {ref}&gt;{alt}</td></tr>
<tr><td>Gene</td><td>{metrics.get("gene_symbol","N/A")}</td></tr>
<tr><td>Assembly</td><td>GRCh38</td></tr>
<tr><td>Max Impact (|Δ|)</td><td>{abs(metrics.get("max_delta",0)):.4f}</td></tr>
<tr><td>Z-Score</td><td>{metrics.get("z_score",0):.2f}</td></tr>
<tr><td>Percentile</td><td>{metrics.get("percentile",0):.1f}%</td></tr>
<tr><td>gnomAD Allele Frequency</td><td>{metrics.get("gnomad_freq",0):.2e}</td></tr>
<tr><td>Model Confidence</td><td>{metrics.get("confidence",0):.0f}%</td></tr>
</table>
<h2>Gene Information</h2>
<p><b>Name:</b> {gene.get("name","N/A")}<br>
<b>Biotype:</b> {gene.get("biotype","N/A")}<br>
<b>Description:</b> {gene.get("description","N/A")}</p>
<h2>Disclaimer</h2>
<p>This report is generated by CardioVar for <strong>research purposes only</strong>. It does not constitute clinical advice. Always consult a certified clinical geneticist for patient-facing interpretation.</p>
<div class='footer'>CardioVar v1.2.0 · Enformer · GRCh38 · Generated {time.strftime("%Y-%m-%d")}</div>
</body></html>"""
            b64 = base64.b64encode(html.encode()).decode()
            href = f'<a href="data:text/html;base64,{b64}" download="cardiovar_report_{metrics.get("gene_symbol","variant")}.html" style="text-decoration:none;background:{COLORS["primary"]};color:white;padding:10px 22px;border-radius:9px;font-weight:600;">⬇ Download HTML Report</a>'
            st.markdown(href, unsafe_allow_html=True)
    else:
        st.info("No analysis results. Please run an analysis in the **Inference** tab first.")
    st.markdown("</div>", unsafe_allow_html=True)


# ══════════════════════════════════════════════════════════════════════════════
# Page: Models / About
# ══════════════════════════════════════════════════════════════════════════════
if selected_page == "Models":
    st.markdown("<div class='section-title'>About CardioVar</div>", unsafe_allow_html=True)
    st.markdown("<div class='section-sub'>Precision genomics platform for cardiovascular disease variant interpretation.</div>", unsafe_allow_html=True)

    ab1, ab2 = st.columns(2)
    with ab1:
        st.markdown("<div class='cv-card'>", unsafe_allow_html=True)
        st.markdown("### 📘 How to Use")
        st.markdown("""
**1. Inference** — Enter a GRCh38 variant or load a demo. Enformer predicts regulatory impact across 5,313 genomic assay tracks.

**2. Gene Panel** — Click any cardio gene to pre-fill a canonical pathogenic variant; analysis starts immediately.

**3. Benchmarking** — Upload a CSV/VCF for batch processing. Results downloadable as CSV.

**4. MSA Explorer** — Upload a FASTA alignment to visualise conservation and per-position Shannon entropy.

**5. Report** — Generate a downloadable HTML report after running inference.
        """)
        st.markdown("</div>", unsafe_allow_html=True)

    with ab2:
        st.markdown("<div class='cv-card'>", unsafe_allow_html=True)
        st.markdown("### 🧠 Models & Methods")
        st.markdown("""
| Component | Method |
|:---|:---|
| **Impact model** | Enformer (Avsec et al. 2021, *Nature Genetics*) |
| **Population freq** | gnomAD v4 (Karczewski et al. 2020) |
| **Conservation** | PhyloP 100-way vertebrate (UCSC) |
| **Z-Score** | Background distribution per gene |
| **Entropy** | Shannon entropy from FASTA MSA |
| **Pathogenicity** | Multi-evidence heuristic |
        """)
        st.markdown("</div>", unsafe_allow_html=True)


# ══════════════════════════════════════════════════════════════════════════════
# Page: Settings
# ══════════════════════════════════════════════════════════════════════════════
if selected_page == "Settings":
    st.markdown("<div class='section-title'>Settings</div>", unsafe_allow_html=True)
    st.markdown("<div class='cv-card'>", unsafe_allow_html=True)
    st.markdown("### Backend Status")
    if st.button("🔗 Check API Connectivity"):
        try:
            resp = requests.get(f"{API_URL}/system-status", timeout=5)
            if resp.status_code == 200:
                st.success("API Online ✓")
                st.json(resp.json())
            else:
                st.error(f"API Error: HTTP {resp.status_code}")
        except Exception as e:
            st.error(f"Connection failed: {e}")
    st.markdown(f"API URL: `{API_URL}`")
    st.markdown("</div>", unsafe_allow_html=True)
