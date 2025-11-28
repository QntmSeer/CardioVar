import warnings
# Suppress numpy binary incompatibility warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

import streamlit as st
import pandas as pd
import numpy as np
import requests
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.colors as mcolors
from matplotlib.patches import Patch
import io

# --- Configuration & Constants ---
import time

# --- Configuration & Constants ---
st.set_page_config(
    page_title="CardioVar",
    page_icon="❤️",
    layout="wide",
    initial_sidebar_state="expanded"
)

import os

# API Configuration
API_URL = os.getenv("API_URL", "http://localhost:8000")

# Color Palette (Tiger Flame / Midnight Violet Theme)
COLORS = {
    "tiger_flame": "#F46036",
    "dusty_denim": "#5B85AA",
    "twilight_indigo": "#414770",
    "midnight_violet": "#171123",
    "dark_amethyst": "#372248",
    "gunmetal": "#303633",
    "aquamarine": "#8BE8CB",
    "cool_steel": "#5B85AA",
    "lavender_grey": "#888DA7",
    "dusty_mauve": "#9C7A97"
}

# Seaborn Theme (Light Mode)
sns.set_theme(style="white", context="notebook")
plt.rcParams['figure.facecolor'] = 'none'
plt.rcParams['axes.facecolor'] = 'none'
plt.rcParams['text.color'] = COLORS['gunmetal']
plt.rcParams['axes.labelcolor'] = COLORS['gunmetal']
plt.rcParams['xtick.color'] = COLORS['gunmetal']
plt.rcParams['ytick.color'] = COLORS['gunmetal']
plt.rcParams['axes.edgecolor'] = COLORS['cool_steel']

# Custom CSS for Fonts & Styling (Light Mode)
st.markdown("""
    <style>
        @import url('https://fonts.googleapis.com/css2?family=Outfit:wght@300;400;600;700&display=swap');
        
        html, body, [class*="css"] {
            font-family: 'Outfit', sans-serif;
        }
        
        h1, h2, h3 {
            font-family: 'Outfit', sans-serif;
            font-weight: 700;
            color: #2C5F7C !important;
        }
        
        .stButton>button {
            background-color: #F46036;
            color: white;
            border-radius: 8px;
            font-weight: 600;
            border: none;
            transition: background-color 0.3s ease;
        }
        .stButton>button:hover {
            background-color: #D94E28; /* Darker shade of Tiger Flame */
            color: white;
            border: none;
        }
        
        /* Light mode adjustments */
        [data-testid="stSidebar"] {
            background-color: #F8F9FA;
        }
    </style>
""", unsafe_allow_html=True)

# --- Helper Functions ---
def plot_deltas_from_api(data, chrom, pos, ref, alt, line_color, highlight_color):
    """
    Plot using data returned from API.
    """
    curve = data["curve"]
    metrics = data["metrics"]
    tracks = data["tracks"]
    
    x = np.array(curve["x"])
    y = np.array(curve["y"])
    
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 9), sharex=True, gridspec_kw={'height_ratios': [3, 0.5, 0.8]})
    
    # 1. Main Delta Plot
    sns.lineplot(x=x, y=y, color=line_color, linewidth=2.5, ax=ax1, label='$\Delta$ RNA-seq')
    ax1.axhline(0, color=COLORS['lavender_grey'], linestyle='--', alpha=0.5)
    ax1.axvline(0, color=highlight_color, linestyle=':', alpha=0.8)
    
    # Highlight Max
    max_pos = metrics["max_pos_rel"]
    max_val = metrics["max_delta"]
    ax1.scatter(max_pos, max_val, color=highlight_color, s=150, zorder=5, edgecolor='white', linewidth=1.5, label='Max Impact')
    ax1.annotate(f'{max_val:.2f}', (max_pos, max_val), xytext=(10, 10), textcoords='offset points',
                 bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=highlight_color, alpha=0.9),
                 arrowprops=dict(arrowstyle='->', connectionstyle="arc3,rad=.2", color=highlight_color))
    
    ax1.set_ylabel("Impact Score")
    ax1.set_title(f"Variant Impact", fontsize=16, fontweight='bold', loc='left', color=COLORS['gunmetal'])
    sns.despine(ax=ax1, left=True, bottom=True)
    ax1.grid(axis='y', linestyle=':', alpha=0.3)
    
    # Legend below plot
    ax1.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=2, frameon=False)
    
    # 2. Gene Structure Track
    exons = tracks.get("exons", [])
    ax2.plot([x[0], x[-1]], [0, 0], color=COLORS['lavender_grey'], linewidth=1)
    
    exon_color = COLORS['cool_steel']
    
    if not exons:
        ax2.text(0, 0, "Non-coding Region", ha='center', va='center', fontsize=10, color=COLORS['lavender_grey'], style='italic')
    else:
        for exon in exons:
            start, end = exon["start"], exon["end"]
            # Clip to window
            start = max(x[0], start)
            end = min(x[-1], end)
            
            if end > start:
                rect = plt.Rectangle((start, -0.4), end-start, 0.8, facecolor=exon_color, alpha=0.8, edgecolor='none')
                ax2.add_patch(rect)
                # Label with Exon ID (shortened)
                label = exon.get("id", "Exon")
                if len(label) > 10: label = "Exon"
                ax2.text((start+end)/2, 0, label, ha='center', va='center', color='white', fontweight='bold', fontsize=8)
    
    ax2.set_yticks([])
    ax2.set_ylabel("Gene")
    sns.despine(ax=ax2, left=True, bottom=True)
    
    cons = np.array(tracks["conservation"])
    cons_pos_color = COLORS['aquamarine']
    cons_neg_color = COLORS['lavender_grey']
    
    ax3.fill_between(x, cons, 0, where=(cons>0), color=cons_pos_color, alpha=0.8, label='Conserved')
    ax3.fill_between(x, cons, 0, where=(cons<0), color=cons_neg_color, alpha=0.4, label='Accelerated')
    ax3.set_ylabel("Conservation")
    ax3.set_xlabel("Distance (bp)")
    sns.despine(ax=ax3, left=True, bottom=False)
    
    # Legend for conservation below
    ax3.legend(loc='upper center', bbox_to_anchor=(0.5, -0.5), ncol=2, frameon=False)
    
    plt.tight_layout()
    return fig

# --- Main UI ---

# Sidebar
with st.sidebar:
    # Minimal text-based logo
    st.markdown(f"""
    <div style="text-align: center; padding: 20px 0; border-bottom: 1px solid #e0e0e0; margin-bottom: 20px;">
        <h1 style="font-family: 'Helvetica Neue', sans-serif; font-weight: 300; font-size: 32px; margin: 0; letter-spacing: -1px; color: {COLORS['gunmetal']};">
            CardioVar
        </h1>
        <p style="font-family: 'Helvetica Neue', sans-serif; font-weight: 300; font-size: 11px; margin: 5px 0 0 0; letter-spacing: 2px; color: {COLORS['lavender_grey']}; text-transform: uppercase;">
            Variant Explorer
        </p>
    </div>
    """, unsafe_allow_html=True)
    st.header("Variant Configuration")
    
    # Genome Build Display (Static)
    st.markdown("**Genome Build:** GRCh38 / hg38")
    assembly_code = "GRCh38"
    st.caption("ℹ️ Only GRCh38 is currently supported.")
    
    chrom = st.selectbox("Chromosome", ["chr22", "chr1", "chr2", "chr3"])
    position = st.number_input("Position", min_value=0, value=36191400, step=100)
    col1, col2 = st.columns(2)
    with col1:
        ref = st.text_input("Reference", "A")
    with col2:
        alt = st.text_input("Alternate", "C")
        
    st.markdown("---")
    st.header("🎨 Customization")
    
    line_color = st.color_picker("Signal Color", COLORS['cool_steel'])
    highlight_color = st.color_picker("Highlight Color", COLORS['dusty_mauve'])
    
    st.markdown("---")
    st.header("⚙️ Advanced")
    force_live = st.checkbox("Force live API calls (debug)", value=False, help="Bypass local cache and force real API requests. May be slower.")
    
    st.markdown("---")
    st.header("🖥️ System Health")
    if st.button("Check Status"):
        try:
            sys_resp = requests.get(f"{API_URL}/system-status")
            if sys_resp.status_code == 200:
                sys_data = sys_resp.json()
                st.metric("CPU Usage", f"{sys_data['cpu_percent']}%")
                st.metric("RAM Usage", f"{sys_data['memory_percent']}%")
                st.caption(f"Used: {sys_data['memory_used_gb']} GB / {sys_data['memory_total_gb']} GB")
            else:
                st.error("Failed to fetch system status")
        except Exception as e:
            st.error(f"Connection error: {e}")
    
    # Model info moved to About tab
    
    run_btn = st.button("Run Analysis", type="primary", use_container_width=True)

# Main Page Tabs
tab1, tab2, tab3, tab4, tab5 = st.tabs(["Variant Explorer", "Annotations", "Related Data", "Batch Analysis", "About"])
    
if run_btn:
    # Call API
    payload = {
        "assembly": assembly_code, 
        "chrom": chrom, 
        "pos": position, 
        "ref": ref, 
        "alt": alt,
        "force_live": force_live
    }
    
    try:
        with st.spinner("Querying Variant Engine..."):
            resp = requests.post(f"{API_URL}/variant-impact", json=payload)
            resp.raise_for_status()
            data = resp.json()
            g_data = data.get('gene', {})

        # Data Provenance Expander in Variant Explorer (tab1)
        with tab1:
            if 'data_sources' in data:
                with st.expander("📊 Data Sources & Provenance"):
                    st.markdown("**Data Provenance**")
                    df_sources = pd.DataFrame(list(data['data_sources'].items()), columns=["Data Element", "Source"])
                    st.table(df_sources)
            
            # Plot the variant impact (Delta Plot)
            # We need to ensure the plot is rendered in tab1 as well, 
            # assuming the original code had it there. 
            # Looking at previous context, plot_deltas_from_api was used.
            if 'curve' in data and 'metrics' in data:
                fig = plot_deltas_from_api(data, chrom, position, ref, alt, line_color, highlight_color)
                st.pyplot(fig)

        # Annotations tab (tab2)
        with tab2:
            # Cached data warning
            if data.get('fallback_used', False):
                st.warning("ℹ️ Some annotation data loaded from cache (gene expression/protein domains may be cached). Core variant data is live.")

            # Expression Plot
            if 'expression' in g_data:
                st.markdown("### Baseline Expression")
                expr_df = pd.DataFrame(g_data['expression'])
                expr_df['is_cardio'] = expr_df['tissue'].apply(lambda x: any(t in x for t in ['Heart', 'Aorta', 'Coronary']))
                fig_expr, ax_expr = plt.subplots(figsize=(10, 5))
                colors = [COLORS['dusty_mauve'] if c else COLORS['cool_steel'] for c in expr_df['is_cardio']]
                sns.barplot(data=expr_df, x='tissue', y='tpm', palette=colors, ax=ax_expr, edgecolor='none', alpha=0.9)
                ax_expr.set_ylabel('TPM')
                ax_expr.set_xlabel('')
                ax_expr.set_title('Tissue Expression Profile', loc='left', color=COLORS['gunmetal'])
                ax_expr.tick_params(axis='x', rotation=45)
                sns.despine(left=True, bottom=True)
                ax_expr.grid(axis='y', linestyle=':', alpha=0.3)
                legend_elements = [Patch(facecolor=COLORS['dusty_mauve'], label='Cardiovascular'),
                                   Patch(facecolor=COLORS['cool_steel'], label='Other')]
                ax_expr.legend(handles=legend_elements, loc='upper center', bbox_to_anchor=(0.5, -0.3), ncol=2, frameon=False)
                plt.tight_layout()
                st.pyplot(fig_expr)

            # Protein Domains
            if 'protein_domains' in g_data and 'protein_length' in g_data:
                st.markdown("### Protein Architecture")
                domains = g_data['protein_domains']
                prot_len = g_data['protein_length']
                fig_prot, ax_prot = plt.subplots(figsize=(10, 3))
                ax_prot.plot([0, prot_len], [0, 0], color=COLORS['lavender_grey'], linewidth=4, solid_capstyle='round')
                domain_colors = [COLORS['aquamarine'], COLORS['cool_steel'], COLORS['dusty_mauve']]
                for i, domain in enumerate(domains):
                    color = domain_colors[i % len(domain_colors)]
                    rect = plt.Rectangle((domain['start'], -0.3), domain['end']-domain['start'], 0.6,
                                         facecolor=color, edgecolor='none', alpha=0.9)
                    ax_prot.add_patch(rect)
                    if (domain['end'] - domain['start']) > (prot_len * 0.05):
                        ax_prot.text((domain['start']+domain['end'])/2, 0, domain['name'],
                                     ha='center', va='center', color='white', fontsize=8, fontweight='bold')
                sns.despine(left=True, bottom=True)
                ax_prot.set_yticks([])
                ax_prot.set_xlabel("Amino Acid Position")
                legend_elements = [Patch(facecolor=COLORS['aquamarine'], label='Domain Type A'),
                                   Patch(facecolor=COLORS['cool_steel'], label='Domain Type B'),
                                   Patch(facecolor=COLORS['dusty_mauve'], label='Domain Type C')]
                ax_prot.legend(handles=legend_elements, loc='upper center', bbox_to_anchor=(0.5, -0.5), ncol=3, frameon=False, title="Domain Types (Representative)")
                plt.tight_layout()
                st.pyplot(fig_prot)

            # Pathways, Disease Associations, External Links
            if 'pathways' in g_data:
                st.markdown("### Pathways")
                for p in g_data.get('pathways', []):
                    st.markdown(f"- {p}")
            if 'disease_associations' in g_data:
                st.markdown("### Disease Associations")
                for d in g_data.get('disease_associations', []):
                    st.markdown(f"- {d}")
            if 'links' in g_data:
                st.markdown("### External Links")
                links = g_data.get('links', {})
                if links:
                    for name, url in links.items():
                        st.markdown(f"- [{name}]({url})")
                else:
                    st.info("No external links available.")

        # Tab 3: Related Data
        with tab3:
            st.subheader("Known Associations")
            try:
                r_resp = requests.post(f"{API_URL}/related-data", json=payload)
                if r_resp.status_code == 200:
                    r_data = r_resp.json()
                    if r_data:
                        st.dataframe(pd.DataFrame(r_data), use_container_width=True)
                    else:
                        st.info("No known ClinVar or GWAS associations found for this specific variant.")
                else:
                    st.error("Failed to fetch related data.")
            except Exception as e:
                st.error(f"Error fetching related data: {e}")

    except requests.exceptions.ConnectionError:
        st.error("❌ Could not connect to backend API. Is `api.py` running?")
    except requests.exceptions.HTTPError as e:
        if e.response.status_code == 400:
            error_detail = e.response.json().get("detail", "Unknown error")
            st.warning(f"⚠️ {error_detail}")
        else:
            st.error(f"❌ Server error: {e}")
    except Exception as e:
        st.error(f"An error occurred: {e}")

# Tab 4: Batch Analysis (Independent)
with tab4:
    st.header("Batch Variant Analysis")
    st.markdown("Upload a CSV with columns: `chrom`, `pos`, `ref`, `alt`")
    
    uploaded_file = st.file_uploader("Choose a CSV file", type="csv")
    
    if uploaded_file is not None:
        if st.button("Run Batch Prediction"):
            try:
                df = pd.read_csv(uploaded_file)
                variants = df.to_dict(orient="records")
                
                # Validate columns
                if not all(k in variants[0] for k in ["chrom", "pos", "ref", "alt"]):
                    st.error("CSV must contain chrom, pos, ref, alt columns.")
                else:
                    # Start Batch
                    with st.spinner("Initializing batch..."):
                        start_resp = requests.post(f"{API_URL}/batch-start", json={"variants": variants})
                        start_resp.raise_for_status()
                        batch_data = start_resp.json()
                        batch_id = batch_data["batch_id"]
                    
                    st.success(f"Batch started! ID: {batch_id}")
                    
                    # Progress Bar
                    progress_bar = st.progress(0)
                    status_text = st.empty()
                    
                    # Poll for completion
                    while True:
                        status_resp = requests.get(f"{API_URL}/batch-status/{batch_id}")
                        if status_resp.status_code != 200:
                            st.error("Failed to get batch status.")
                            break
                        
                        status_data = status_resp.json()
                        state = status_data["status"]
                        processed = status_data["processed"]
                        total = status_data["total"]
                        
                        # Update progress
                        if total > 0:
                            progress = min(processed / total, 1.0)
                            progress_bar.progress(progress)
                            status_text.text(f"Processed {processed}/{total} variants...")
                        
                        if state in ["completed", "failed"]:
                            break
                        
                        time.sleep(0.5)
                    
                    if state == "completed":
                        results = status_data["results"]
                        res_df = pd.DataFrame(results)
                        st.success(f"Processed {len(results)} variants.")
                        
                        # Display Table
                        st.dataframe(res_df.style.applymap(lambda x: "color: red" if x == "High" else "color: orange" if x == "Medium" else "color: green", subset=["priority"]), use_container_width=True)
                        
                        # Download
                        csv = res_df.to_csv(index=False).encode('utf-8')
                        st.download_button("💾 Download Batch Results", data=csv, file_name="batch_results.csv", mime="text/csv")
                    else:
                        st.error(f"Batch failed: {status_data.get('error', 'Unknown error')}")
                    
            except Exception as e:
                st.error(f"Batch processing failed: {e}")

# Tab 5: About
with tab5:
    st.header("About CardioVar")
    
    st.markdown("""
    **CardioVar** is a precision medicine tool for predicting the functional impact of genetic variants 
    in cardiovascular disease genes.
    """)
    
    st.subheader("🧠 Deep Learning Model")
    st.markdown("""
    **Enformer** - Sequence-based Transformer model
    
    - **Architecture:** Transformer-based neural network
    - **Context Window:** 196,608 bp genomic sequence
    - **Output Tracks:** 5,313 epigenomic features
    - **Hardware:** GPU-accelerated (CUDA)
    - **Reference:** [Avsec et al. 2021, Nature Genetics](https://www.nature.com/articles/s41588-021-00782-6)
    
    The model predicts functional impact by analyzing chromatin accessibility, 
    histone modifications, and transcription factor binding patterns across the genome.
    """)
    
    st.subheader("📊 Data Sources")
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("""
        **Live APIs:**
        - ✅ ClinVar (NCBI E-utilities)
        - ✅ gnomAD v4 (population frequencies)
        - ✅ Ensembl REST API (gene annotations)
        - ✅ GTEx v8 (tissue expression)
        - ✅ UCSC (PhyloP conservation)
        """)
    
    with col2:
        st.markdown("""
        **Fallback Data:**
        - 📁 Gene annotations (JSON)
        - 📁 Protein domains (JSON)
        - 📁 Expression profiles (TSV)
        - 📁 Conservation scores (NumPy)
        - 📁 Gene structure (JSON)
        """)
    
    st.subheader("⚙️ Technical Details")
    st.markdown(f"""
    - **Genome Build:** GRCh38 / hg38
    - **API Caching:** 24-hour TTL
    - **Rate Limiting:** 3 calls/second (NCBI)
    - **Model:** {metrics.get('model_used', 'Enformer (Deep Learning)') if 'metrics' in locals() else 'Enformer (Deep Learning)'}
    """)
    
    st.subheader("⚠️ Disclaimer")
    st.warning("""
    **RESEARCH USE ONLY** — This tool is for research and educational purposes only. 
    It is NOT intended for clinical diagnosis or treatment decisions. Always consult 
    with qualified healthcare professionals for medical advice.
    """)
    
    st.subheader("📄 License")
    st.markdown("""
    Apache License 2.0 - Open Source
    
    © 2024 CardioVar Project
    """)
