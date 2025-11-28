import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import io
import json
import os

# Import the variant engine directly (no API needed)
from variant_engine import compute_variant_impact

# Load data files
def load_gene_annotations():
    with open("data/gene_annotations.json", "r") as f:
        return json.load(f)

def load_related_variants():
    with open("data/related_variants.json", "r") as f:
        return json.load(f)

GENE_DATA = load_gene_annotations()
RELATED_DATA = load_related_variants()

# Page Config
st.set_page_config(
    page_title="CardioVar Explorer",
    page_icon="favicon.svg",
    layout="wide"
)

# --- Helper Functions ---
def get_gene_annotation(gene_symbol):
    for record in GENE_DATA:
        if record["symbol"].upper() == gene_symbol.upper():
            return record
    return {"symbol": gene_symbol, "note": "No detailed annotations found."}

def get_related_data(chrom, pos, ref, alt):
    key = f"{chrom}:{pos}:{ref}:{alt}"
    return RELATED_DATA.get(key, [])

def plot_deltas_from_data(data, chrom, pos, ref, alt, line_color, highlight_color):
    """Plot using data returned from variant engine."""
    curve = data["curve"]
    metrics = data["metrics"]
    tracks = data["tracks"]
    
    x = np.array(curve["x"])
    y = np.array(curve["y"])
    
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 8), sharex=True, gridspec_kw={'height_ratios': [3, 0.5, 0.5]})
    
    # 1. Main Delta Plot
    sns.lineplot(x=x, y=y, color=line_color, linewidth=2.5, ax=ax1, label='$\\Delta$ RNA-seq')
    ax1.axhline(0, color='gray', linestyle='--', alpha=0.3)
    ax1.axvline(0, color=highlight_color, linestyle=':', alpha=0.8)
    
    # Highlight Max
    max_pos = metrics["max_pos_rel"]
    max_val = metrics["max_delta"]
    ax1.scatter(max_pos, max_val, color=highlight_color, s=150, zorder=5, edgecolor='white', linewidth=1.5)
    ax1.annotate(f'Max: {max_val:.2f}', (max_pos, max_val), xytext=(10, 10), textcoords='offset points',
                 bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=highlight_color, alpha=0.9),
                 arrowprops=dict(arrowstyle='->', connectionstyle="arc3,rad=.2", color=highlight_color))
    
    ax1.set_ylabel("$\\Delta$ RNA-seq")
    ax1.set_title(f"Variant Impact: {chrom}:{pos} {ref}→{alt}", fontsize=14, fontweight='bold')
    ax1.legend(loc='upper right', frameon=False)
    sns.despine(ax=ax1, trim=True)
    
    # 2. Gene Structure Track
    ax2.plot([x[0], x[-1]], [0, 0], color='black', linewidth=1)
    for exon in tracks["exons"]:
        rect = plt.Rectangle((exon["start"], -0.4), exon["end"]-exon["start"], 0.8, facecolor='#3498db', alpha=0.7)
        ax2.add_patch(rect)
        ax2.text((exon["start"]+exon["end"])/2, 0, exon["label"], ha='center', va='center', color='white', fontsize=8)
    ax2.set_yticks([])
    ax2.set_ylabel("Gene")
    sns.despine(ax=ax2, left=True, bottom=True)
    
    # 3. Conservation Track
    cons = np.array(tracks["conservation"])
    ax3.fill_between(x, cons, 0, where=(cons>0), color='#27ae60', alpha=0.6)
    ax3.fill_between(x, cons, 0, where=(cons<0), color='#95a5a6', alpha=0.3)
    ax3.set_ylabel("PhyloP")
    ax3.set_xlabel("Relative Genomic Coordinate (bp)")
    sns.despine(ax=ax3, bottom=False)
    
    plt.tight_layout()
    return fig

# --- Main UI ---

# Sidebar
with st.sidebar:
    # Minimal text-based logo
    st.markdown("""
    <div style="text-align: center; padding: 20px 0; border-bottom: 1px solid #e0e0e0; margin-bottom: 20px;">
        <h1 style="font-family: 'Helvetica Neue', sans-serif; font-weight: 300; font-size: 32px; margin: 0; letter-spacing: -1px; color: #000;">
            CardioVar
        </h1>
        <p style="font-family: 'Helvetica Neue', sans-serif; font-weight: 300; font-size: 11px; margin: 5px 0 0 0; letter-spacing: 2px; color: #666; text-transform: uppercase;">
            Variant Explorer
        </p>
    </div>
    """, unsafe_allow_html=True)
    st.header("Variant Configuration")
    
    # Genome Build Selector
    assembly = st.selectbox(
        "Genome Build",
        ["GRCh38 / hg38", "GRCh37 / hg19"],
        index=0,
        help="All positions are interpreted in the selected genome build"
    )
    # Extract just the assembly code
    assembly_code = "GRCh38" if "38" in assembly else "GRCh37"
    st.caption("ℹ️ All positions are interpreted in the selected genome build")
    
    chrom = st.selectbox("Chromosome", ["chr22", "chr1", "chr2", "chr3"])
    position = st.number_input("Position", min_value=0, value=36191400, step=100)
    col1, col2 = st.columns(2)
    with col1:
        ref = st.text_input("Reference", "A")
    with col2:
        alt = st.text_input("Alternate", "C")
        
    st.markdown("---")
    st.header("🎨 Customization")
    line_color = st.color_picker("Signal Color", "#4C72B0")
    highlight_color = st.color_picker("Highlight Color", "#DD8452")
    
    run_btn = st.button("Run Analysis", type="primary")

# Title
st.title("🧬 CardioVar: CVD Variant Impact Explorer")

# Intro
st.markdown("""
**CardioVar** is a precision medicine tool designed to predict the functional impact of genetic variants in cardiovascular disease genes. 
By leveraging deep learning models, it estimates how a specific variant alters RNA-seq expression profiles, potentially disrupting gene regulation.

**Key Features:**
- 📉 **Predict Impact**: Visualize $\\Delta$ RNA-seq effects of single nucleotide variants.
- 🧬 **Genomic Context**: View gene structure (exons) and evolutionary conservation.
- 🏥 **Clinical Relevance**: Cross-reference with ClinVar and GWAS Catalog.
- 🚀 **Batch Processing**: Analyze multiple variants simultaneously.
""")

# Tabs
tab1, tab2, tab3, tab4 = st.tabs(["📊 Variant Explorer", "🧬 Gene Annotations", "🗂️ Related Data", "🚀 Batch Analysis"])

if run_btn:
    try:
        with st.spinner("Running Variant Analysis..."):
            # Call variant engine directly with assembly
            data = compute_variant_impact(chrom, position, ref, alt, assembly_code)
            
        metrics = data["metrics"]
        
        # Tab 1: Explorer
        with tab1:
            # Metrics - use containers to prevent jittering
            st.markdown("---")
            c1, c2, c3 = st.columns(3)
            with c1:
                st.metric("Max Delta", f"{metrics['max_delta']:.2f}")
            with c2:
                st.metric("gnomAD Freq", f"{metrics['gnomad_freq']:.5f}")
            with c3:
                st.metric("Gene", metrics['gene_symbol'])
            st.markdown("---")
            
            # Plot
            fig = plot_deltas_from_data(data, chrom, position, ref, alt, line_color, highlight_color)
            st.pyplot(fig)
            
            # Additional Plots Row
            col_a, col_b = st.columns(2)
            
            with col_a:
                st.subheader("Tissue-Specific Impact")
                tissue_df = pd.DataFrame(data["tissue_effects"])
                tissue_df['is_cardio'] = tissue_df['tissue'].apply(lambda x: 'Cardiovascular' if any(t in x for t in ['Heart', 'Aorta', 'Coronary']) else 'Other')
                
                fig_tissue, ax_tissue = plt.subplots(figsize=(6, 4))
                colors = ['#E74C3C' if t == 'Cardiovascular' else '#95A5A6' for t in tissue_df['is_cardio']]
                ax_tissue.barh(tissue_df['tissue'], tissue_df['delta'], color=colors)
                ax_tissue.set_xlabel('|Δ RNA-seq|')
                ax_tissue.set_title('Predicted Impact Across Tissues')
                sns.despine()
                plt.tight_layout()
                st.pyplot(fig_tissue)
                st.caption("🔴 Cardiovascular tissues highlighted")
            
            with col_b:
                st.subheader("Variant Percentile")
                bg_data = data["background_distribution"]
                bg_deltas = np.abs(bg_data["background_deltas"])
                var_delta = bg_data["variant_delta"]
                percentile = metrics['percentile']
                
                fig_dist, ax_dist = plt.subplots(figsize=(6, 4))
                ax_dist.hist(bg_deltas, bins=30, color='#95A5A6', alpha=0.6, edgecolor='black')
                ax_dist.axvline(var_delta, color='#E74C3C', linewidth=3, linestyle='--', label=f'This variant (Top {100-percentile:.1f}%)')
                ax_dist.set_xlabel('|Δ RNA-seq|')
                ax_dist.set_ylabel('Frequency')
                ax_dist.set_title(f'Distribution in {metrics["gene_symbol"]}')
                ax_dist.legend()
                sns.despine()
                plt.tight_layout()
                st.pyplot(fig_dist)
                st.caption(f"📊 This variant is in the **top {100-percentile:.1f}%** of predicted impact")
            
            # Export Plot
            fn = f"cardiovar_plot_{chrom}_{position}.png"
            img = io.BytesIO()
            fig.savefig(img, format='png')
            st.download_button("📸 Download Plot", data=img, file_name=fn, mime="image/png")
            
            # Export Data
            csv_data = pd.DataFrame({
                "x": data["curve"]["x"],
                "y": data["curve"]["y"]
            }).to_csv(index=False).encode('utf-8')
            st.download_button("💾 Download Data (CSV)", data=csv_data, file_name=f"variant_{chrom}_{position}.csv", mime="text/csv")

        # Tab 2: Gene Annotations
        with tab2:
            gene_sym = metrics['gene_symbol']
            st.subheader(f"Annotations for {gene_sym}")
            
            g_data = get_gene_annotation(gene_sym)
            st.write(f"**Name:** {g_data.get('name', 'N/A')}")
            
            # Expression Plot
            if 'expression' in g_data:
                st.markdown("### Baseline Expression Across Tissues")
                expr_df = pd.DataFrame(g_data['expression'])
                expr_df['is_cardio'] = expr_df['tissue'].apply(lambda x: any(t in x for t in ['Heart', 'Aorta', 'Coronary']))
                
                fig_expr, ax_expr = plt.subplots(figsize=(10, 5))
                colors = ['#E74C3C' if c else '#4C72B0' for c in expr_df['is_cardio']]
                ax_expr.bar(expr_df['tissue'], expr_df['tpm'], color=colors, edgecolor='black', alpha=0.8)
                ax_expr.set_ylabel('TPM (Transcripts Per Million)')
                ax_expr.set_title(f'{gene_sym} Expression (GTEx-style)')
                ax_expr.tick_params(axis='x', rotation=45)
                sns.despine()
                plt.tight_layout()
                st.pyplot(fig_expr)
                st.caption("🔴 Cardiovascular tissues | 🔵 Other tissues")
            
            # Protein Domains
            if 'protein_domains' in g_data and 'protein_length' in g_data:
                st.markdown("### Protein Domain Architecture")
                domains = g_data['protein_domains']
                prot_len = g_data['protein_length']
                
                fig_prot, ax_prot = plt.subplots(figsize=(10, 2))
                ax_prot.plot([0, prot_len], [0, 0], color='black', linewidth=2)
                
                domain_colors = ['#3498db', '#9b59b6', '#e67e22', '#1abc9c', '#34495e']
                for i, domain in enumerate(domains):
                    color = domain_colors[i % len(domain_colors)]
                    rect = plt.Rectangle((domain['start'], -0.3), domain['end']-domain['start'], 0.6, 
                                        facecolor=color, edgecolor='black', alpha=0.7)
                    ax_prot.add_patch(rect)
                    ax_prot.text((domain['start']+domain['end'])/2, 0, domain['name'], 
                               ha='center', va='center', fontsize=8, color='white', fontweight='bold')
                
                var_aa_pos = int(prot_len * 0.4)
                ax_prot.plot([var_aa_pos, var_aa_pos], [0.6, 1.2], color='#E74C3C', linewidth=2)
                ax_prot.scatter([var_aa_pos], [1.2], color='#E74C3C', s=100, zorder=5, marker='v')
                ax_prot.text(var_aa_pos, 1.4, 'Variant', ha='center', fontsize=9, color='#E74C3C', fontweight='bold')
                
                ax_prot.set_xlim(0, prot_len)
                ax_prot.set_ylim(-0.5, 1.6)
                ax_prot.set_yticks([])
                ax_prot.set_xlabel('Amino Acid Position')
                ax_prot.set_title(f'{gene_sym} Protein Domains (Length: {prot_len} aa)')
                sns.despine(left=True)
                plt.tight_layout()
                st.pyplot(fig_prot)
                st.caption("🔻 Approximate variant position (mock)")
            
            st.markdown("### Pathways")
            for p in g_data.get("pathways", []):
                st.markdown(f"- {p}")
                
            st.markdown("### Disease Associations")
            for d in g_data.get("disease_associations", []):
                st.markdown(f"- {d}")
                
            st.markdown("### External Links")
            links = g_data.get("links", {})
            for k, v in links.items():
                st.markdown(f"[{k}]({v})")

        # Tab 3: Related Data
        with tab3:
            st.subheader("Known Associations")
            r_data = get_related_data(chrom, position, ref, alt)
            if r_data:
                st.dataframe(pd.DataFrame(r_data), use_container_width=True)
            else:
                st.info("No known ClinVar or GWAS associations found for this specific variant.")
                
    except ValueError as e:
        # Handle assembly validation errors
        st.warning(f"⚠️ {str(e)}")
    except Exception as e:
        st.error(f"An error occurred: {e}")

# Tab 4: Batch Analysis
with tab4:
    st.header("Batch Variant Analysis")
    st.markdown("Upload a CSV with columns: `chrom`, `pos`, `ref`, `alt`")
    
    uploaded_file = st.file_uploader("Choose a CSV file", type="csv")
    
    if uploaded_file is not None:
        if st.button("Run Batch Prediction"):
            try:
                df = pd.read_csv(uploaded_file)
                
                if not all(k in df.columns for k in ["chrom", "pos", "ref", "alt"]):
                    st.error("CSV must contain chrom, pos, ref, alt columns.")
                else:
                    results = []
                    progress_bar = st.progress(0)
                    
                    for idx, row in df.iterrows():
                        res = compute_variant_impact(row['chrom'], row['pos'], row['ref'], row['alt'])
                        metrics = res["metrics"]
                        
                        priority = "Low"
                        if abs(metrics["max_delta"]) > 3.0:
                            priority = "High"
                        elif abs(metrics["max_delta"]) > 1.5:
                            priority = "Medium"
                            
                        results.append({
                            "variant_id": res["variant_id"],
                            "gene": metrics["gene_symbol"],
                            "max_delta": metrics["max_delta"],
                            "gnomad_freq": metrics["gnomad_freq"],
                            "priority": priority
                        })
                        
                        progress_bar.progress((idx + 1) / len(df))
                    
                    res_df = pd.DataFrame(results)
                    st.success(f"Processed {len(results)} variants.")
                    
                    st.dataframe(res_df, use_container_width=True)
                    
                    csv = res_df.to_csv(index=False).encode('utf-8')
                    st.download_button("💾 Download Batch Results", data=csv, file_name="batch_results.csv", mime="text/csv")
                    
            except Exception as e:
                st.error(f"Batch processing failed: {e}")
