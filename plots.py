"""
plots.py — CardioVar Production Visualization Module
All 8 scientific-grade interactive Plotly figures.
"""
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import scipy.stats as stats

# ── Colour palette (matches dashboard.py COLORS) ──────────────────────────────
C = {
    "primary":     "#F46036",
    "secondary":   "#5B85AA",
    "accent":      "#9C7A97",
    "text_dark":   "#171123",
    "text_light":  "#888DA7",
    "success":     "#2A9D8F",
    "warning":     "#E9C46A",
    "danger":      "#C44027",
    "bg":          "rgba(0,0,0,0)",
}

_LAYOUT = dict(
    plot_bgcolor=C["bg"],
    paper_bgcolor=C["bg"],
    font=dict(family="Inter, sans-serif", size=13, color=C["text_dark"]),
    margin=dict(l=20, r=20, t=50, b=20),
    hovermode="x unified",
)


def _unavailable_figure(title: str, message: str, height: int = 360) -> go.Figure:
    """Return a clean placeholder figure when real data is unavailable."""
    fig = go.Figure()
    fig.add_annotation(
        text=message,
        x=0.5,
        y=0.5,
        xref="paper",
        yref="paper",
        showarrow=False,
        font=dict(size=14, color=C["text_light"]),
    )
    fig.update_layout(
        title=title,
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
        height=height,
        **_LAYOUT,
    )
    return fig


# ── Plot A: Multi-Track Genomic Browser ───────────────────────────────────────
def plot_variant_impact_profile(data: dict) -> go.Figure:
    """
    3-track figure: (1) ΔRNA-seq impact, (2) Gene structure, (3) PhyloP conservation.
    Draws ±1 SD ribbon only when real pre-computed background exists, and marks PhyloP
    as unavailable if UCSC lookup failed.
    """
    curve   = data.get("curve", {})
    metrics = data.get("metrics", {})
    tracks  = data.get("tracks", {})
    bg      = data.get("background_distribution", {})

    x     = np.array(curve.get("x", []))
    y     = np.array(curve.get("y", []))
    raw_cons = tracks.get("conservation")
    cons  = np.array(raw_cons) if raw_cons is not None else None
    exons = tracks.get("exons", [])

    bg_deltas = bg.get("background_deltas")

    fig = make_subplots(
        rows=3, cols=1,
        shared_xaxes=True,
        vertical_spacing=0.07,
        row_heights=[0.55, 0.20, 0.25],
        subplot_titles=(
            "Enformer Variant Impact (ΔRNA-seq signal)",
            "Gene Structure",
            "Conservation Score (PhyloP, 100-way)",
        ),
    )

    # — Track 1: optional real background ribbon then signal —————
    if bg_deltas and len(bg_deltas) > 1 and len(x):
        bg_mean   = float(np.mean(bg_deltas))
        bg_std    = float(np.std(bg_deltas))
        ribbon_hi = np.full_like(y, bg_mean + bg_std)
        ribbon_lo = np.full_like(y, -(bg_mean + bg_std))
        fig.add_trace(go.Scatter(
            x=np.concatenate([x, x[::-1]]),
            y=np.concatenate([ribbon_hi, ribbon_lo[::-1]]),
            fill="toself", fillcolor="rgba(91,133,170,0.12)",
            line=dict(width=0), name="±1 SD background",
            hoverinfo="skip", showlegend=True,
        ), row=1, col=1)

    fig.add_trace(go.Scatter(
        x=x, y=y, mode="lines",
        name="ΔRNA-seq",
        line=dict(color=C["secondary"], width=2),
        hovertemplate="Position: %{x} bp<br>Impact: %{y:.4f}<extra></extra>",
    ), row=1, col=1)

    max_pos = metrics.get("max_pos_rel", 0)
    max_val = metrics.get("max_delta", 0)
    fig.add_trace(go.Scatter(
        x=[max_pos], y=[max_val], mode="markers",
        name="Peak impact",
        marker=dict(color=C["primary"], size=12, symbol="diamond",
                    line=dict(color="white", width=2)),
        hovertemplate=f"Peak: {max_val:.4f} at {max_pos} bp<extra></extra>",
    ), row=1, col=1)

    fig.add_hline(y=0, line_dash="dash", line_color=C["text_light"], opacity=0.6, row=1, col=1)
    fig.add_vline(x=0, line_dash="dot", line_color=C["primary"], opacity=0.8,
                  annotation_text="Variant site", annotation_position="top right", row=1, col=1)

    # — Track 2: gene structure ———————————————————————————————
    if exons:
        for ex in exons:
            s = max(x[0] if len(x) else -100, ex["start"])
            e = min(x[-1] if len(x) else 100, ex["end"])
            if e > s:
                fig.add_shape(type="rect",
                    x0=s, y0=-0.4, x1=e, y1=0.4,
                    fillcolor=C["secondary"], opacity=0.85,
                    line=dict(color=C["secondary"]), row=2, col=1)
    else:
        fig.add_annotation(text="No exon annotations returned in window",
                           x=0, y=0, showarrow=False,
                           font=dict(color=C["text_light"]), row=2, col=1)

    fig.update_yaxes(showticklabels=False, range=[-1, 1], row=2, col=1)

    # — Track 3: PhyloP conservation ————————————————————————————
    if cons is not None and len(x) == len(cons) and len(cons) > 0:
        pos_mask = cons >= 0
        neg_mask = cons < 0
        fig.add_trace(go.Bar(
            x=x[pos_mask], y=cons[pos_mask],
            name="Conserved (+)", marker_color=C["success"], opacity=0.8,
            hovertemplate="PhyloP: %{y:.2f}<extra></extra>",
        ), row=3, col=1)
        fig.add_trace(go.Bar(
            x=x[neg_mask], y=cons[neg_mask],
            name="Fast-evolving (−)", marker_color=C["danger"], opacity=0.8,
            hovertemplate="PhyloP: %{y:.2f}<extra></extra>",
        ), row=3, col=1)
        fig.add_hline(y=0, line_dash="solid", line_color=C["text_light"], opacity=0.5, row=3, col=1)
    else:
        fig.add_annotation(text="PhyloP conservation data unavailable (UCSC API failed)",
                           x=0, y=0, showarrow=False,
                           font=dict(color=C["text_light"]), row=3, col=1)
        fig.update_yaxes(showticklabels=False, range=[-1, 1], row=3, col=1)

    fig.update_layout(
        height=720,
        barmode="relative",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        **_LAYOUT,
    )
    fig.update_yaxes(title_text="Impact (Δ)", row=1, col=1)
    fig.update_yaxes(title_text="PhyloP", row=3, col=1)
    fig.update_xaxes(title_text="Distance from variant site (bp)", row=3, col=1)
    return fig


# ── Plot B: Background Distribution KDE ───────────────────────────────────────
def plot_background_kde(data: dict) -> go.Figure:
    """KDE of gene background distribution with the variant's delta marked."""
    bg      = data.get("background_distribution", {})
    metrics = data.get("metrics", {})

    raw_deltas  = bg.get("background_deltas")
    deltas      = np.array(raw_deltas) if raw_deltas is not None else np.array([])
    var_delta   = abs(metrics.get("max_delta", 0))
    z_score     = metrics.get("z_score")
    percentile  = metrics.get("percentile")

    if len(deltas) < 2 or z_score is None or percentile is None:
        return _unavailable_figure(
            "Variant Impact vs. Gene Background Distribution",
            "Unavailable — no pre-computed gene background distribution",
            height=400,
        )

    kde_x = np.linspace(0, max(deltas.max(), var_delta) * 1.3, 300)
    kde   = stats.gaussian_kde(deltas, bw_method=0.4)
    kde_y = kde(kde_x)

    # Shade the tail beyond the variant
    tail_x = kde_x[kde_x >= var_delta]
    tail_y = kde_y[kde_x >= var_delta]

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=kde_x, y=kde_y, mode="lines",
        name="Background (KDE)",
        line=dict(color=C["secondary"], width=2.5),
        fill="tozeroy", fillcolor="rgba(91,133,170,0.15)",
        hovertemplate="Δ: %{x:.3f}<extra></extra>",
    ))
    if len(tail_x):
        fig.add_trace(go.Scatter(
            x=np.concatenate([[tail_x[0]], tail_x]),
            y=np.concatenate([[0], tail_y]),
            fill="tozeroy", fillcolor="rgba(196,64,39,0.25)",
            line=dict(width=0), name=f"Top {100-percentile:.1f}% tail",
            hoverinfo="skip",
        ))

    fig.add_vline(x=var_delta, line_dash="dash", line_color=C["primary"], line_width=2.5,
                  annotation_text=f"This variant<br>Z={z_score:.2f} | P{percentile:.0f}",
                  annotation_position="top right",
                  annotation=dict(font_color=C["primary"]))

    fig.update_layout(
        title="Variant Impact vs. Gene Background Distribution",
        xaxis_title="Absolute Impact Score (|Δ|)",
        yaxis_title="Density",
        height=400,
        showlegend=True,
        legend=dict(orientation="h", yanchor="bottom", y=1.01, xanchor="right", x=1),
        **_LAYOUT,
    )
    return fig


# ── Plot C: Pathogenicity Radar ────────────────────────────────────────────────
def plot_pathogenicity_radar(data: dict) -> go.Figure:
    """
    Radar chart of available normalized evidence dimensions (only includes non-missing metrics).
    """
    metrics = data.get("metrics", {})
    cons_arr = data.get("tracks", {}).get("conservation")
    gnomad_raw = metrics.get("gnomad_freq")
    z_raw = metrics.get("z_score")

    dims, vals = [], []
    if cons_arr is not None and len(cons_arr) > 0:
        dims.append("Conservation")
        vals.append(float(np.clip(np.mean(cons_arr) / 4.0, 0, 1)))
    if gnomad_raw is not None:
        dims.append("Population Rarity")
        vals.append(float(np.clip(1 - np.log10(max(gnomad_raw, 1e-8)) / -8, 0, 1)))
    dims.append("Model Impact")
    vals.append(float(np.clip(abs(metrics.get("max_delta", 0)) / 1.0, 0, 1)))
    if z_raw is not None:
        dims.append("Z-Score")
        vals.append(float(np.clip((z_raw + 3) / 6, 0, 1)))
    dims.append("Data Completeness")
    vals.append(float(np.clip(metrics.get("confidence", 0) / 100.0, 0, 1)))

    if len(dims) < 3:
        return _unavailable_figure(
            "Multi-Evidence Profile",
            "Unavailable — fewer than 3 real evidence dimensions available",
            height=420,
        )

    fig = go.Figure()
    fig.add_trace(go.Scatterpolar(
        r=vals + [vals[0]],
        theta=dims + [dims[0]],
        fill="toself",
        name="This Variant",
        fillcolor="rgba(244,96,54,0.5)",
        line=dict(color=C["primary"], width=2),
    ))

    fig.update_layout(
        polar=dict(radialaxis=dict(visible=True, range=[0, 1], tickfont_size=10)),
        title="Multi-Evidence Profile (Available Sources Only)",
        height=420,
        showlegend=False,
        **_LAYOUT,
    )
    return fig


# ── Plot D: Evidence Stack + Classification Badge ─────────────────────────────
def plot_evidence_stack(data: dict) -> tuple[go.Figure, str, str]:
    """
    Horizontal stacked bar of available evidence components; returns (fig, label, badge_color).
    """
    metrics = data.get("metrics", {})
    gnomad_raw = metrics.get("gnomad_freq")
    conserve_raw = data.get("tracks", {}).get("conservation")

    components = []
    if gnomad_raw is not None:
        rarity_score = float(np.clip(1 - np.log10(max(gnomad_raw, 1e-8)) / -8, 0, 1)) * 100
        components.append(("Population Rarity", rarity_score, C["primary"]))

    impact_score = float(np.clip(abs(metrics.get("max_delta", 0)) / 1.0, 0, 1)) * 100
    components.append(("Model Impact", impact_score, C["secondary"]))

    if conserve_raw is not None and len(conserve_raw) > 0:
        cons_score = float(np.clip(np.mean(conserve_raw) / 4.0, 0, 1)) * 100
        components.append(("Conservation", cons_score, C["success"]))

    total = sum(s for _, s, _ in components)
    raw_score = total / len(components) if components else 0.0
    if len(components) < 2:
        label = "INCOMPLETE DATA"
        badge_color = C["text_light"]
    elif raw_score >= 65:
        label = "HIGH PRIORITY"
        badge_color = C["danger"]
    elif raw_score >= 40:
        label = "MODERATE / VUS"
        badge_color = C["warning"]
    else:
        label = "LOW PRIORITY"
        badge_color = C["success"]

    fig = go.Figure()
    for name, score, col in components:
        pct = (score / total * 100.0) if total > 0 else (100.0 / len(components))
        fig.add_trace(go.Bar(orientation="h", x=[pct], y=["Evidence"],
                             name=name, marker_color=col,
                             hovertemplate=f"{name}: {pct:.1f}%<extra></extra>"))

    fig.update_layout(
        barmode="stack",
        title="Available Evidence Breakdown",
        xaxis_title="Relative Evidence Contribution (%)",
        height=200,
        xaxis=dict(range=[0, 100]),
        showlegend=True,
        legend=dict(orientation="h", y=-0.5, x=0),
        **_LAYOUT,
    )
    return fig, label, badge_color


# ── Plot E: gnomAD Allele Frequency Context ────────────────────────────────────
def plot_gnomad_context(data: dict) -> go.Figure:
    """
    Plots the variant's real gnomAD allele frequency against standard rarity thresholds
    (no simulated gnomAD background distribution).
    """
    metrics = data.get("metrics", {})
    var_freq = metrics.get("gnomad_freq")

    if var_freq is None:
        return _unavailable_figure(
            "Allele Frequency Context (gnomAD)",
            "Unavailable — gnomAD and MyVariant.info allele frequency lookups returned no data",
            height=380,
        )

    safe_freq = max(float(var_freq), 1e-8)
    var_log = float(np.log10(safe_freq))

    tiers = [
        (-8, -4, "Ultra-rare (<1e-4)", "rgba(196,64,39,0.16)"),
        (-4, -3, "Rare (1e-4–1e-3)", "rgba(244,96,54,0.14)"),
        (-3, -2, "Low-freq (1e-3–1e-2)", "rgba(233,196,106,0.18)"),
        (-2, 0,  "Common (≥1e-2)", "rgba(42,157,143,0.16)"),
    ]

    fig = go.Figure()
    for x0, x1, tier_label, fill in tiers:
        fig.add_vrect(
            x0=x0, x1=x1, fillcolor=fill, line_width=0,
            annotation_text=tier_label, annotation_position="bottom",
            annotation_font_size=10,
        )

    fig.add_trace(go.Scatter(
        x=[var_log], y=[0.5], mode="markers+text",
        marker=dict(color=C["primary"], size=16, symbol="diamond", line=dict(color="white", width=2)),
        text=[f"AF = {var_freq:.2e}"],
        textposition="top center",
        name="This variant",
        hovertemplate=f"gnomAD AF: {var_freq:.3e}<br>log₁₀(AF): {var_log:.2f}<extra></extra>",
    ))

    fig.update_layout(
        title="Allele Frequency vs. Standard Population Rarity Tiers (gnomAD)",
        xaxis=dict(title="log₁₀(Allele Frequency)", range=[-8, 0]),
        yaxis=dict(visible=False, range=[0, 1]),
        height=380,
        showlegend=False,
        **_LAYOUT,
    )
    return fig


# ── Plot F: ClinVar Lollipop Plot ─────────────────────────────────────────────
def plot_clinvar_lollipop(data: dict, chrom: str, pos: int) -> go.Figure:
    """
    Lollipop plot of real nearby ClinVar variants (±50 kb, GRCh38) from NCBI ClinVar API.
    Displays an explicit unavailable/empty state when no real records are returned.
    """
    window = 50000
    clinvar_region = data.get("clinvar_region")

    if clinvar_region is None:
        return _unavailable_figure(
            f"ClinVar Variants in Region ({chrom}:{pos-window:,}–{pos+window:,})",
            "Unavailable — NCBI ClinVar region lookup failed",
            height=400,
        )

    if len(clinvar_region) == 0:
        return _unavailable_figure(
            f"ClinVar Variants in Region ({chrom}:{pos-window:,}–{pos+window:,})",
            "No ClinVar variants reported in this ±50 kb GRCh38 window",
            height=400,
        )

    SIG_MAP = {
        "pathogenic": (3, C["danger"], "Pathogenic"),
        "likely pathogenic": (2, "#E07050", "Likely pathogenic"),
        "uncertain": (1, C["warning"], "VUS"),
        "conflicting": (1, C["warning"], "Conflicting / VUS"),
        "likely benign": (-1, "#80B0D0", "Likely benign"),
        "benign": (-2, C["success"], "Benign"),
    }

    def _classify(sig_text: str):
        s = (sig_text or "").lower()
        for key in ("likely pathogenic", "pathogenic", "likely benign", "benign", "conflicting", "uncertain"):
            if key in s:
                return SIG_MAP[key]
        return (0, C["text_light"], sig_text or "Other")

    fig = go.Figure()
    for entry in clinvar_region:
        vpos = entry.get("pos")
        if vpos is None:
            continue
        sig_raw = entry.get("clinical_significance", "Uncertain significance")
        title = entry.get("title", "")
        y_val, col, _ = _classify(sig_raw)
        is_query = (vpos == pos)
        fig.add_shape(type="line", x0=vpos, x1=vpos, y0=0, y1=y_val, line=dict(color=col, width=1.5))
        fig.add_trace(go.Scatter(
            x=[vpos], y=[y_val], mode="markers",
            marker=dict(color=col, size=15 if is_query else 10,
                        symbol="star" if is_query else "circle",
                        line=dict(color="white", width=1.5)),
            name=sig_raw,
            hovertemplate=f"{title}<br>Pos: {vpos:,}<br>ClinVar: {sig_raw}<extra></extra>",
            showlegend=False,
        ))

    fig.add_vline(x=pos, line_dash="dot", line_color=C["primary"], opacity=0.8,
                  annotation_text="Query site", annotation_position="top right")
    fig.add_hline(y=0, line_color=C["text_light"], opacity=0.5)

    for y_val, col, label in [
        (3, C["danger"], "Pathogenic"),
        (2, "#E07050", "Likely pathogenic"),
        (1, C["warning"], "VUS / Conflicting"),
        (-1, "#80B0D0", "Likely benign"),
        (-2, C["success"], "Benign"),
    ]:
        fig.add_annotation(x=pos - window * 0.95, y=y_val,
                           text=label, showarrow=False,
                           font=dict(color=col, size=10), xanchor="left")

    fig.update_layout(
        title=f"ClinVar Variants in Region ({chrom}:{pos-window:,}–{pos+window:,}) — NCBI Live Data",
        xaxis_title="Genomic Position (GRCh38)",
        yaxis=dict(visible=False, range=[-3.5, 4]),
        showlegend=False,
        height=400,
        **_LAYOUT,
    )
    return fig


# ── Plot G: Tissue Heatmap ────────────────────────────────────────────────────
def plot_tissue_heatmap(tissue_effects: list) -> go.Figure:
    """
    Horizontal gradient bar chart of tissue-specific impact, sorted by magnitude.
    """
    if not tissue_effects:
        fig = go.Figure()
        fig.add_annotation(text="No tissue data available", x=0.5, y=0.5,
                           showarrow=False, xref="paper", yref="paper")
        return fig

    sorted_te = sorted(tissue_effects, key=lambda x: x["delta"], reverse=True)
    tissues   = [t["tissue"] for t in sorted_te]
    deltas    = [t["delta"]  for t in sorted_te]
    max_d     = max(deltas) if deltas else 1

    # Color gradient: orange (high) → blue (low)
    colors = [
        f"rgba({int(244*(d/max_d))},{int(96 + 50*(1-d/max_d))},{int(54*(1-d/max_d))},0.85)"
        for d in deltas
    ]

    fig = go.Figure(go.Bar(
        x=deltas,
        y=tissues,
        orientation="h",
        marker_color=colors,
        hovertemplate="%{y}: Δ=%{x:.4f}<extra></extra>",
        name="Tissue Impact",
    ))

    fig.update_layout(
        title="Predicted Tissue-Specific Impact (Enformer)",
        xaxis_title="Impact Score (|Δ|)",
        yaxis=dict(autorange="reversed"),
        height=max(300, 50 + 35 * len(tissues)),
        showlegend=False,
        **_LAYOUT,
    )
    return fig


# ── Plot H: Enformer Top Tracks Heatmap ──────────────────────────────────────
def plot_enformer_tracks(data: dict) -> go.Figure:
    """
    Heatmap of the real top-20 Enformer output tracks with highest |delta| at the variant bin.
    Displays an explicit unavailable state when real per-track outputs are absent.
    """
    top_tracks = data.get("enformer_top_tracks")
    if not top_tracks:
        return _unavailable_figure(
            "Top 20 Enformer Tracks — Differential Signal",
            "Unavailable — real per-track Enformer deltas were not returned",
            height=400,
        )

    y_labels = [t["track"] for t in top_tracks[:20]]
    signed_deltas = [float(t["delta"]) for t in top_tracks[:20]]
    z_vals = [[v] for v in signed_deltas]

    fig = go.Figure(go.Heatmap(
        z=z_vals,
        x=["Δ Signal"],
        y=y_labels,
        colorscale=[
            [0.0, "#2166AC"], [0.4, "#D1E5F0"],
            [0.5, "#f7f7f7"],
            [0.6, "#FDDBC7"], [1.0, C["danger"]],
        ],
        zmid=0,
        colorbar=dict(title="Δ Signal", thickness=12, len=0.8),
        hovertemplate="Track: %{y}<br>Δ: %{z:.4f}<extra></extra>",
        text=[[f"{v:+.4f}"] for v in signed_deltas],
        texttemplate="%{text}",
        textfont=dict(size=10),
    ))

    fig.update_layout(
        title="Top 20 Enformer Human Tracks — Real Differential Signal",
        xaxis=dict(visible=False),
        yaxis=dict(tickfont=dict(size=10), autorange="reversed"),
        height=600,
        **_LAYOUT,
    )
    return fig
