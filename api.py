from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
from typing import List, Optional
import json
import os
import pandas as pd
from variant_engine import compute_variant_impact

app = FastAPI(title="CardioVar API", version="1.0")

# --- Models ---
class VariantRequest(BaseModel):
    chrom: str
    pos: int
    ref: str
    alt: str
    window_size: Optional[int] = 100

class BatchRequest(BaseModel):
    variants: List[VariantRequest]

# --- Data Loading ---
def load_json(filename):
    path = os.path.join("data", filename)
    if os.path.exists(path):
        with open(path, "r") as f:
            return json.load(f)
    return {}

GENE_DATA = load_json("gene_annotations.json")
RELATED_DATA = load_json("related_variants.json")

# --- Endpoints ---

@app.post("/variant-impact")
def get_variant_impact(req: VariantRequest):
    """
    Compute impact for a single variant.
    """
    try:
        result = compute_variant_impact(req.chrom, req.pos, req.ref, req.alt, req.window_size)
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/gene-annotations")
def get_gene_annotations(gene: str):
    """
    Get annotations for a specific gene.
    """
    # Simple linear search for mock data
    for record in GENE_DATA:
        if record["symbol"].upper() == gene.upper():
            return record
    return {"symbol": gene, "note": "No detailed annotations found."}

@app.post("/related-data")
def get_related_data(req: VariantRequest):
    """
    Get ClinVar/GWAS data for a variant.
    """
    key = f"{req.chrom}:{req.pos}:{req.ref}:{req.alt}"
    return RELATED_DATA.get(key, [])

@app.post("/batch-impact")
def batch_impact(req: BatchRequest):
    """
    Process a batch of variants.
    """
    results = []
    for v in req.variants:
        res = compute_variant_impact(v.chrom, v.pos, v.ref, v.alt)
        metrics = res["metrics"]
        
        # Determine Priority
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
    return results

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
