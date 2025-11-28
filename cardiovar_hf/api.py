import warnings
# Suppress numpy binary incompatibility warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

from fastapi import FastAPI, HTTPException, Body
from pydantic import BaseModel
from typing import List, Optional
import json
import os
from variant_engine import compute_variant_impact

app = FastAPI(title="CardioVar API", version="1.0")

# --- Models ---
class VariantRequest(BaseModel):
    assembly: str = "GRCh38"  # Genome build
    chrom: str
    pos: int
    ref: str
    alt: str
    window_size: Optional[int] = 100
    
    class Config:
        # Validate assembly values
        schema_extra = {
            "example": {
                "assembly": "GRCh38",
                "chrom": "chr22",
                "pos": 36191400,
                "ref": "A",
                "alt": "C",
                "window_size": 100
            }
        }

class BatchRequest(BaseModel):
    variants: List[VariantRequest]

# --- Endpoints ---

@app.post("/variant-impact")
def get_variant_impact(req: VariantRequest):
    """
    Compute impact for a single variant.
    """
    try:
        result = compute_variant_impact(req.chrom, req.pos, req.ref, req.alt, req.assembly, req.window_size)
        return result
    except ValueError as e:
        # Return 400 for validation errors (e.g., unsupported assembly)
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/gene-annotations")
def get_gene_annotations(gene: str):
    """
    Get gene annotations from Ensembl API (with fallback to local data).
    """
    from api_integrations import fetch_ensembl_gene, fetch_gtex_expression, load_fallback_gene_data
    
    try:
        # Try real APIs first
        ensembl_data = fetch_ensembl_gene(gene)
        gtex_data = fetch_gtex_expression(gene)
        
        if ensembl_data:
            # Merge with local data for additional fields
            local_data = load_fallback_gene_data(gene)
            
            # Prepare base response
            response = {
                "symbol": gene,
                "name": ensembl_data.get("name"),
                "ensembl_id": ensembl_data.get("ensembl_id"),
                "chromosome": ensembl_data.get("chromosome"),
                "start": ensembl_data.get("start"),
                "end": ensembl_data.get("end"),
                "biotype": ensembl_data.get("biotype"),
                "source": "Ensembl API"
            }
            
            # Add GTEx data if available
            if gtex_data:
                response["expression"] = gtex_data
                response["source"] += " + GTEx API"
            
            # Add local data (domains, pathways, etc.)
            if local_data:
                response.update({k: v for k, v in local_data.items() if k not in response})
                if "expression" not in response: # Use local expression if GTEx failed
                    response["expression"] = local_data.get("expression")
                response["source"] += " + Local"
                
            return response
        
        # Fallback to local data if Ensembl fails
        local_data = load_fallback_gene_data(gene)
        if local_data:
            return {**local_data, "source": "Local (API unavailable)"}
        
        raise HTTPException(status_code=404, detail=f"Gene {gene} not found")
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/related-data")
def get_related_data(req: VariantRequest):
    """
    Get related variant data (ClinVar, GWAS) with fallback to local data.
    """
    from api_integrations import fetch_clinvar_data, load_fallback_related_data
    
    try:
        # Try real API (placeholder for now)
        clinvar_data = fetch_clinvar_data(req.chrom, req.pos, req.ref, req.alt)
        
        # Always use local data for now (ClinVar API is complex)
        local_data = load_fallback_related_data(req.chrom, req.pos, req.ref, req.alt)
        
        # Combine results
        all_data = clinvar_data + local_data
        
        return all_data if all_data else []
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

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
