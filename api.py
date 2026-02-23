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
    force_live: Optional[bool] = False
    
    class Config:
        # Validate assembly values
        schema_extra = {
            "example": {
                "assembly": "GRCh38",
                "chrom": "chr22",
                "pos": 36191400,
                "ref": "A",
                "alt": "C",
                "window_size": 100,
                "force_live": False
            }
        }

class BatchRequest(BaseModel):
    variants: List[VariantRequest]

# --- Endpoints ---

import uuid
from fastapi import BackgroundTasks
import psutil

# --- Global State ---
# --- Global State ---
BATCH_JOBS = {}
SINGLE_JOBS = {}  # Store full results for single variant analysis

# --- Models ---
class BatchResponse(BaseModel):
    batch_id: str
    message: str

class JobResponse(BaseModel):
    job_id: str
    message: str

# --- Helper Functions ---
def process_single_variant_task(job_id: str, req: VariantRequest):
    """
    Background task to process a single variant with FULL details.
    """
    try:
        SINGLE_JOBS[job_id]["status"] = "processing"
        
        # Compute full impact with all tracks and curves
        result = compute_variant_impact(
            req.chrom, 
            req.pos, 
            req.ref, 
            req.alt, 
            req.assembly, 
            req.window_size,
            force_live=req.force_live
        )
        
        SINGLE_JOBS[job_id]["result"] = result
        SINGLE_JOBS[job_id]["status"] = "completed"
        
    except Exception as e:
        SINGLE_JOBS[job_id]["status"] = "failed"
        SINGLE_JOBS[job_id]["error"] = str(e)

def process_batch_task(batch_id: str, variants: List[VariantRequest]):
    """
    Background task to process variants.
    """
    try:
        BATCH_JOBS[batch_id]["status"] = "processing"
        results = []
        
        for i, v in enumerate(variants):
            try:
                # Process variant
                res = compute_variant_impact(v.chrom, v.pos, v.ref, v.alt, force_live=False)
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
                    "priority": priority,
                    "status": "success"
                })
            except Exception as e:
                # Handle individual variant failure
                results.append({
                    "variant_id": f"{v.chrom}:{v.pos}:{v.ref}:{v.alt}",
                    "status": "failed",
                    "error": str(e)
                })
            
            # Update progress
            BATCH_JOBS[batch_id]["processed"] = i + 1
            
        BATCH_JOBS[batch_id]["results"] = results
        BATCH_JOBS[batch_id]["status"] = "completed"
        
    except Exception as e:
        BATCH_JOBS[batch_id]["status"] = "failed"
        BATCH_JOBS[batch_id]["error"] = str(e)

# --- Endpoints ---

@app.post("/variant-impact")
def get_variant_impact(req: VariantRequest):
    """
    Compute impact for a single variant (Synchronous - Legacy).
    """
    try:
        result = compute_variant_impact(
            req.chrom, 
            req.pos, 
            req.ref, 
            req.alt, 
            req.assembly, 
            req.window_size,
            force_live=req.force_live
        )
        return result
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/single-cell-expression")
def get_single_cell_expression(payload: dict = Body(...)):
    """
    Fetch single-cell expression data for a given gene.
    Payload: {"gene_symbol": "MYH9"}
    """
    from api_integrations import fetch_single_cell_expression
    gene_symbol = payload.get("gene_symbol")
    if not gene_symbol:
        raise HTTPException(status_code=400, detail="gene_symbol is required")
        
    try:
        data = fetch_single_cell_expression(gene_symbol)
        if data is None:
            # Return empty list instead of 404 to handle gracefully on frontend
            return []
        return data
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/variant-impact-async", response_model=JobResponse)
def get_variant_impact_async(req: VariantRequest, background_tasks: BackgroundTasks):
    """
    Start async impact computation for a single variant.
    """
    job_id = str(uuid.uuid4())
    SINGLE_JOBS[job_id] = {
        "status": "pending",
        "result": None
    }
    background_tasks.add_task(process_single_variant_task, job_id, req)
    return {"job_id": job_id, "message": "Analysis started"}

@app.get("/job-status/{job_id}")
def get_job_status(job_id: str):
    """
    Get status and result of a single variant job.
    """
    if job_id not in SINGLE_JOBS:
        raise HTTPException(status_code=404, detail="Job ID not found")
    return SINGLE_JOBS[job_id]

@app.get("/gene-annotations")
def get_gene_annotations(gene: str, force_live: bool = False):
    """
    Get gene annotations from gnomAD/Ensembl API (with fallback to local data).
    """
    from api_integrations import (
        fetch_gnomad_gene,  # Try gnomAD first - faster!
        fetch_ensembl_gene, 
        fetch_gtex_expression, 
        fetch_protein_domains,
        load_fallback_gene_data
    )
    
    try:
        # Try gnomAD first (faster than Ensembl)
        ensembl_data = fetch_gnomad_gene(gene, force_live=force_live)
        print(f"DEBUG: ensembl_data type: {type(ensembl_data)}")
        
        gtex_data = fetch_gtex_expression(gene, force_live=force_live)
        print(f"DEBUG: gtex_data type: {type(gtex_data)}")
        
        protein_data = fetch_protein_domains(gene, force_live=force_live)
        print(f"DEBUG: protein_data type: {type(protein_data)}")
        if protein_data:
             print(f"DEBUG: protein_data content: {protein_data}")
             if isinstance(protein_data, list):
                 print("ERROR: protein_data is a list! Expected dict.")
                 protein_data = None # Prevent crash
        
        if ensembl_data:
            # Merge with local data for additional fields
            local_data = load_fallback_gene_data(gene)
            
            # Prepare base response
            response = {
                "symbol": gene,
                "name": ensembl_data.get("description"),
                "ensembl_id": ensembl_data.get("id"),
                "chromosome": ensembl_data.get("seq_region_name"),
                "start": ensembl_data.get("start"),
                "end": ensembl_data.get("end"),
                "biotype": ensembl_data.get("biotype"),
                "source": "Ensembl API"
            }
            
            # Add links
            if "links" in ensembl_data:
                response["links"] = ensembl_data["links"]
            
            # Add GTEx data
            if gtex_data:
                response["expression"] = gtex_data
                response["source"] += " + GTEx API"
            
            # Add protein domain data
            if protein_data:
                response["protein_length"] = protein_data["protein_length"]
                response["protein_domains"] = protein_data["protein_domains"]
                response["source"] += " + Ensembl Protein Features"
            
            # Add local data
            if local_data:
                for key in ["pathways", "disease_associations"]:
                    if key in local_data:
                        response[key] = local_data[key]
                
                if "expression" not in response and "expression" in local_data:
                    response["expression"] = local_data["expression"]
                
                if "protein_domains" not in response:
                    if "protein_domains" in local_data:
                        response["protein_domains"] = local_data["protein_domains"]
                    if "protein_length" in local_data:
                        response["protein_length"] = local_data["protein_length"]
                
                response["source"] += " + Local"
                
            return response
        
        # Fallback to local data
        local_data = load_fallback_gene_data(gene)
        if local_data:
            return {**local_data, "source": "Local (API unavailable)"}
        
        raise HTTPException(status_code=404, detail=f"Gene {gene} not found")
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/related-data")
def get_related_data(chrom: str, pos: int, ref: str, alt: str, force_live: bool = False):
    """
    Get related variant data (ClinVar, dbSNP) with fallback to local data.
    """
    from api_integrations import (
        fetch_clinvar_variants, 
        fetch_dbsnp_variants,
        load_fallback_related_data
    )
    
    try:
        # Fetch data from new wrappers
        clinvar = fetch_clinvar_variants(chrom, pos, ref, alt, force_live=force_live)
        dbsnp = fetch_dbsnp_variants(chrom, pos, ref, alt, force_live=force_live)
        
        # Load local fallback data (if any)
        try:
            local_data = load_fallback_related_data(chrom, pos, ref, alt)
        except ImportError:
            local_data = []
            
        # Return structured response
        return {
            "clinvar": clinvar,
            "dbsnp": dbsnp,
            "local_fallback": local_data
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

# --- Admin / Cache Management ---

class CacheInvalidateRequest(BaseModel):
    key: Optional[str] = None
    pattern: Optional[str] = None

@app.post("/admin/cache/invalidate")
def invalidate_cache(req: CacheInvalidateRequest):
    """
    Invalidate cache entries by key or pattern.
    """
    from api_integrations import cache
    
    if req.key:
        cache.invalidate(req.key)
        return {"message": f"Invalidated key: {req.key}"}
    elif req.pattern:
        cache.invalidate_pattern(req.pattern)
        return {"message": f"Invalidated pattern: {req.pattern}"}
    else:
        raise HTTPException(status_code=400, detail="Must provide 'key' or 'pattern'")

@app.post("/batch-start", response_model=BatchResponse)
def start_batch(req: BatchRequest, background_tasks: BackgroundTasks):
    """
    Start a batch processing job in the background.
    """
    batch_id = str(uuid.uuid4())
    BATCH_JOBS[batch_id] = {
        "status": "pending",
        "total": len(req.variants),
        "processed": 0,
        "results": []
    }
    
    background_tasks.add_task(process_batch_task, batch_id, req.variants)
    
    return {"batch_id": batch_id, "message": "Batch processing started"}

@app.get("/batch-status/{batch_id}")
def get_batch_status(batch_id: str):
    """
    Get the status and results of a batch job.
    """
    if batch_id not in BATCH_JOBS:
        raise HTTPException(status_code=404, detail="Batch ID not found")
    
    return BATCH_JOBS[batch_id]

@app.get("/system-status")
def get_system_status():
    """
    Get current system resource usage.
    """
    try:
        cpu_percent = psutil.cpu_percent(interval=None)
        memory = psutil.virtual_memory()
        return {
            "cpu_percent": cpu_percent,
            "memory_percent": memory.percent,
            "memory_used_gb": round(memory.used / (1024**3), 2),
            "memory_total_gb": round(memory.total / (1024**3), 2)
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
