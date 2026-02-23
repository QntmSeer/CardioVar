"""
API integration utilities for CardioVar.
All external calls include a User-Agent header, retry logic, smart caching with TTL, and graceful fallbacks to locally cached reference data when the public APIs are unreachable.
"""

import os
import json
import time
import logging
import requests
from typing import Optional, Dict, Any, List
import pandas as pd
import numpy as np
from api_cache import APICache

# Global HTTP session for connection reuse
SESSION = requests.Session()
SESSION.headers.update({"User-Agent": "CardioVar/1.0 (+https://github.com/yourorg/CardioVar)"})

# ---------------------------------------------------------------------------
# Configuration & Globals
# ---------------------------------------------------------------------------
logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')

# Initialize cache (24‑hour TTL by default)
cache = APICache(db_path="data/api_cache.db", default_ttl_hours=24)

# Request headers
HEADERS = {"User-Agent": "CardioVar/1.0 (+https://github.com/yourorg/CardioVar)"}

# API endpoints
GNOMAD_API = "https://gnomad.broadinstitute.org/api"
ENSEMBL_API = "https://rest.ensembl.org"
UCSC_API = "https://api.genome.ucsc.edu"
GTEX_API = "https://gtexportal.org/rest/v1"

# Fallback data directory (relative to this file)
FALLBACK_DIR = os.path.join(os.path.dirname(__file__), "data")

# Global flag indicating whether a fallback was used in the most recent request
fallback_used = False

def reset_fallback_flag():
    """Reset the global ``fallback_used`` flag.
    Call this before a series of API calls if you want to know whether any of them fell back.
    """
    global fallback_used
    fallback_used = False

# ---------------------------------------------------------------------------
# Helper: retry decorator
# ---------------------------------------------------------------------------
def retry_on_failure(retries: int = 3, backoff: float = 2.0):
    """Retry a function on exception or ``None`` result.
    Args:
        retries: Number of attempts.
        backoff: Multiplier for sleep between attempts.
    """
    def decorator(func):
        def wrapper(*args, **kwargs):
            delay = backoff
            for attempt in range(retries):
                try:
                    result = func(*args, **kwargs)
                    if result is not None:
                        return result
                except Exception as e:
                    logging.debug(f"{func.__name__} attempt {attempt+1} failed: {e}")
                if attempt < retries - 1:
                    time.sleep(delay)
                    delay *= backoff
            return None
        return wrapper
    return decorator

# ---------------------------------------------------------------------------
# Fallback helpers
# ---------------------------------------------------------------------------
def _load_json_fallback(filename: str) -> List[Dict[str, Any]]:
    path = os.path.join(FALLBACK_DIR, filename)
    if not os.path.exists(path):
        logging.debug(f"Fallback file not found: {path}")
        return []
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)

def load_fallback_gene_data(gene_symbol: str) -> Optional[Dict[str, Any]]:
    # Fix: Use correct filename 'gene_annotations.json' instead of 'gene_annotations_fallback.json'
    data = _load_json_fallback("gene_annotations.json")
    for entry in data:
        if entry.get("symbol") == gene_symbol:
            return entry
    return None

# ... (lines 94-368 omitted) ...

@retry_on_failure(retries=3, backoff=2.0)
def fetch_protein_domains(gene_symbol: str, force_live: bool = False) -> Optional[Dict[str, Any]]:
    """Fetch protein domain information from Ensembl with caching and JSON fallback."""
    global fallback_used
    cache_key = f"protein_domains:{gene_symbol}"
    cached = cache.get(cache_key)
    if cached is not None:
        logging.debug(f"Cache hit for {cache_key}")
        return cached

    try:
        gene_url = f"{ENSEMBL_API}/lookup/symbol/homo_sapiens/{gene_symbol}?expand=1"
        gene_resp = SESSION.get(gene_url, timeout=15)
        if gene_resp.status_code != 200:
            raise RuntimeError("Gene lookup failed")
        gene_data = gene_resp.json()
        canonical = gene_data.get("canonical_transcript")
        if not canonical:
            raise RuntimeError("No canonical transcript")
        protein_url = f"{ENSEMBL_API}/overlap/translation/{canonical}?feature=protein_feature;content-type=application/json"
        protein_resp = SESSION.get(protein_url, timeout=15)
        if protein_resp.status_code != 200:
            raise RuntimeError("Protein features fetch failed")
        features = protein_resp.json()
        domains = []
        protein_length = 0
        for feat in features:
            if feat.get("type") in ["domain", "pfam"]:
                domains.append({
                    "name": feat.get("id", "Unknown"),
                    "start": feat.get("start", 0),
                    "end": feat.get("end", 0),
                    "type": feat.get("type", "domain"),
                })
                protein_length = max(protein_length, feat.get("end", 0))
        if domains:
            result = {"protein_length": protein_length, "protein_domains": domains}
            cache.set(cache_key, result)
            return result
    except Exception as e:
        logging.debug(f"Protein domains API error: {e}")

    if force_live:
        logging.warning("Force live Protein Domains failed; not attempting fallback.")
        fallback_used = True
        return None

    # Fallback JSON
    try:
        path = os.path.join(FALLBACK_DIR, "protein_domains_fallback.json")
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
        
        # Fix: Filter by gene_symbol
        result = None
        for item in data:
            if item.get("gene_symbol") == gene_symbol:
                result = {
                    "protein_length": item.get("protein_length"),
                    "protein_domains": item.get("domains")
                }
                break
        
        if result:
            fallback_used = True
            cache.set(cache_key, result)
            return result
            
    except Exception as e:
        logging.debug(f"Protein domains fallback load error: {e}")
        
    fallback_used = True
    return None

# ---------------------------------------------------------------------------
# API Functions (each uses cache + fallback)
# ---------------------------------------------------------------------------
@retry_on_failure(retries=3, backoff=2.0)
def fetch_gnomad_frequency(chrom: str, pos: int, ref: str, alt: str, force_live: bool = False) -> Optional[float]:
    """Fetch allele frequency from gnomAD v2 REST API with caching and fallback.
    Returns ``None`` if no data is available.
    """
    global fallback_used
    cache_key = f"gnomad:{chrom}:{pos}:{ref}:{alt}"
    cached = cache.get(cache_key)
    if cached is not None:
        logging.debug(f"Cache hit for {cache_key}")
        return cached

    try:
        url = f"{GNOMAD_API}/variant/{chrom}-{pos}-{ref}-{alt}"
        resp = SESSION.get(url, timeout=10)
        if resp.status_code == 200:
            data = resp.json()
            freq = data.get("variant", {}).get("exome", {}).get("ac", 0) / max(
                data.get("variant", {}).get("exome", {}).get("an", 1), 1
            )
            cache.set(cache_key, freq)
            return freq
    except Exception as e:
        logging.debug(f"gnomAD API error: {e}")

    if force_live:
        logging.warning("Force live gnomAD failed; not attempting fallback.")
        fallback_used = True
        return None

    # Fallback JSON
    fallback = _load_json_fallback("gnomad_fallback.json")
    
    # Check if fallback is a list (old format) or dict (new format)
    if isinstance(fallback, list):
        for entry in fallback:
            if (
                entry.get("chrom") == chrom
                and entry.get("pos") == pos
                and entry.get("ref") == ref
                and entry.get("alt") == alt
            ):
                fallback_used = True
                cache.set(cache_key, entry.get("af"))
                return entry.get("af")
    elif isinstance(fallback, dict):
        # Construct key: chrom-pos-ref-alt (chrom without 'chr')
        clean_chrom = chrom.replace("chr", "")
        key = f"{clean_chrom}-{pos}-{ref}-{alt}"
        if key in fallback:
            fallback_used = True
            val = fallback[key]
            cache.set(cache_key, val)
            return val
            
    fallback_used = True
    return None

@retry_on_failure(retries=3, backoff=2.0)
def fetch_ensembl_gene(gene_symbol: str, force_live: bool = False) -> Optional[Dict[str, Any]]:
    """Fetch gene information from Ensembl REST API with caching and fallback."""
    global fallback_used
    cache_key = f"ensembl_gene:{gene_symbol}"
    cached = cache.get(cache_key)
    if cached is not None:
        logging.debug(f"Cache hit for {cache_key}")
        return cached

    url = f"{ENSEMBL_API}/lookup/symbol/homo_sapiens/{gene_symbol}?expand=1"
    try:
        resp = SESSION.get(url, timeout=15)
        if resp.status_code == 200:
            data = resp.json()
            # Add external links
            links = {
                "GeneCards": f"https://www.genecards.org/cgi-bin/carddisp.pl?gene={gene_symbol}",
                "UniProt": f"https://www.uniprot.org/uniprot/?query={gene_symbol}&sort=score",
                "gnomAD": f"https://gnomad.broadinstitute.org/gene/{data.get('id')}?dataset=gnomad_r4",
                "Ensembl": f"https://www.ensembl.org/Homo_sapiens/Gene/Summary?g={data.get('id')}",
                "OMIM": f"https://www.omim.org/search?index=entry&start=1&limit=10&sort=score+desc%2C+prefix_sort+desc&search={gene_symbol}",
            }
            data["links"] = links
            data["source"] = "Ensembl"
            cache.set(cache_key, data)
            return data
    except Exception as e:
        logging.debug(f"Ensembl API error: {e}")

    if force_live:
        logging.warning("Force live Ensembl Gene failed; not attempting fallback.")
        fallback_used = True
        return None

    # Fallback JSON
    fallback_data = load_fallback_gene_data(gene_symbol)
    if fallback_data:
        # Add links to fallback data too
        links = {
            "GeneCards": f"https://www.genecards.org/cgi-bin/carddisp.pl?gene={gene_symbol}",
            "UniProt": f"https://www.uniprot.org/uniprot/?query={gene_symbol}&sort=score",
            "gnomAD": f"https://gnomad.broadinstitute.org/gene/{fallback_data.get('ensembl_id')}?dataset=gnomad_r4",
            "Ensembl": f"https://www.ensembl.org/Homo_sapiens/Gene/Summary?g={fallback_data.get('ensembl_id')}",
            "OMIM": f"https://www.omim.org/search?index=entry&start=1&limit=10&sort=score+desc%2C+prefix_sort+desc&search={gene_symbol}",
        }
        fallback_data["links"] = links
        fallback_data["source"] = "Local Fallback"
        fallback_used = True
        cache.set(cache_key, fallback_data)
        return fallback_data

    fallback_used = True
    return None


@retry_on_failure(retries=2, backoff=1.5)
def fetch_gnomad_gene(gene_symbol: str, force_live: bool = False) -> Optional[Dict[str, Any]]:
    """Fetch gene information from gnomAD GraphQL API - faster alternative to Ensembl."""
    global fallback_used
    cache_key = f"gnomad_gene:{gene_symbol}"
    cached = cache.get(cache_key)
    if cached is not None:
        logging.debug(f"Cache hit for {cache_key}")
        return cached

    graphql_url = "https://gnomad.broadinstitute.org/api"
    
    # GraphQL query for gene information
    query = """
    query GeneInfo($geneSymbol: String!, $referenceGenome: ReferenceGenomeId!) {
      gene(gene_symbol: $geneSymbol, reference_genome: $referenceGenome) {
        gene_id
        gene_version
        symbol
        name
        chrom
        start
        stop
        strand
        canonical_transcript_id
        mane_select_transcript {
          ensembl_id
          ensembl_version
          refseq_id
          refseq_version
        }
        pext {
          tissues {
            tissue
            mean
          }
        }
      }
    }
    """
    
    variables = {
        "geneSymbol": gene_symbol,
        "referenceGenome": "GRCh38"
    }
    
    try:
        resp = SESSION.post(
            graphql_url,
            json={"query": query, "variables": variables},
            timeout=10,  # Shorter timeout than Ensembl
            headers={"Content-Type": "application/json"}
        )
        
        if resp.status_code == 200:
            data = resp.json()
            
            if "errors" in data:
                logging.debug(f"gnomAD GraphQL errors: {data['errors']}")
                raise RuntimeError("GraphQL query failed")
            
            gene_data = data.get("data", {}).get("gene")
            if gene_data:
                # Format to match Ensembl-style response
                result = {
                    "id": gene_data.get("gene_id"),
                    "symbol": gene_data.get("symbol"),
                    "description": gene_data.get("name", f"Gene {gene_symbol}"),
                    "seq_region_name": gene_data.get("chrom"),
                    "start": gene_data.get("start"),
                    "end": gene_data.get("stop"),
                    "strand": 1 if gene_data.get("strand") == "+" else -1,
                    "canonical_transcript": gene_data.get("canonical_transcript_id"),
                    "mane_select": gene_data.get("mane_select_transcript"),
                    "source": "gnomAD",
                    "links": {
                        "GeneCards": f"https://www.genecards.org/cgi-bin/carddisp.pl?gene={gene_symbol}",
                        "UniProt": f"https://www.uniprot.org/uniprot/?query={gene_symbol}&sort=score",
                        "gnomAD": f"https://gnomad.broadinstitute.org/gene/{gene_data.get('gene_id')}?dataset=gnomad_r4",
                        "Ensembl": f"https://www.ensembl.org/Homo_sapiens/Gene/Summary?g={gene_data.get('gene_id')}",
                        "OMIM": f"https://www.omim.org/search?index=entry&start=1&limit=10&sort=score+desc%2C+prefix_sort+desc&search={gene_symbol}",
                    }
                }
                
                cache.set(cache_key, result)
                logging.info(f"gnomAD gene data retrieved for {gene_symbol}")
                return result
            else:
                logging.debug(f"No gnomAD data for gene {gene_symbol}")
                
    except Exception as e:
        logging.debug(f"gnomAD gene API error: {e}")

    if force_live:
        fallback_used = True
        return None

    # Try Ensembl as fallback
    return fetch_ensembl_gene(gene_symbol, force_live=False)

@retry_on_failure(retries=3, backoff=2.0)
def fetch_ucsc_phylop(chrom: str, start: int, end: int, force_live: bool = False) -> Optional[List[float]]:
    """Fetch PhyloP scores from UCSC with caching and fallback (numpy file)."""
    global fallback_used
    cache_key = f"phylop:{chrom}:{start}:{end}"
    cached = cache.get(cache_key)
    if cached is not None:
        logging.debug(f"Cache hit for {cache_key}")
        return cached

    if not chrom.startswith("chr"):
        chrom = f"chr{chrom}"
    url = (
        f"{UCSC_API}/getData/track?genome=hg38;track=phyloP100way;"
        f"chrom={chrom};start={start};end={end}"
    )
    try:
        resp = SESSION.get(url, timeout=20)
        if resp.status_code == 200:
            data = resp.json()
            if "phyloP100way" in data:
                scores = [float(item["value"]) for item in data["phyloP100way"] if "value" in item]
                expected_len = end - start
                if scores and len(scores) != expected_len:
                    # Interpolate missing positions
                    x_old = np.linspace(0, 1, len(scores))
                    x_new = np.linspace(0, 1, expected_len)
                    scores = np.interp(x_new, x_old, scores).tolist()
                if len(scores) == expected_len:
                    cache.set(cache_key, scores)
                    return scores
    except Exception as e:
        logging.debug(f"UCSC API error: {e}")

    if force_live:
        logging.warning("Force live PhyloP failed; not attempting fallback.")
        fallback_used = True
        return None

    # Fallback: numpy cache file
    try:
        path = os.path.join(FALLBACK_DIR, "phylop_fallback.npy")
        np_cache = np.load(path, allow_pickle=True).item()
        key = (chrom, start, end)
        val = np_cache.get(key)
        if val is not None:
            fallback_used = True
            cache.set(cache_key, val)
            return val
    except Exception as e:
        logging.debug(f"PhyloP fallback load error: {e}")
    fallback_used = True
    return None

@retry_on_failure(retries=3, backoff=2.0)
def fetch_gtex_expression(gene_symbol: str, force_live: bool = False) -> Optional[List[Dict[str, Any]]]:
    """Fetch GTEx expression data with caching and TSV fallback."""
    global fallback_used
    cache_key = f"gtex:{gene_symbol}"
    cached = cache.get(cache_key)
    if cached is not None:
        logging.debug(f"Cache hit for {cache_key}")
        return cached

    url = f"{GTEX_API}/expression/medianGeneExpression"
    params = {"geneId": gene_symbol, "datasetId": "gtex_v8"}
    try:
        resp = SESSION.get(url, params=params, timeout=15)
        if resp.status_code == 200:
            data = resp.json()
            if "medianGeneExpression" in data:
                result = [
                    {"tissue": item["tissueSiteDetailId"], "tpm": item["median"]}
                    for item in data["medianGeneExpression"]
                ]
                cache.set(cache_key, result)
                return result
    except Exception as e:
        logging.debug(f"GTEx API error: {e}")

    if force_live:
        logging.warning("Force live GTEx failed; not attempting fallback.")
        fallback_used = True
        return None

    # Fallback TSV
    try:
        path = os.path.join(FALLBACK_DIR, "gtex_expression.tsv")
        df = pd.read_csv(path, sep="\t")
        rows = df[df["gene_symbol"].str.upper() == gene_symbol.upper()]
        if not rows.empty:
            result = rows[["tissue", "tpm"]].to_dict(orient="records")
            fallback_used = True
            cache.set(cache_key, result)
            return result
    except Exception as e:
        logging.debug(f"GTEx fallback load error: {e}")
    fallback_used = True
    return None

@retry_on_failure(retries=3, backoff=2.0)
def fetch_genomic_sequence(chrom: str, start: int, end: int, force_live: bool = False) -> str:
    """Fetch raw DNA sequence from UCSC with caching and fallback (empty string)."""
    global fallback_used
    cache_key = f"sequence:{chrom}:{start}:{end}"
    cached = cache.get(cache_key)
    if cached is not None:
        logging.debug(f"Cache hit for {cache_key}")
        return cached

    if not chrom.startswith("chr"):
        chrom = f"chr{chrom}"
    url = f"{UCSC_API}/getData/sequence?genome=hg38&chrom={chrom}&start={start}&end={end}"
    try:
        resp = requests.get(url, timeout=15, headers=HEADERS)
        if resp.status_code == 200:
            data = resp.json()
            seq = data.get("dna", "").upper()
            cache.set(cache_key, seq)
            return seq
        else:
            logging.debug(f">> Warning: UCSC API failed: {resp.status_code}")
    except Exception as e:
        logging.debug(f">> Error fetching sequence: {e}")
    fallback_used = True
    return ""

@retry_on_failure(retries=3, backoff=2.0)
def fetch_gene_structure(chrom: str, pos: int, window_size: int = 100, force_live: bool = False) -> List[dict]:
    """Fetch exon coordinates from Ensembl with caching and JSON fallback."""
    global fallback_used
    cache_key = f"gene_structure:{chrom}:{pos}:{window_size}"
    cached = cache.get(cache_key)
    if cached is not None:
        logging.debug(f"Cache hit for {cache_key}")
        return cached

    clean_chrom = chrom.replace("chr", "")
    start = pos - window_size
    end = pos + window_size
    url = f"{ENSEMBL_API}/overlap/region/human/{clean_chrom}:{start}-{end}?feature=exon;content-type=application/json"
    try:
        resp = SESSION.get(url, timeout=15)
        if resp.status_code == 200:
            features = resp.json()
            exons = []
            for f in features:
                rel_start = f["start"] - pos
                rel_end = f["end"] - pos
                exons.append({
                    "start": rel_start,
                    "end": rel_end,
                    "id": f.get("id"),
                    "gene_id": f.get("Parent"),
                })
            cache.set(cache_key, exons)
            return exons
    except Exception as e:
        logging.debug(f"Error fetching gene structure: {e}")

    if force_live:
        logging.warning("Force live Gene Structure failed; not attempting fallback.")
        fallback_used = True
        return []

    # Fallback JSON (demo data)
    fallback_path = os.path.join(FALLBACK_DIR, "gene_structure_fallback.json")
    if os.path.exists(fallback_path):
        try:
            with open(fallback_path, "r", encoding="utf-8") as f:
                cache_data = json.load(f)
            demo_genes = {"chr22": "MYH9", "chr1": "PCSK9", "chr2": "APOB", "chr3": "MYL3"}
            c_key = chrom if chrom.startswith("chr") else f"chr{chrom}"
            gene_sym = demo_genes.get(c_key)
            if gene_sym and gene_sym in cache_data:
                result = cache_data[gene_sym]
                fallback_used = True
                cache.set(cache_key, result)
                return result
        except Exception as e:
            logging.debug(f"Gene structure fallback load error: {e}")
    fallback_used = True
    return []

def load_fallback_related_data(chrom: str, pos: int, ref: str, alt: str) -> List[Dict[str, Any]]:
    """Load related data from fallback JSON if available."""
    try:
        path = os.path.join(FALLBACK_DIR, "related_data_fallback.json")
        if not os.path.exists(path):
            return []
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
        # Simple filter
        results = []
        for item in data:
            if (item.get("chrom") == chrom and 
                item.get("pos") == pos and 
                item.get("ref") == ref and 
                item.get("alt") == alt):
                results.append(item)
        return results
    except Exception as e:
        logging.debug(f"Related data fallback load error: {e}")
        return []

# ---------------------------------------------------------------------------
# Rate Limiting
# ---------------------------------------------------------------------------
def rate_limit(calls: int = 1, period: float = 1.0):
    """Simple rate limiter decorator.
    Ensures that the decorated function is not called more than `calls` times in `period` seconds.
    """
    min_interval = period / calls
    last_call_time = 0.0

    def decorator(func):
        def wrapper(*args, **kwargs):
            nonlocal last_call_time
            elapsed = time.time() - last_call_time
            if elapsed < min_interval:
                time.sleep(min_interval - elapsed)
            result = func(*args, **kwargs)
            last_call_time = time.time()
            return result
        return wrapper
    return decorator

# ---------------------------------------------------------------------------
# Backup APIs (MyGene.info / MyVariant.info)
# ---------------------------------------------------------------------------
def fetch_mygene_info(gene_symbol):
    """
    Fetch gene metadata from MyGene.info (v3).
    Acts as a backup when Ensembl is unavailable.
    """
    url = f"https://mygene.info/v3/query?q=symbol:{gene_symbol}&species=human&fields=symbol,name,summary,genomic_pos,type_of_gene,alias"
    try:
        resp = requests.get(url, timeout=5)
        if resp.status_code == 200:
            data = resp.json()
            if data.get("hits"):
                hit = data["hits"][0]
                # Transform to match Ensembl format
                return {
                    "id": hit.get("_id"),
                    "display_name": hit.get("symbol"),
                    "description": hit.get("name"),
                    "seq_region_name": hit.get("genomic_pos", {}).get("chr"),
                    "start": hit.get("genomic_pos", {}).get("start"),
                    "end": hit.get("genomic_pos", {}).get("end"),
                    "biotype": hit.get("type_of_gene"),
                    "summary": hit.get("summary"),
                    "source": "MyGene.info"
                }
    except Exception as e:
        print(f"MyGene.info fetch failed: {e}")
    return None

def fetch_myvariant_info(chrom, pos, ref, alt):
    """
    Fetch variant annotations from MyVariant.info (v1).
    Acts as a backup for ClinVar and gnomAD frequencies.
    Uses hg38 coordinates via the query endpoint.
    """
    # Remove 'chr' prefix for query if present (MyVariant expects '1', 'X', etc. for hg38.chr)
    c = chrom.replace("chr", "")
    
    # Construct query for hg38
    # q=hg38.chr:22 AND hg38.start:36305975 AND vcf.ref:G AND vcf.alt:A
    query = f"hg38.chr:{c} AND hg38.start:{pos} AND vcf.ref:{ref} AND vcf.alt:{alt}"
    
    url = "https://myvariant.info/v1/query"
    params = {
        "q": query,
        "fields": "clinvar,gnomad_exome,gnomad_genome,dbsnp",
        "size": 1
    }
    
    try:
        resp = requests.get(url, params=params, timeout=5)
        if resp.status_code == 200:
            data = resp.json()
            if data.get("hits"):
                return data["hits"][0]
            else:
                # Debug: print if no hits found
                # print(f"MyVariant.info: No hits for {query}")
                return None
        else:
            print(f"MyVariant.info query failed: {resp.status_code} for {resp.url}")
    except Exception as e:
        print(f"MyVariant.info fetch failed: {e}")
    return None


def fetch_single_cell_expression(gene_symbol: str) -> Optional[List[Dict[str, Any]]]:
    """Fetch single-cell expression data from local curated dataset."""
    global fallback_used
    cache_key = f"single_cell:{gene_symbol}"
    cached = cache.get(cache_key)
    if cached is not None:
        logging.debug(f"Cache hit for {cache_key}")
        return cached

    # Load from JSON
    try:
        path = os.path.join(FALLBACK_DIR, "single_cell_data.json")
        if not os.path.exists(path):
            logging.warning(f"Single cell data file not found at {path}")
            return None
            
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
            
        for entry in data:
            if entry.get("symbol") == gene_symbol:
                expr_data = entry.get("expression", [])
                cache.set(cache_key, expr_data)
                return expr_data
                
    except Exception as e:
        logging.debug(f"Single cell data load error: {e}")
        
    return None

# ---------------------------------------------------------------------------
# New API Wrappers (ClinVar & dbSNP)
# ---------------------------------------------------------------------------

@rate_limit(calls=3, period=1.0)
@retry_on_failure(retries=3, backoff=2.0)
def fetch_clinvar_variants(chrom: str, pos: int, ref: str, alt: str, force_live: bool = False) -> Optional[Dict[str, Any]]:
    """Fetch ClinVar data using NCBI ClinVar Variation API (better than E-utilities)."""
    global fallback_used
    cache_key = f"clinvar:{chrom}:{pos}:{ref}:{alt}"
    cached = cache.get(cache_key)
    if cached is not None:
        logging.debug(f"Cache hit for {cache_key}")
        return cached

    clean_chrom = chrom.replace("chr", "")
    
    # Use NCBI ClinVar Variation Viewer API
    # This provides structured JSON data for variants
    try:
        # Build SPDI notation: NC_000001.11:12345:A:T
        # Map chromosome to RefSeq accession
        chrom_to_refseq = {
            "1": "NC_000001.11", "2": "NC_000002.12", "3": "NC_000003.12",
            "4": "NC_000004.12", "5": "NC_000005.10", "6": "NC_000006.12",
            "7": "NC_000007.14", "8": "NC_000008.11", "9": "NC_000009.12",
            "10": "NC_000010.11", "11": "NC_000011.10", "12": "NC_000012.12",
            "13": "NC_000013.11", "14": "NC_000014.9", "15": "NC_000015.10",
            "16": "NC_000016.10", "17": "NC_000017.11", "18": "NC_000018.10",
            "19": "NC_000019.10", "20": "NC_000020.11", "21": "NC_000021.9",
            "22": "NC_000022.11", "X": "NC_000023.11", "Y": "NC_000024.10"
        }
        
        refseq_acc = chrom_to_refseq.get(clean_chrom)
        if not refseq_acc:
            logging.debug(f"No RefSeq accession for chromosome {clean_chrom}")
            return None
        
        # Try ClinVar API with genomic coordinates (GRCh38/hg38)
        # Use esearch to find variants at this position
        search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        search_term = f"{clean_chrom}[CHR] AND {pos}[CHRPOS38]"  # Use GRCh38 coordinates
        
        search_params = {
            "db": "clinvar",
            "term": search_term,
            "retmode": "json",
            "retmax": 10
        }
        
        logging.info(f"ClinVar search: {search_term}")
        search_resp = SESSION.get(search_url, params=search_params, timeout=15)
        
        if search_resp.status_code != 200:
            logging.warning(f"ClinVar search failed: {search_resp.status_code}")
            raise RuntimeError(f"ClinVar search failed: {search_resp.status_code}")
        
        search_data = search_resp.json()
        id_list = search_data.get("esearchresult", {}).get("idlist", [])
        
        if not id_list:
            logging.info(f"No ClinVar entries found for {clean_chrom}:{pos} (GRCh38)")
            return None
        
        # Get detailed information using esummary
        summary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        summary_params = {
            "db": "clinvar",
            "id": ",".join(id_list[:5]),  # Limit to first 5
            "retmode": "json"
        }
        
        summary_resp = SESSION.get(summary_url, params=summary_params, timeout=15)
        
        if summary_resp.status_code != 200:
            raise RuntimeError(f"ClinVar summary failed: {summary_resp.status_code}")
        
        summary_data = summary_resp.json()
        result_data = summary_data.get("result", {})
        
        # Find the best matching variant
        best_match = None
        for uid in id_list[:5]:
            item = result_data.get(uid)
            if not item:
                continue
            
            # Check if this variant matches our ref/alt
            # ClinVar data structure varies, so we'll take the first reasonable match
            variation_set = item.get("variation_set", [])
            if variation_set:
                # Found a match, extract data
                best_match = {
                    "variation_id": uid,
                    "title": item.get("title", ""),
                    "clinical_significance": item.get("clinical_significance", {}).get("description", ""),
                    "review_status": item.get("clinical_significance", {}).get("review_status", ""),
                    "last_evaluated": item.get("clinical_significance", {}).get("last_evaluated", ""),
                    "germline_classification": item.get("germline_classification", {}).get("description", ""),
                    "variation_type": item.get("variation_type", ""),
                    "molecular_consequence": item.get("molecular_consequence", []),
                    "gene_symbol": item.get("genes", [{}])[0].get("symbol") if item.get("genes") else None,
                    "accession": item.get("accession", ""),
                    "source": "NCBI ClinVar"
                }
                break
        
        if best_match:
            cache.set(cache_key, best_match)
            return best_match
        else:
            logging.debug(f"No matching ClinVar variant found for {clean_chrom}:{pos}:{ref}>{alt}")
            return None
            
    except Exception as e:
        logging.debug(f"ClinVar API error: {e}")

    if force_live:
        logging.warning("Force live ClinVar failed; not attempting fallback.")
        fallback_used = True
        return None

    # Fallback to local data if available
    fallback_used = True
    return None

@rate_limit(calls=3, period=1.0)
@retry_on_failure(retries=3, backoff=2.0)
def fetch_dbsnp_variants(chrom: str, pos: int, ref: str, alt: str, force_live: bool = False) -> Optional[Dict[str, Any]]:
    """Fetch dbSNP rsID and info using NCBI E-utilities."""
    global fallback_used
    cache_key = f"dbsnp:{chrom}:{pos}:{ref}:{alt}"
    cached = cache.get(cache_key)
    if cached is not None:
        logging.debug(f"Cache hit for {cache_key}")
        return cached

    clean_chrom = chrom.replace("chr", "")
    term = f"{clean_chrom}[CHR] AND {pos}[POS]"
    
    try:
        # 1. Search
        search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        params = {
            "db": "snp",
            "term": term,
            "retmode": "json",
            "retmax": 5
        }
        resp = SESSION.get(search_url, params=params, timeout=10)
        if resp.status_code == 200:
            search_data = resp.json()
            id_list = search_data.get("esearchresult", {}).get("idlist", [])
            
            if id_list:
                # 2. Summary
                summary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
                sum_params = {
                    "db": "snp",
                    "id": ",".join(id_list),
                    "retmode": "json"
                }
                sum_resp = SESSION.get(summary_url, params=sum_params, timeout=10)
                if sum_resp.status_code == 200:
                    sum_data = sum_resp.json()
                    result = sum_data.get("result", {})
                    for uid in id_list:
                        item = result.get(uid)
                        if item:
                            data = {
                                "rsid": f"rs{uid}",
                                "global_maf": item.get("global_maf"),
                                "clinical_significance": item.get("clinical_significance"),
                                "gene": item.get("genes", [{}])[0].get("name") if item.get("genes") else None
                            }
                            cache.set(cache_key, data)
                            return data
    except Exception as e:
        logging.debug(f"dbSNP API error: {e}")

    if force_live:
        logging.warning("Force live dbSNP failed; not attempting fallback.")
        fallback_used = True
        return None

    fallback_used = True
    return None
