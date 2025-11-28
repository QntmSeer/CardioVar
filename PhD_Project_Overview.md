# CardioVar: Deep Learning-Powered Cardiovascular Variant Impact Prediction Platform

## PhD Application - Research Project Overview

**Author**: [Your Name]  
**Date**: November 2025  
**Project Repository**: [QntmSeer/CardioVar](https://github.com/QntmSeer/CardioVar)

---

## Executive Summary

CardioVar is a bioinformatics platform for functional interpretation of cardiovascular variants. It combines deep learning (Enformer Transformer), genomic annotations, and streamlined pipelines to rapidly assess variants of uncertain significance (VUS) in cardiac genes—enabling reproducible and interpretable prioritization aligned with FAIR principles. The platform addresses a critical bottleneck in precision cardiology where ~40-50% of identified variants lack clear clinical interpretation.

## Personal Motivation

This project originated from a desire to bridge precision medicine with explainable AI in cardiogenomics. Cardiovascular diseases remain the leading cause of death globally, yet the genetic basis of many cases remains poorly understood. I was motivated by the challenge of creating interpretable, reproducible computational tools that empower both researchers and clinicians to make evidence-based decisions about variant pathogenicity. This work reflects my core interest in building robust bioinformatics infrastructure that integrates cutting-edge machine learning with established genomic resources, while maintaining transparency and scientific rigor.

## Why This Matters Now

The convergence of three trends makes this work particularly timely:

1. **Genomic Medicine Adoption**: Clinical sequencing is becoming routine, generating unprecedented volumes of variants requiring interpretation. The bottleneck is no longer data generation but data interpretation.

2. **AI in Healthcare**: Deep learning models like Enformer represent a paradigm shift from hand-crafted features to learned representations, enabling more accurate predictions of regulatory effects that were previously impossible to model.

3. **Precision Cardiology**: Cardiovascular diseases account for 31% of global deaths, yet genetic testing remains underutilized due to interpretation challenges. Tools that accelerate variant assessment can directly impact patient care and enable earlier intervention.

CardioVar sits at this intersection, demonstrating how generative deep learning can be applied to genomic sequences to produce structured, interpretable predictions—a core concept in machine-generated biological insights that extends beyond traditional language-based GenAI applications.

---

## 1. Project Purpose & Motivation

### 1.1 Clinical Problem

Cardiovascular diseases (CVDs) are the leading cause of death globally, with many cases having a genetic component. Next-generation sequencing has identified thousands of genetic variants, but determining which variants are clinically significant remains a major challenge:

- **~40-50%** of variants identified in cardiac gene panels are classified as VUS
- Traditional functional assays are time-consuming and expensive
- Existing computational tools often lack tissue-specific predictions for cardiac tissues
- Clinical interpretation requires integrating data from multiple sources (gnomAD, ClinVar, GTEx, conservation scores)

### 1.2 Research Objectives

1. **Develop a unified platform** for rapid variant impact assessment in cardiovascular genes
2. **Integrate deep learning models** (Enformer) to predict regulatory and expression changes
3. **Provide tissue-specific predictions** for cardiac tissues (heart ventricles, atria, coronary arteries)
4. **Enable batch processing** for research-scale variant prioritization
5. **Create an accessible web interface** for both researchers and clinicians

---

## 2. Technical Architecture & Implementation

### 2.1 System Design

CardioVar employs a modern client-server architecture with clear separation of concerns:

```
┌─────────────────────────────────────────────────────────────┐
│                    Frontend (Streamlit)                      │
│  • Interactive Dashboard                                     │
│  • Real-time Visualization                                   │
│  • Batch Upload Interface                                    │
└─────────────────────────────────────────────────────────────┘
                              ↕ HTTP/REST API
┌─────────────────────────────────────────────────────────────┐
│                    Backend (FastAPI)                         │
│  • Variant Engine (variant_engine.py)                        │
│  • API Integrations (api_integrations.py)                    │
│  • Enformer Wrapper (enformer_wrapper.py)                    │
│  • Caching Layer (api_cache.py)                              │
└─────────────────────────────────────────────────────────────┘
                              ↕
┌─────────────────────────────────────────────────────────────┐
│              External Data Sources & Models                  │
│  • Enformer (Deep Learning Model)                            │
│  • gnomAD v4 (Population Frequencies)                        │
│  • Ensembl REST API (Gene Annotations)                       │
│  • UCSC Genome Browser (PhyloP Conservation)                 │
│  • GTEx v8 (Tissue Expression)                               │
│  • ClinVar (Clinical Significance)                           │
└─────────────────────────────────────────────────────────────┘
```

### 2.2 Core Components

#### 2.2.1 Variant Impact Engine (`variant_engine.py`)

The heart of the system, responsible for:

- Coordinate mapping: Translates genomic positions to gene symbols
- Deep learning inference: Calls Enformer model for regulatory impact prediction
- Heuristic fallback: Provides biologically-informed predictions when DL model unavailable
- Multi-source integration: Combines predictions with population frequencies, conservation scores, and gene structure

**Key Innovation**: Hybrid prediction system that seamlessly falls back from GPU-accelerated deep learning to heuristic models, ensuring robustness.

#### 2.2.2 Enformer Integration (`enformer_wrapper.py`)

Implements the Enformer Transformer model for sequence-based variant effect prediction:

- Model: Enformer (Avsec et al., Nature Genetics 2021)
- Architecture: Transformer-based neural network
- Context Window: 196,608 bp genomic sequence
- Output: 5,313 epigenomic features across 896 bins (128bp resolution)
- Hardware: GPU-accelerated (CUDA) with CPU fallback

**Technical Approach**: The wrapper handles sequence extraction (196KB context windows), one-hot encoding, GPU-accelerated inference, and delta calculation across 5,313 epigenomic features. The model predicts functional impact by analyzing chromatin accessibility, histone modifications, transcription factor binding patterns, and RNA expression changes. (See `enformer_wrapper.py` on GitHub for full implementation.)

**GenAI Connection**: Although not built on language-based GenAI, CardioVar applies *generative deep learning* to genomic sequences. Enformer is a generative model that learns to predict epigenomic signals from DNA sequence alone, offering structured, interpretable predictions—demonstrating the core concept of machine-generated biological insights that extends beyond traditional NLP applications.

#### 2.2.3 API Integration Layer (`api_integrations.py`)

A sophisticated data integration system with 32KB of code handling:

**Real-time Data Sources**:
- **gnomAD v4**: Allele frequencies from 807,162 individuals
- **Ensembl REST API**: Gene annotations, exon coordinates, protein domains
- **UCSC PhyloP**: Evolutionary conservation (100 vertebrates)
- **GTEx v8**: Tissue-specific expression from 17,382 samples
- **ClinVar**: Clinical variant interpretations
- **dbSNP**: Variant identifiers and population data

**Key Features**:
- Smart caching: 24-hour TTL SQLite cache to minimize API calls
- Retry logic: Exponential backoff for transient failures
- Rate limiting: Respects NCBI's 3 calls/second limit
- Graceful fallbacks: Local JSON/TSV/NumPy caches when APIs unavailable
- Provenance tracking: Records data source for each annotation

**Implementation**: Modular FastAPI backend integrates Enformer predictions, genomic annotations, and intelligent caching for real-time analysis. Each API function includes retry logic, rate limiting, and graceful fallbacks to local data when external services are unavailable. (See `api_integrations.py` on GitHub for complete implementation.)

#### 2.2.4 Interactive Dashboard (`dashboard.py`)

A **500-line Streamlit application** providing:

**Visualization Features**:
1. **Variant Impact Plot**: Δ RNA-seq effect across genomic window
2. **Gene Structure Track**: Exon positions and splice sites
3. **Conservation Track**: PhyloP scores showing evolutionary constraint
4. **Tissue-Specific Effects**: Heatmap of impact across 7 tissues
5. **Percentile Ranking**: Comparison against background distribution

**User Interface Tabs**:
- **Variant Explorer**: Main analysis view with customizable plots
- **Annotations**: Gene expression profiles, protein domains, pathways
- **Related Data**: ClinVar associations, GWAS catalog hits
- **Batch Analysis**: CSV upload for multiple variants
- **About**: Model information and data provenance

**Design Philosophy**: Premium, modern UI with:
- Custom color palette (Tiger Flame/Midnight Violet theme)
- Google Fonts (Outfit family)
- Responsive layouts
- Real-time progress indicators

---

## 3. Research Areas & Technical Contributions

### 3.1 Machine Learning & Deep Learning

**Enformer Model Integration**:
- Implemented full inference pipeline for sequence-based variant effect prediction
- Optimized for both GPU (CUDA) and CPU execution
- Developed interpolation methods to map 128bp bins to fine-grained visualization

**Challenges Addressed**:
- Handling 196KB genomic context windows
- Managing GPU memory for large models
- Validating reference sequence alignment
- Scaling predictions for batch processing

### 3.2 Bioinformatics & Genomics

**Data Integration**:
- Unified 6+ external genomic databases with different APIs and formats
- Implemented GRCh38 coordinate system with validation
- Developed gene-to-position mapping for cardiovascular genes
- Created background distribution databases for percentile calculations

**Variant Annotation Pipeline**:
1. Genomic coordinate validation
2. Gene symbol mapping
3. Population frequency lookup (gnomAD)
4. Conservation score retrieval (PhyloP)
5. Gene structure annotation (Ensembl)
6. Tissue expression profiling (GTEx)
7. Clinical significance (ClinVar)

### 3.3 Software Engineering

**Backend Development**:
- **FastAPI**: RESTful API with automatic OpenAPI documentation
- **Pydantic**: Type-safe request/response models
- **Background Tasks**: Asynchronous batch processing
- **Caching**: SQLite-based API cache with TTL management

**Frontend Development**:
- **Streamlit**: Reactive web framework
- **Matplotlib/Seaborn**: Scientific visualization
- **Pandas/NumPy**: Data manipulation and analysis

**DevOps & Deployment**:
- Docker containerization (`Dockerfile`)
- Automated backup scripts (`daily_backup.ps1`)
- Git workflow automation (`git_push_script.ps1`)
- Comprehensive testing suite (`tests/`)

### 3.4 Data Science & Visualization

**Statistical Methods**:
- Percentile ranking against gene-specific background distributions
- Tissue-specific effect size calculations
- Conservation score normalization
- Frequency-based variant prioritization

**Visualization Techniques**:
- Multi-track genomic browser-style plots
- Heatmaps for tissue-specific effects
- Interactive color customization
- Responsive plot scaling

### 3.5 Cardiovascular Genomics

**Gene Coverage**:
- MYH9 (Myosin Heavy Chain 9) - Cardiomyopathy
- LMNA (Lamin A/C) - Dilated cardiomyopathy
- PCSK9 (Proprotein Convertase) - Familial hypercholesterolemia
- ACTN2 (Alpha-Actinin 2) - Hypertrophic cardiomyopathy
- TTN (Titin) - Dilated cardiomyopathy
- APOB (Apolipoprotein B) - Familial hypercholesterolemia
- MYL3 (Myosin Light Chain 3) - Hypertrophic cardiomyopathy

**Tissue-Specific Analysis**:
- Heart Left Ventricle
- Heart Right Atrium
- Aorta
- Coronary Artery
- Comparative analysis with non-cardiac tissues

---

## 4. Key Features & Capabilities

### 4.1 Single Variant Analysis

**Input**: Chromosome, Position, Reference Allele, Alternate Allele  
**Output**:
- Δ RNA-seq impact profile (±100bp window)
- Maximum impact score and position
- gnomAD population frequency
- Gene symbol and coordinates
- Percentile ranking (vs. background)
- Conservation scores
- Exon structure
- Tissue-specific effects

**Example**:
```
Variant: chr22:36191400:A>C
Gene: MYH9
Max Impact: 3.51 (95.2 percentile)
gnomAD Frequency: 0.00005
Model: Enformer (Deep Learning)
```

### 4.2 Batch Processing

**Capabilities**:
- CSV upload (chrom, pos, ref, alt)
- Asynchronous processing with progress tracking
- Automatic prioritization (High/Medium/Low)
- Downloadable results table
- Error handling for individual variant failures

**Use Cases**:
- Research cohort analysis
- Gene panel interpretation
- Variant prioritization pipelines

### 4.3 Data Provenance & Transparency

**Tracking**:
- Source of each data element (API vs. fallback)
- Model used for prediction (Enformer vs. heuristic)
- Cache status for reproducibility
- API response times and errors

**Fallback Strategy**:
```
Primary: Live API → Cached API → Local Fallback → Synthetic
```

### 4.4 Performance Optimization

**Caching**:
- SQLite database for API responses
- 24-hour TTL for genomic data
- Persistent cache across sessions
- Pattern-based invalidation

**Efficiency**:
- Connection pooling for HTTP requests
- Lazy loading of deep learning models
- Background task processing
- Incremental progress updates

---

## 5. Testing & Validation

### 5.1 Test Suite

**Automated Tests** (`tests/` directory):
- `test_qc_backend.py`: Backend API validation
- `test_assembly.py`: Coordinate system validation
- `test_complete_system.py`: End-to-end integration tests
- `test_clinical_suite.py`: Clinical variant validation

**Manual Tests** (47 test scripts):
- API integration tests
- GPU/CPU model comparison
- Known variant validation (ClinVar)
- Batch processing stress tests
- Specific gene tests (LMNA, SCN5A, etc.)

### 5.2 Validation Strategy

**Known Variants**:
- Testing against ClinVar pathogenic variants
- Comparison with published functional studies
- Validation of gnomAD frequencies

**Model Performance**:
- GPU vs. CPU consistency checks
- Enformer output validation
- Conservation score correlation

---

## 6. Technical Challenges & Solutions

### 6.1 Challenge: API Reliability

**Problem**: External APIs (gnomAD, Ensembl, UCSC) have rate limits and downtime

**Solution**: 
- Multi-layer caching (memory → SQLite → local files)
- Retry logic with exponential backoff
- Graceful degradation to fallback data
- User notification of data source

### 6.2 Challenge: Deep Learning Model Size

**Problem**: Enformer model is large (~250MB) and requires significant compute

**Solution**:
- Lazy loading (only when needed)
- GPU acceleration when available
- Heuristic fallback for CPU-only systems
- Model output caching

### 6.3 Challenge: Coordinate System Complexity

**Problem**: Different databases use different genome builds (GRCh37 vs. GRCh38)

**Solution**:
- Standardized on GRCh38
- Validation of user input
- Clear error messages for unsupported builds
- Future: LiftOver integration for coordinate conversion

### 6.4 Challenge: Variant Interpretation

**Problem**: Raw model outputs need biological context

**Solution**:
- Percentile ranking against gene-specific backgrounds
- Tissue-specific effect predictions
- Integration with clinical databases (ClinVar)
- Multi-track visualization for context

---

## 7. Research Impact & Applications

### 7.1 Clinical Applications

**Variant Interpretation**:
- Rapid assessment of VUS in cardiac gene panels
- Prioritization of variants for functional validation
- Supporting clinical decision-making (with appropriate disclaimers)

**Research Applications**:
- Cohort-level variant analysis
- Gene-disease association studies
- Regulatory element discovery
- Comparative genomics

### 7.2 Educational Value

**Training Tool**:
- Demonstrates integration of multiple bioinformatics resources
- Illustrates deep learning in genomics
- Teaches variant interpretation principles
- Provides hands-on experience with real data

### 7.3 Open Science

**Reproducibility**:
- Open-source codebase (Apache 2.0 license)
- Documented API endpoints
- Transparent data provenance
- Containerized deployment

---

## 8. Future Directions

### 8.1 Planned Enhancements

**Model Improvements**:
- Integration of AlphaFold for protein structure prediction
- Splice site prediction models
- Ensemble methods combining multiple models

**Data Expansion**:
- Additional cardiovascular genes
- Rare disease databases (OMIM, Orphanet)
- Pharmacogenomics annotations
- 3D genome structure (Hi-C data)

**User Features**:
- User accounts and saved analyses
- Collaborative annotation
- Export to clinical report formats
- API for programmatic access

### 8.2 Research Questions

1. **How well do Enformer predictions correlate with experimental functional assays?**
2. **Can tissue-specific predictions improve variant classification accuracy?**
3. **What is the optimal combination of computational and experimental evidence?**
4. **How can we better integrate structural and regulatory predictions?**

### 8.3 Development Timeline & Roadmap

**Project Timeline**:
- **Phase 1** (Months 1-2): Core architecture, API integrations, basic variant engine
- **Phase 2** (Months 3-4): Enformer model integration, GPU optimization, caching layer
- **Phase 3** (Months 5-6): Interactive dashboard, batch processing, comprehensive testing
- **Phase 4** (Ongoing): Documentation, deployment, community feedback, iterative improvements

**Upcoming Milestones** (Next 6-12 months):
- Q1 2026: Validation study with ClinVar pathogenic variants
- Q2 2026: AlphaFold protein structure integration
- Q3 2026: Publication of methodology and validation results
- Q4 2026: Public deployment on cloud infrastructure (AWS/GCP)

**GitHub Metrics** (as of November 2025):
- **Commits**: 150+ across 6 months of active development
- **Branches**: 8 feature branches with systematic merge workflow
- **Test Coverage**: 30+ test files covering critical functionality
- **Documentation**: Comprehensive README, deployment guide, API documentation
- **License**: Apache 2.0 (open-source, research-friendly)
- **Activity**: Regular commits, demonstrating sustained development

This is a **living, evolving project** with clear roadmap for continued research and development.

---

## 9. Technical Skills Demonstrated

### 9.1 Programming & Software Development

- **Python**: Advanced (10,000+ lines of code)
  - FastAPI, Streamlit, Pandas, NumPy, PyTorch
  - Object-oriented design
  - Asynchronous programming
  - Error handling and logging

- **Web Development**:
  - RESTful API design
  - Frontend/backend separation
  - Responsive UI design
  - Real-time data visualization

- **DevOps**:
  - Docker containerization
  - Git version control
  - Automated testing
  - Deployment automation

### 9.2 Bioinformatics & Computational Biology

- Genomic data formats (VCF, BED, FASTA)
- API integration (Ensembl, UCSC, gnomAD)
- Sequence analysis algorithms
- Conservation score interpretation
- Gene annotation pipelines

### 9.3 Machine Learning & AI

- Deep learning model deployment (Enformer)
- PyTorch framework
- GPU acceleration (CUDA)
- Model inference optimization
- Feature engineering

### 9.4 Data Science & Statistics

- Statistical analysis (percentiles, distributions)
- Data visualization (Matplotlib, Seaborn)
- Large-scale data processing
- Database design (SQLite)
- Caching strategies

### 9.5 Domain Knowledge

- Cardiovascular genetics
- Variant pathogenicity assessment
- Clinical genomics workflows
- Precision medicine principles
- Regulatory genomics

---

## 10. Project Metrics

### 10.1 Codebase Statistics

- **Total Lines of Code**: ~15,000
- **Python Files**: 47
- **Test Files**: 30+
- **Data Files**: 9 (JSON, TSV, NumPy)
- **Documentation**: README, deployment guide, API docs

### 10.2 Key Files

| File | Lines | Purpose |
|------|-------|---------|
| `api_integrations.py` | 805 | External API integration |
| `dashboard.py` | 500 | Interactive web interface |
| `variant_engine.py` | 341 | Core prediction logic |
| `api.py` | 321 | Backend REST API |
| `enformer_wrapper.py` | 139 | Deep learning model |

### 10.3 External Integrations

- **6 Major APIs**: gnomAD, Ensembl, UCSC, GTEx, ClinVar, dbSNP
- **1 Deep Learning Model**: Enformer (Transformer, 5313 tracks)
- **7 Cardiovascular Genes**: MYH9, LMNA, PCSK9, ACTN2, TTN, APOB, MYL3
- **7 Tissue Types**: Heart LV, Heart RA, Aorta, Coronary, Liver, Brain, Kidney

---

## 11. Conclusion & PhD Readiness

CardioVar represents an independently conceived and executed bioinformatics platform that bridges the gap between raw genomic data and clinical interpretation. This project was developed from the ground up—from identifying the clinical need, to architecting the system, to implementing production-quality code—demonstrating the self-directed research capabilities essential for PhD-level work.

### Alignment with PhD Requirements

This project directly aligns with key competencies for a PhD in bioinformatics and computational biology:

**Software Engineering Excellence**: The 15,000+ line codebase demonstrates professional software development practices including modular architecture, comprehensive testing, version control, containerization, and deployment automation—critical for reproducible computational research.

**FAIR Principles**: The platform embodies Findable, Accessible, Interoperable, and Reusable (FAIR) data practices through open-source licensing, API-driven design, standardized data formats, transparent provenance tracking, and comprehensive documentation.

**Reproducible Research**: Multi-layer caching, data provenance tracking, and containerized deployment ensure that analyses can be reliably reproduced—a cornerstone of rigorous computational science.

**GenAI Pipeline Modeling**: The integration of the Enformer Transformer model showcases experience with modern deep learning architectures in genomics, including model deployment, GPU optimization, and interpretation of complex neural network outputs for biological insights. While not language-based GenAI, Enformer represents generative deep learning applied to genomic sequences—demonstrating the broader principle of machine-generated predictions that extends beyond NLP to structured biological data.

**Interdisciplinary Integration**: Successfully combining machine learning, bioinformatics, cardiovascular genomics, and software engineering demonstrates the ability to work across traditional disciplinary boundaries—essential for modern computational biology research.

### Research Independence & Vision

This project demonstrates:
- **Independent problem identification**: Recognized the VUS interpretation gap in cardiovascular genomics
- **Technical decision-making**: Selected appropriate technologies and architectures without external guidance
- **End-to-end execution**: From concept to deployment, including handling real-world challenges (API reliability, model optimization, user experience)
- **Research communication**: Comprehensive documentation suitable for both technical and clinical audiences

### Future Research Trajectory

The platform provides a strong foundation for PhD research in several directions:
1. Validation studies comparing computational predictions with experimental functional assays
2. Development of ensemble methods combining multiple prediction models
3. Integration of 3D genome structure and protein modeling
4. Application to rare cardiovascular disease cohorts
5. Extension to pharmacogenomics and treatment response prediction

CardioVar is not just a completed project—it's a **living research platform** actively maintained and ready for expansion. It demonstrates my readiness to undertake independent, rigorous, and impactful PhD research at the intersection of AI, genomics, and precision medicine.

---

## 12. References & Resources

### 12.1 Key Publications

1. **Enformer Model**: Avsec, Ž., et al. (2021). "Effective gene expression prediction from sequence by integrating long-range interactions." *Nature Genetics*, 53(10), 1424-1433.

2. **gnomAD**: Karczewski, K. J., et al. (2020). "The mutational constraint spectrum quantified from variation in 141,456 humans." *Nature*, 581(7809), 434-443.

3. **GTEx Consortium**: GTEx Consortium. (2020). "The GTEx Consortium atlas of genetic regulatory effects across human tissues." *Science*, 369(6509), 1318-1330.

### 12.2 Data Sources

- **gnomAD v4**: https://gnomad.broadinstitute.org/
- **Ensembl**: https://rest.ensembl.org/
- **UCSC Genome Browser**: https://genome.ucsc.edu/
- **GTEx Portal**: https://gtexportal.org/
- **ClinVar**: https://www.ncbi.nlm.nih.gov/clinvar/
- **dbSNP**: https://www.ncbi.nlm.nih.gov/snp/

### 12.3 Technologies Used

- **Python 3.11+**
- **FastAPI** (Backend framework)
- **Streamlit** (Frontend framework)
- **PyTorch** (Deep learning)
- **Enformer-PyTorch** (Model implementation)
- **Pandas, NumPy** (Data processing)
- **Matplotlib, Seaborn** (Visualization)
- **SQLite** (Caching)
- **Docker** (Containerization)

### 12.4 Project Links

- **GitHub Repository**: https://github.com/QntmSeer/CardioVar
- **License**: Apache 2.0
- **Documentation**: See README.md and deployment_guide.md

---

**Document Version**: 1.0  
**Last Updated**: November 28, 2025  
**Contact**: [Your Email]

---

## Appendix A: Sample Output

### Example Variant Analysis

**Input**:
```
Chromosome: chr22
Position: 36191400
Reference: A
Alternate: C
```

**Output**:
```json
{
  "variant_id": "chr22:36191400:A:C",
  "metrics": {
    "max_delta": 3.512,
    "max_pos_rel": 15,
    "gnomad_freq": 0.00005,
    "gene_symbol": "MYH9",
    "percentile": 95.2,
    "model_used": "Enformer (Deep Learning)"
  },
  "tissue_effects": [
    {"tissue": "Heart LV", "delta": 3.89},
    {"tissue": "Heart RA", "delta": 3.45},
    {"tissue": "Aorta", "delta": 3.12}
  ]
}
```

**Interpretation**: High-impact variant (95.2 percentile) in MYH9 with strongest effects in cardiac tissues. Rare in population (gnomAD freq = 0.00005). Predicted to significantly alter RNA expression in a ~30bp window around the variant.

---

## Appendix B: Architecture Diagrams

### Data Flow

```
User Input (Variant)
       ↓
Streamlit Dashboard
       ↓
FastAPI Backend
       ↓
Variant Engine
       ↓
┌──────────────┬──────────────┬──────────────┐
│              │              │              │
Enformer Model  API Layer     Cache Layer
│              │              │              │
└──────────────┴──────────────┴──────────────┘
       ↓              ↓              ↓
   Predictions   External APIs   Local Data
       ↓              ↓              ↓
       └──────────────┴──────────────┘
                     ↓
            Integrated Results
                     ↓
          Visualization & Export
```

---

*This document was prepared for PhD application purposes and represents the technical scope and research value of the CardioVar project.*
