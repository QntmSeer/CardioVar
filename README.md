# CardioVar: CVD Variant Impact Explorer

![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)
![Python](https://img.shields.io/badge/python-3.11-blue.svg)

**CardioVar** is a precision medicine tool designed to predict the functional impact of genetic variants in cardiovascular disease genes. By leveraging deep learning models, it estimates how a specific variant alters RNA-seq expression profiles, potentially disrupting gene regulation.

## 🌟 Features

- **📉 Predict Impact**: Visualize Δ RNA-seq effects of single nucleotide variants
- **🧬 Genomic Context**: View gene structure (exons) and evolutionary conservation
- **🏥 Clinical Relevance**: Cross-reference with ClinVar and GWAS Catalog
- **🚀 Batch Processing**: Analyze multiple variants simultaneously
- **🧪 Tissue-Specific Analysis**: See how variants affect different cardiovascular tissues
- **📊 Percentile Ranking**: Compare variant impact against background distribution

## 🏗️ Architecture

CardioVar uses a **Client-Server Architecture**:

- **Backend (FastAPI)**: REST API serving variant predictions and annotations
- **Frontend (Streamlit)**: Interactive dashboard for visualization and analysis
- **Data Layer**: JSON-based mock data for genes and variants

```
cardiovar/
├── api.py                  # FastAPI backend
├── variant_engine.py       # Core prediction logic
├── dashboard.py            # Streamlit frontend
├── data/
│   ├── gene_annotations.json
│   └── related_variants.json
├── tests/
│   └── test_qc_backend.py
└── requirements.txt
```

## 🚀 Quick Start

### Prerequisites

- Python 3.11+
- pip

### Installation

1. **Clone the repository**
```bash
git clone https://github.com/YOUR_USERNAME/cardiovar.git
cd cardiovar
```

2. **Install dependencies**
```bash
pip install -r requirements.txt
```

3. **Set up environment variables**
Create a `.env` file:
```bash
ALPHAGENOME_API_KEY=your_api_key_here
```

### Running the Application

You need to run **two** processes:

**Terminal 1: Start the Backend**
```bash
python -m uvicorn api:app --host 0.0.0.0 --port 8000 --reload
```

**Terminal 2: Start the Dashboard**
```bash
streamlit run dashboard.py
```

Then open your browser to **http://localhost:8501**

## 📊 Usage

### Single Variant Analysis

1. Enter variant details in the sidebar:
   - Chromosome (e.g., `chr22`)
   - Position (e.g., `36191400`)
   - Reference allele (e.g., `A`)
   - Alternate allele (e.g., `C`)

2. Click **Run Analysis**

3. Explore the results across 4 tabs:
   - **Variant Explorer**: Main Δ RNA-seq plot, tissue-specific impact, percentile ranking
   - **Gene Annotations**: Expression profiles, protein domains, pathways
   - **Related Data**: ClinVar and GWAS associations
   - **Batch Analysis**: Upload CSV for multiple variants

### Batch Analysis

Upload a CSV file with columns: `chrom`, `pos`, `ref`, `alt`

Example:
```csv
chrom,pos,ref,alt
chr22,36191400,A,C
chr1,12345678,G,T
```

## 🧪 Testing

Run the automated test suite:

```bash
pytest tests/test_qc_backend.py -v
```

## 📦 API Endpoints

### `POST /variant-impact`
Compute impact for a single variant.

**Request:**
```json
{
  "chrom": "chr22",
  "pos": 36191400,
  "ref": "A",
  "alt": "C"
}
```

**Response:**
```json
{
  "variant_id": "chr22:36191400:A:C",
  "metrics": {
    "max_delta": 3.512,
    "gnomad_freq": 0.00005,
    "gene_symbol": "MYH9",
    "percentile": 95.2
  },
  "curve": { "x": [...], "y": [...] },
  "tissue_effects": [...],
  "background_distribution": {...}
}
```

### `GET /gene-annotations?gene=MYH9`
Get detailed gene annotations including expression and protein domains.

### `POST /batch-impact`
Process multiple variants in batch.

## 🎨 Customization

- **Plot Colors**: Adjust signal and highlight colors in the sidebar
- **Theme**: Modify `.streamlit/config.toml` for custom themes

## 📝 Data Sources

Currently uses **mock data** for demonstration. To integrate real data:

1. Replace `variant_engine.py` logic with actual deep learning model
2. Connect to real APIs:
   - ClinVar: https://www.ncbi.nlm.nih.gov/clinvar/
   - gnomAD: https://gnomad.broadinstitute.org/
   - GTEx: https://gtexportal.org/

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## 📄 License

This project is licensed under the Apache License 2.0 - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- Mock gene data inspired by GTEx and Ensembl
- Visualization design influenced by UCSC Genome Browser and IGV

## 📧 Contact

For questions or feedback, please open an issue on GitHub.

---

**Note**: This is a research prototype using synthetic data. Not for clinical use.
