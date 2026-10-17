# OpenROH-Filter: Clinical Screening Engine for Consanguineous Marriage Genetic Profiles

Developed as an academic project for Columbia University, **OpenROH-Filter** is an open-source workflow engine and decision-support tool designed for medical geneticists and clinical researchers. 

The software scans patient genomic data (**Variant Call Format - VCF**) to accurately map physical **Runs of Homozygosity (ROH)** and isolate high-risk autosomal recessive mutations caught within those blocks. This targeted screening pipeline assists in identifying latent genetic risks in consanguineous reproductive counseling.

---

## 🔬 Scientific & Clinical Context

In population genetics and clinical counseling, consanguinity (cousin marriage) increases the coefficient of inbreeding (F), leading to elongated chromosomal segments where both maternal and paternal alleles are identical by descent (IBD). These regions manifest as continuous **Runs of Homozygosity (ROH)**.

While an ROH block itself is a structural genomic feature, it creates a high-risk zone for health complications. If a benign or pathogenic mutation is inherited recessively within an ROH zone, the individual will present as homozygous alternative (`1/1`), completely removing the protective effect of a wild-type allele. OpenROH-Filter automates the scanning of these regions to isolate critical medical culprits (e.g., pathogenic variants in genes like *GJB2*, *MEFV*, or *CFTR*) before clinical presentation.

---

## 🧮 Algorithmic Architecture (How it Works)

The engine utilizes a custom, two-pass sliding window architecture based on physical chromosomal distance rather than index indexing, preventing false-positive flags from varying variant densities.

### Pass 1: Physical Base-Pair Windows
The software slides along chromosomal coordinates (POS) evaluating physical windows up to **5 Megabases (Mb)** (5,000,000 bp). 
* A window is initialized at a variant position and extended down the chromosome.
* The window continuously counts heterozygous loci (`0/1` or `1/0`).
* If the number of heterozygotes exceeds the strict threshold (≤ 1), the window is terminated.
* **Density Filter:** To ensure structural integrity, a region is only classified as a valid ROH block if it contains a dense cluster of **≥ 3 consecutive variants** meeting the hom-alt requirement.

### Pass 2: Pathogenic Isolate Extraction
Once the structural boundaries of the ROH zones are locked, the engine scans the encapsulated loci. It triggers a severe clinical warning *only* when a variant satisfies two simultaneous conditions:
1. It sits entirely inside a validated ROH coordinate block.
2. It is classified as homozygous alternative (`1/1`) **and** carries a `HIGH` or `PATHOGENIC` clinical impact score.

---

## 💻 Repository Structure

```text
openroh-filter/
├── main.py            # Core genomic engine backend (Base-pair calculation script)
├── app.py             # Mobile-responsive web dashboard (Streamlit user interface)
├── requirements.txt   # Server deployment dependencies (Streamlit framework)
└── README.md          # Project documentation, scientific framing, and guide
```

---

## 🔒 Zero-Storage Data Privacy Protocol

Because patient genomic data represents highly sensitive Protected Health Information (PHI), OpenROH-Filter is architected with a strict **zero-storage security policy**:
* **100% Client-Side Processing:** The application processes VCF streams entirely within the local browser memory via `io.BytesIO`.
* No data is written to disk, cached, or transferred to any external cloud database or third-party server.
* This execution structure makes the application safely align with modern digital health privacy philosophies.

---

## ⚙️ Local Deployment & Execution

To test the interactive clinical interface locally on a machine, execute the following commands in your terminal:

```bash
# 1. Install the lightweight framework dependencies
pip install -r requirements.txt

# 2. Launch the responsive web dashboard
streamlit run app.py
```

*Disclaimer: This software is designed as a decision-support tool for academic and clinical research validation. All automated outputs should be cross-verified against official curated databases (like ClinVar) by a board-certified medical geneticist or genetic counselor.*
