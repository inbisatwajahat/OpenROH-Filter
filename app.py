import streamlit as st
import io
from datetime import datetime

# 1. Page Configuration
st.set_page_config(
    page_title="OpenROH-Filter Clinical Engine",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="collapsed"
)

# 2. Main Title and Clinical Context
st.title("🧬 OpenROH-Filter Workflow Engine")
st.caption("Columbia University Project • Consanguineous Marriage Genetic Screening Tool")

st.markdown("""
This tool scans patient genetic files (**VCF format**) to identify continuous **Runs of Homozygosity (ROH)**. 
It isolates high-risk autosomal recessive mutations caught inside these regions to assist clinical diagnostics.
""")

# 3. Privacy Guarantee Banner
st.info("🔒 **Data Privacy Notice:** All file processing happens entirely inside your browser/local machine. No patient genetic data is ever uploaded or saved to an external server.")

# 4. Core Genomic Engine Logic
def parse_clinical_vcf(vcf_stream):
    variant_records = []
    lines = vcf_stream.getvalue().decode("utf-8").splitlines()
    
    for line in lines:
        clean_line = line.strip()
        if not clean_line or clean_line.startswith("#"):
            continue
        fields = clean_line.split("\t")
        if len(fields) < 10:
            continue
            
        chrom = fields[0]            
        position = int(fields[1])    
        variant_id = fields[2]        
        gene_name = fields[3]         
        clinical_impact = fields[4].strip().upper() 
        genotype_data = fields[9]     
        
        if "0/1" in genotype_data or "1/0" in genotype_data:
            numerical_genotype = 1    
        elif "1/1" in genotype_data:
            numerical_genotype = 2    
        else:
            numerical_genotype = 0    
            
        variant_records.append({
            "chrom": chrom,
            "position": position,
            "id": variant_id,
            "gene": gene_name,
            "impact": clinical_impact,
            "genotype": numerical_genotype
        })
    return variant_records

def analyze_patient_genome_by_distance(parsed_variants, max_window_bp=5000000, max_het_allowed=1):
    detected_roh_indices = set()
    total_loci = len(parsed_variants)
    
    for start_idx in range(total_loci):
        current_chrom = parsed_variants[start_idx]["chrom"]
        start_bp = parsed_variants[start_idx]["position"]
        end_idx = start_idx
        het_count = 0
        
        while end_idx < total_loci and parsed_variants[end_idx]["chrom"] == current_chrom:
            if parsed_variants[end_idx]["position"] - start_bp > max_window_bp:
                break
            if parsed_variants[end_idx]["genotype"] == 1:
                het_count += 1
            if het_count > max_het_allowed:
                break
            end_idx += 1
            
        if (end_idx - start_idx) >= 3:
            for idx in range(start_idx, end_idx):
                detected_roh_indices.add(idx)
                
    isolated_targets = []
    for idx, variant in enumerate(parsed_variants):
        is_inside_roh_block = idx in detected_roh_indices
        is_high_risk_recessive = (variant["genotype"] == 2) and (variant["impact"] in ["HIGH", "PATHOGENIC"])
        if is_inside_roh_block and is_high_risk_recessive:
            isolated_targets.append(variant)
            
    return detected_roh_indices, isolated_targets

# 5. Drag-and-Drop File Upload UI Component
uploaded_file = st.file_uploader("📤 Drag and drop your patient VCF file here", type=["vcf", "txt"])

# 6. Interactive Sample Data Button
mock_vcf_data = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tGENE\tIMPACT\tQUAL\tFILTER\tINFO\tFORMAT\tPATIENT_001
chr1\t1000\trs11\tBRCA2\tMODERATE\t60\tPASS\tDP=30\tGT\t0/1
chr1\t2000\trs12\tGJB2\tHIGH\t60\tPASS\tDP=32\tGT\t1/1
chr1\t3000\trs13\tPAH\tLOW\t60\tPASS\tDP=28\tGT\t1/1
chr1\t4000\trs14\tMEFV\tHIGH\t60\tPASS\tDP=35\tGT\t1/1
chr1\t5000\trs15\tCFTR\tHIGH\t60\tPASS\tDP=31\tGT\t1/1
chr1\t6000\trs16\tALDH2\tMODERATE\t60\tPASS\tDP=29\tGT\t0/1"""

if st.button("💡 Click here to load sample data and see how it works"):
    uploaded_file = io.BytesIO(mock_vcf_data.encode("utf-8"))

# 7. Processing and Displaying Results
if uploaded_file is not None:
    st.markdown("---")
    st.subheader("📋 Final Clinical Diagnostic Report")
    
    # Process the data
    with st.spinner("Analyzing genomic structure..."):
        clinical_records = parse_clinical_vcf(uploaded_file)
        roh_blocks, clinical_culprits = analyze_patient_genome_by_distance(clinical_records)
    
    # Generate the text string for the clinical download file
    report_text = f"=========================================================\n"
    report_text += f"         OPENROH-FILTER CLINICAL DIAGNOSTIC REPORT       \n"
    report_text += f"=========================================================\n"
    report_text += f"Generated On: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
    report_text += f"Data Privacy Status: Locally Processed (Secure)\n"
    report_text += f"---------------------------------------------------------\n\n"
    
    if clinical_culprits:
        st.error(f"⚠️ {len(clinical_culprits)} High-Risk Pathogenic Target(s) Isolated!")
        report_text += f"STATUS: {len(clinical_culprits)} PATHOGENIC TARGET(S) ISOLATED\n\n"
        
        for culprit in clinical_culprits:
            # Display visual cards on the webpage
            with st.container():
                st.markdown(f"### 🔥 PATHOGENIC TARGET ISOLATED: **{culprit['id']}**")
                col1, col2 = st.columns(2)
                with col1:
                    st.write(f"🧬 **Impacted Gene Locus:** {culprit['gene']}")
                    st.write(f"📍 **Chromosomal Position:** Locus {culprit['position']} bp")
                with col2:
                    st.write(f"⚠️ **Clinical Severity:** {culprit['impact']} (Pathogenic)")
                st.markdown("---")
            
            # Append item data to the downloadable report text
            report_text += f"🔥 PATHOGENIC TARGET: {culprit['id']}\n"
            report_text += f"🧬 Gene Locus             : {culprit['gene']}\n"
            report_text += f"📍 Map Position          : {culprit['position']} bp\n"
            report_text += f"⚠️ Clinical Severity Score: {culprit['impact']} (Pathogenic)\n"
            report_text += f"---------------------------------------------------------\n"
    else:
        st.success("✅ No high-risk pathogenic targets identified inside homozygosity blocks.")
        report_text += "STATUS: CLEAR. No high-risk pathogenic targets identified inside homozygosity blocks.\n"
        
    report_text += "\n\nDisclaimer: This report is generated automatically for research use. Findings must be cross-verified by a certified medical geneticist or genetic counselor before clinical action."

    # 8. Download Button Interface Component
    st.download_button(
        label="📥 Download Clinical Report (.txt)",
        data=report_text,
        file_name=f"OpenROH_Report_{datetime.now().strftime('%Y%md_%H%M%S')}.txt",
        mime="text/plain"
    )
        
    # Patient Educational Disclaimer
    st.caption("ℹ️ *Disclaimer: This report is generated automatically for research use. Findings must be cross-verified by a certified medical geneticist or genetic counselor before clinical action.*")
