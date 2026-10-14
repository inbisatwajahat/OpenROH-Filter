import io

def parse_clinical_vcf(vcf_stream):
    variant_records = []
    for line in vcf_stream:
        clean_line = line.strip()
        if not clean_line or clean_line.startswith("#"):
            continue
        fields = clean_line.split("\t")
        position = int(fields[1])      
        variant_id = fields[2]        
        gene_name = fields[3]         
        clinical_impact = fields[4].strip() 
        genotype_data = fields[9]     
        if "0/1" in genotype_data or "1/0" in genotype_data:
            numerical_genotype = 1    
        elif "1/1" in genotype_data:
            numerical_genotype = 2    
        else:
            numerical_genotype = 0    
        variant_records.append({
            "position": position,
            "id": variant_id,
            "gene": gene_name,
            "impact": clinical_impact,
            "genotype": numerical_genotype
        })
    return variant_records

def analyze_patient_genome(parsed_variants, window_size=3, max_het_allowed=1):
    genotype_array = [v["genotype"] for v in parsed_variants]
    total_loci = len(genotype_array)
    detected_roh_indices = set()
    for start in range(0, total_loci - window_size + 1):
        end = start + window_size
        window_slice = genotype_array[start:end]
        het_count = window_slice.count(1)
        if het_count <= max_het_allowed:
            for idx in range(start, end):
                detected_roh_indices.add(idx)
    isolated_targets = []
    for idx, variant in enumerate(parsed_variants):
        is_inside_roh_block = idx in detected_roh_indices
        is_high_risk_recessive = (variant["genotype"] == 2) and (variant["impact"] == "HIGH")
        if is_inside_roh_block and is_high_risk_recessive:
            isolated_targets.append(variant)
    return detected_roh_indices, isolated_targets

mock_patient_vcf = """##fileformat=VCFv4.2
##source=AyubMedicalCollegeGenomicsPlatform
#CHROM	POS	ID	GENE	IMPACT	QUAL	FILTER	INFO	FORMAT	PATIENT_001
chr1	1000	rs11	BRCA2	MODERATE	60	PASS	DP=30	GT	0/1
chr1	2000	rs12	GJB2	HIGH    	60	PASS	DP=32	GT	1/1
chr1	3000	rs13	PAH  	LOW     	60	PASS	DP=28	GT	1/1
chr1	4000	rs14	MEFV 	HIGH    	60	PASS	DP=35	GT	1/1
chr1	5000	rs15	CFTR 	HIGH    	60	PASS	DP=31	GT	1/1
chr1	6000	rs16	ALDH2	MODERATE	60	PASS	DP=29	GT	0/1
"""

print("=========================================================")
print("         RUNNING OPENROH-FILTER WORKFLOW ENGINE          ")
print("=========================================================")
vcf_stream = io.StringIO(mock_patient_vcf)
clinical_records = parse_clinical_vcf(vcf_stream)
roh_blocks, clinical_culprits = analyze_patient_genome(clinical_records, window_size=3, max_het_allowed=1)

print("---------------------------------------------------------")
print("             FINAL CLINICAL DIAGNOSTIC REPORT            ")
print("---------------------------------------------------------")
if clinical_culprits:
    for culprit in clinical_culprits:
        print(f"🔥 PATHOGENIC TARGET ISOLATED: {culprit['id']}")
        print(f"🧬 Impacted Gene Locus      : {culprit['gene']}")
        print(f"📍 Chromosomal Map Position : Locus Position {culprit['position']}")
        print(f"⚠️ Clinical Severity Score  : {culprit['impact']} (Pathogenic)")
        print("---------------------------------------------------------")
