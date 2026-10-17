import io

def parse_clinical_vcf(vcf_stream):
    """
    Parses the clinical VCF stream. Keeps the chromosome and base-pair 
    positions intact for accurate distance calculations.
    """
    variant_records = []
    for line in vcf_stream:
        clean_line = line.strip()
        if not clean_line or clean_line.startswith("#"):
            continue
        fields = clean_line.split("\t")
        if len(fields) < 10:
            continue

        chrom = fields[0]              # Extract chromosome (e.g., 'chr1')
        position = int(fields[1])      # Physical position in base-pairs (bp)
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
    """
    Scans the genome using physical base-pair (bp) distances.
    Default max_window_bp is 5 Megabases (a standard clinical block size).
    """
    detected_roh_indices = set()
    
    # Isolate records by chromosome up front to preserve biological integrity and speed up processing
    chrom_groups = {}
    for global_idx, var in enumerate(parsed_variants):
        chrom_groups.setdefault(var["chrom"], []).append((global_idx, var))
        
    for chrom, indexed_vars in chrom_groups.items():
        # Ensure variants are strictly sorted by linear genomic position
        indexed_vars.sort(key=lambda x: x[1]["position"])
        chrom_loci_count = len(indexed_vars)

        # Slide through every variant on the chromosome as a starting point
        for start_idx in range(chrom_loci_count):
            global_start_idx, start_var = indexed_vars[start_idx]
            start_bp = start_var["position"]

            end_idx = start_idx
            het_count = 0

            # Extend the window physically down the chromosome
            while end_idx < chrom_loci_count:
                global_curr_idx, curr_var = indexed_vars[end_idx]
                
                # Stop if window size exceeds maximum physical distance
                if curr_var["position"] - start_bp > max_window_bp:
                    break

                if curr_var["genotype"] == 1:
                    het_count += 1

                # Terminate window if it contains too many heterozygotes
                if het_count > max_het_allowed:
                    break

                end_idx += 1

            # Clinical Rule: If the window spans a solid cluster (e.g., 3+ variants), mark it as ROH
            if (end_idx - start_idx) >= 3:
                for idx in range(start_idx, end_idx):
                    detected_roh_indices.add(indexed_vars[idx][0])

    # Isolate the high-risk targets trapped inside those validated ROH blocks
    isolated_targets = []
    for idx, variant in enumerate(parsed_variants):
        is_inside_roh_block = idx in detected_roh_indices
        is_high_risk_recessive = (variant["genotype"] == 2) and (variant["impact"] in ["HIGH", "PATHOGENIC"])
        if is_inside_roh_block and is_high_risk_recessive:
            isolated_targets.append(variant)

    return detected_roh_indices, isolated_targets

# Mock Patient VCF Data for validation
mock_patient_vcf = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tGENE\tIMPACT\tQUAL\tFILTER\tINFO\tFORMAT\tPATIENT_001
chr1\t1000\trs11\tBRCA2\tMODERATE\t60\tPASS\tDP=30\tGT\t0/1
chr1\t2000\trs12\tGJB2\tHIGH\t60\tPASS\tDP=32\tGT\t1/1
chr1\t3000\trs13\tPAH\tLOW\t60\tPASS\tDP=28\tGT\t1/1
chr1\t4000\trs14\tMEFV\tHIGH\t60\tPASS\tDP=35\tGT\t1/1
chr1\t5000\trs15\tCFTR\tHIGH\t60\tPASS\tDP=31\tGT\t1/1
chr1\t6000\trs16\tALDH2\tMODERATE\t60\tPASS\tDP=29\tGT\t0/1
"""

if __name__ == "__main__":
    print("=========================================================")
    print("     RUNNING OPENROH-FILTER WORKFLOW ENGINE (v2.0)       ")
    print("=========================================================")
    vcf_stream = io.StringIO(mock_patient_vcf)
    clinical_records = parse_clinical_vcf(vcf_stream)
    roh_blocks, clinical_culprits = analyze_patient_genome_by_distance(clinical_records)

    print("---------------------------------------------------------")
    print("             FINAL CLINICAL DIAGNOSTIC REPORT            ")
    print("---------------------------------------------------------")
    if clinical_culprits:
        for culprit in clinical_culprits:
            print(f"🔥 PATHOGENIC TARGET ISOLATED: {culprit['id']}")
            print(f"🧬 Impacted Gene Locus      : {culprit['gene']}")
            print(f"📍 Chromosomal Map Position : Locus Position {culprit['position']} bp")
            print(f"⚠️ Clinical Severity Score  : {culprit['impact']} (Pathogenic)")
            print("---------------------------------------------------------")
    else:
        print("✅ No high-risk pathogenic targets identified.")
