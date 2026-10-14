# OpenROH-Filter: Algorithmic Mapping of Consanguineous Genomic Blind Spots

An open-source Python framework designed to map **Runs of Homozygosity (ROH)** blocks to uncover pathogenic homozygous recessive variants within consanguineous patient cohorts (cousin-marriage populations).

## 🧬 Scientific Background & Clinical Workflow
Standard clinical genomic pipelines are statistically optimized for outbred populations. In cohorts with high rates of consanguinity, such as populations managed by regional centers like **Ayub Medical College in Abbottabad, Pakistan**, these pipelines hit critical blind spots. 

When related parents pass down extensive identical chromosomal segments, the patient's genome is flooded with massive, non-pathogenic homozygous tracts, masking rare, severe recessive mutations under the background noise.

**OpenROH-Filter** solves this by establishing a custom sliding-window matrix. It tracks the depletion boundaries of heterozygous loci, strips away safe parental background noise, and aggressively isolates high-impact recessive variants (like the *GJB2* gene) trapped within those windows.

