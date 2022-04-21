# SARS_CoV19_INHOUSE
Pipeline is used to process SARS-CoV19 WGS data.

## Workflow summary
1. The analysis begins with **quality control of raw fastq files**:
   - Adapter sequences are removed using *cutadapt*[1] (screening for [default illumina adapter sequences](https://support-docs.illumina.com/SHARE/AdapterSeq/Content/SHARE/AdapterSeq/TruSeq/UDIndexes.htm))
   - Reads are trimmed based on base quality (Phred score >= 30) and length (read length >= 50bp) using *fastp*[2]
2. **Reads that passed QC are aligned to the reference** using *bwa*[3].
   - The output from *bwa*[3] is passed to *samtools*[6] to produce binary alignment map (bam) file.
3. **Primer sequences are removed from bam file** using *ivar*[4].
   - Primer sequences should be provided in Browser Extensible Data (bed) format - generated from fasta file using *bwa*[3] & *bedtools*[5]. 
4. To improve alignment quality, **local realignment is performed** on primer-free bam file using *abra*[9].
   - bed file with realignment targets is created from primer-free bam file using *bedtools*[5].
5. **Alignment QC metrics are extracted** from raw and primer-free bam files using *samtools*[6].
6. **Variant-calling is performed** on post-realignment bam file using *freebayes*[7].
7. **Raw variants are filtered** using *vcflib/vcffilter*[10] based on quality ([QUAL](https://samtools.github.io/hts-specs/VCFv4.1.pdf) > 30) and sequencing depth ([DP](https://samtools.github.io/hts-specs/VCFv4.1.pdf) > 15).
8. **Filtered variants are annotated** using *snpEff*[8], using [genbank reference](https://www.ncbi.nlm.nih.gov/nuccore/MN908947).
9. **Consensus sequence is generated** (in fasta format) from post-realignment bam file using *ivar*[4].
    - invalid base (N) is called if coverage is less that 15.
10. **Sample id is added to fasta header and invalid bases are replaced with N** using bash scripts.
11. **Annotated vcf file is converted to csv format and coverage depth plot is generated**
from sequencing depth data using inhouse-developed python scripts.
12. **Temporary files are deleted** after each sample is processed.
13. **Lineage assignment is performed** based on consensus sequence using *Pangolin*[11].
14. Inhouse-developed python & bash scripts are used to control the flow of analysis for multiple samples, generate summary report and visualize the results.
    

## Tools & References
1. cutadapt 2.31 - https://doi.org/10.14806/ej.17.1.200
2. fastp 0.20.1 - https://doi.org/10.1093/bioinformatics/bty560
3. bwa 0.7.17-r1198-dirty - https://arxiv.org/abs/1303.3997
4. ivar 1.3.1 - https://doi.org/10.1186/s13059-018-1618-7
5. bedtools v.2.30.00 - https://doi.org/10.1093/bioinformatics/btq033
6. samtools 1.12 - https://doi.org/10.1093/bioinformatics/btp352
7. freebayes v0.9.21 - https://arxiv.org/abs/1207.3907
8. snpEff 5.0e - https://pcingola.github.io/SnpEff/adds/SnpEff_paper.pdf
9. abra 0.97 - https://doi.org/10.1093/bioinformatics/btu376
10. vcflib 1.0.2 - https://doi.org/10.1101/2021.05.21.445151 
11. Pangolin - https://doi.org/10.1038/s41564-020-0770-5
