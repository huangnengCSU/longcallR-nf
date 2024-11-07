## Quick Run
Before running this pipeline, make sure that [Nextflow](https://www.nextflow.io/docs/latest/install.html) is installed. This workflow contains several steps:

- minimap2 alignment
- calling variants using `longcallR-nn` (Whole chromosome)
- calling variants using `longcallR-nn + longcallR`, phasing (Whole chromosome)
- calling variants using `longcallR-nn + longcallR`, phasing (Exon)
- calling variants only using `longcallR`, phasing (Whole chromosome)
- transcript analysis with IsoQuant
### Clone the repo
```bash
git clone https://github.com/huangnengCSU/longcallR-nf.git
cd longcallR-nf
```
### run from reads file
```bash

nextflow run main.nf \
-with-trace trace.txt \
--sample_name sample1 \   # input sample name
--reads /path/to/flnc.fastq \   # support multiple reads like "reads1.fastq,reads2.fastq"
--ref /path/to/ref.fa \
--platform pb \ # choice: pb/ont
--datatype masseq \ # for pb: masseq/isoseq, for ont: cDNA/dRNA
--outdir /path/to/outdir \
--annotation /path/to/gencode.v44.annotation.gtf \  # only support for gtf file, isoquant has issue with gff3 file
--threads 50 \
--threads_per_job 2 \   # setting of threads pools
--no_cuda true  # true for cpu mode and false for gpu mode

```

### run from bam file
```bash

nextflow run main.nf \
-with-trace trace.txt \
--sample_name sample1 \   # input sample name
--bam /path/to/aln.sort.bam \
--ref /path/to/ref.fa \
--platform pb \ # choice: pb/ont
--datatype masseq \ # for pb: masseq/isoseq, for ont: cDNA/dRNA
--outdir /path/to/outdir \
--annotation /path/to/gencode.v44.annotation.gtf \  # only support for gtf file, isoquant has issue with gff3 file
--threads 50 \
--threads_per_job 2 \   # setting of threads pools
--no_cuda true  # true for cpu mode and false for gpu mode

```