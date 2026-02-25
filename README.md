# oist/luscombeu_rrnascan

[![Open in GitHub Codespaces](https://img.shields.io/badge/Open_In_GitHub_Codespaces-black?labelColor=grey&logo=github)](https://github.com/codespaces/new/oist/luscombeu_rrnascan)
[![GitHub Actions CI Status](https://github.com/oist/luscombeu_rrnascan/actions/workflows/nf-test.yml/badge.svg)](https://github.com/oist/luscombeu_rrnascan/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/oist/luscombeu_rrnascan/actions/workflows/linting.yml/badge.svg)](https://github.com/oist/luscombeu_rrnascan/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.04.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.5.1-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.5.1)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/oist/luscombeu_rrnascan)

## Introduction

**oist/luscombeu_rrnascan** is a comprehensive Nextflow pipeline for detecting, extracting, and analyzing ribosomal RNA (rRNA) and internal transcribed spacer (ITS) sequences from genomic assemblies or sequencing reads. The pipeline supports multiple input types: Illumina paired-end reads (assembled via NOVOPlasty), long-read sequencing data (ONT/PacBio), or pre-assembled genome sequences. It searches for SSU (18S), 5.8S, and LSU (28S) rRNA sequences using Infernal's covariance models and automatically extracts the ITS region spanning ITS1, ITS2 and 5.8S.

## Pipeline Overview

1. **Assembly** (if reads provided):
   - Illumina reads → [`NOVOPlasty`](https://github.com/eversinc33/NOVOPlasty) (de novo circular assembly)
   - ONT/PacBio reads → [`NOVOloci`](https://github.com/eversinc33/NOVOloci) (long-read assembly NOT READY)
   - Pre-assembled genomes → Direct input

2. **rRNA Detection & Extraction**:
   - rRNA discovery ([`Infernal`](http://infernal.janelia.org/) cmsearch with covariance models)
   - Sequence extraction ([`EMBOSS`](http://emboss.sourceforge.net/) seqret)
   - ITS region identification (between rRNA boundaries)

3. **Quality Control**:
   - Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
   - Results summary ([`MultiQC`](http://multiqc.info/))

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

### Prepare Input Samplesheet

Create a CSV file (`samplesheet.csv`) with your input data:

```csv
sample_ID,read_type,read_1,read_2,seed
SAMPLE1,illumina,reads_R1.fastq.gz,reads_R2.fastq.gz,seed.fasta
SAMPLE2,genome,assembly.fasta,,
SAMPLE3,ONT,nanopore_reads.fastq.gz,,
```

**Columns:**
- `sample_ID`: Unique sample identifier
- `read_type`: One of `illumina`, `ONT`, `Pacbio`, or `genome`
- `read_1`: Path to forward reads (FASTQ/FASTA) or assembly (FASTA)
- `read_2`: Path to reverse reads (FASTQ, optional for Illumina)
- `seed`: Optional seed sequence for NOVOPlasty (FASTA)

### Run Pipeline

```bash
nextflow run oist/luscombeu_rrnascan \
   -profile singularity \
   --input samplesheet.csv \
   --outdir results
```

### Key Parameters

```bash
# Assembly parameters
--genome_range '6000-12000'    # Expected genome size range for NOVOPlasty
--kmer 33                       # K-mer size for assembly (default: 33)
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

## Output

Results are saved in the `--outdir` directory with the following structure:

```
results/
├── assembly/          # NOVOPlasty assemblies (if reads provided)
├── cmsearch/          # rRNA search results (tab-separated)
├── extract/      # Extracted rRNA and ITS sequences (FASTA)
│   ├── *_SSU.fasta    # 18S rRNA sequences
│   ├── *_ITS1.fasta   # ITS1 regions
│   ├── *_5.8S.fasta   # 5.8S rRNA sequences
│   ├── *_ITS2.fasta   # ITS2 regions
│   ├── *_LSU.fasta    # 28S rRNA sequences
│   └── *_rrna.tsv     # Detailed rRNA coordinates
├── fastqc/            # FastQC reports
├── multiqc/           # MultiQC summary
└── pipeline_info/     # Pipeline execution details
```

## Credits

oist/luscombeu_rrnascan was originally written by Johannes Nicolaus Wibisana.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).
