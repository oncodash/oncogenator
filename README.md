# Oncodash Annotation Tool

This tool is used for post-annotation of raw VCF files that have been annotated with Annovar. It performs additional annotation based on copy number information, estimates the homogeneity of somatic alterations, and classifies somatic alterations by their functional consequences.

## Description

The tool processes somatic variants and copy number alterations (CNAs) to provide a comprehensive annotation of genomic alterations. It integrates data from various sources, including OncoKB and Cancer Genome Interpreter (CGI), to enhance the annotation with clinical relevance and treatment information.

### Features

- **Somatic Variant Annotation**: Annotates somatic variants with functional consequences and clinical relevance.
- **Copy Number Alteration Annotation**: Annotates CNAs with functional consequences and clinical relevance.
- **Homogeneity Estimation**: Estimates the homogeneity of somatic alterations.
- **Classification**: Classifies somatic alterations by their functional consequences.

### Requirements
- Python 3.x
- OncoKB API token
- Cancer Genome Interpreter API access

### Installation
Clone the repository.
Install the required Python packages:
>pip install -r requirements.txt

Set up the virtual environment and activate it:
> python -m venv venv
source venv/bin/activate

Configure the paths and API tokens in the scripts as needed.

## Usage

### Command Line Options

```text
Usage: main.py [OPTIONS]

Options:
  --annotator <local|external|both>  Choose which annotator to run (required).
  --output <str>                     Path to output file (required).
  --somatic_variants <str>           Path to somatic variants file.
  --copy_number_alterations <str>    Path to copy number alterations file.
  --ascatestimates <str>             Path to ASCAT estimates file (required).
  --cn_annotations <str>             Path to filtered and annotated CNAs.
  --tumortype <str>                  Tumor type identifier (default: HGSOC).
  --refgen <str>                     Reference genome version (default: GRCh38).
  --cores <int>                      Number of cores to use for processing (default: 0).
  --post_annotate                    Perform post-annotation.
  --ploidy_threshold <float>         Threshold for filtering by ploidy (default: 2.5).
  --pid <str>                        Patient ID.
  --homogeneity_threshold <float>    Homogeneity threshold for variant filtering (default: 0.05).
  --rf_score_threshold <float>       Random Forest score threshold for variant filtering (default: 0.95).
  --ada_score_threshold <float>      AdaBoost score threshold for variant filtering (default: 0.95).
  --oncokbcna                        Query OncoKB for copy number alterations.
  --oncokbsnv                        Query OncoKB for somatic mutations.
  --cgiquery                         Query Cancer Genome Interpreter.
  --cgijobid <str>                   Download results from CGI by jobid.
```

### Examples
Local Annotator
1. Annotate somatic variants and CNAs using the local annotator:
> python main.py --annotator local --output path/to/output --somatic_variants path/to/snvs.tsv --ascatestimates path/to/ascat.tsv

2. Annotate CNAs using the local annotator:
> python main.py --annotator local --output path/to/output --copy_number_alterations path/to/cnas.tsv --ascatestimates path/to/ascat.tsv

External Annotator (execution order is important as the CGI annotations are supplementing the OncoKB annotations in default)
3. Annotate CNAs using OncoKB:
> python main.py --annotator external --output path/to/output --oncokbcna --copy_number_alterations path/to/locally_annotated_cnas.tsv

4. Annotate somatic variants using OncoKB:
>python main.py --annotator external --output path/to/output --oncokbsnv --somatic_variants path/to/locally_annotated_snvs.tsv

5. Annotate CNAs using Cancer Genome Interpreter:

>python external_annotator.py --cgiquery --copy_number_alterations path/to/oncokb_annotated_cnas.tsv

6. Annotate somatic variants using Cancer Genome Interpreter:
>python main.py --annotator external --output path/to/output --cgiquery --somatic_variants path/to/oncokb_annotated_snvs.tsv

SLURM Scripts: edit the scripts to set the correct paths and SLURM sbatch parameters.

Submit a batch job to SLURM cluster to annotate on multiple computing nodes:
>./slurm_scripts/annotate_cnas.sh path/to/sample_list.txt

>./slurm_scripts/snv_annotation.sbatch path/to/sample_list.txt
