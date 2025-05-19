import argparse

import pandas as pd
import dask.dataframe as dd
import os

from copy_number_annotator import CopyNumberAnnotator
from utils import df_apply
from somatic_variant_annotator import SomaticVariantAnnotator

'''
    Usage: local_annotator.py [OPTIONS]
    
    Options:
      --output <str>                   Path to output file (required).
      --somatic_variants <str>         Path to somatic variants file.
      --copy_number_alterations <str>  Path to copy number alterations file.
      --ascatestimates <str>           Path to ASCAT estimates file (required).
      --cn_annotations <str>           Path to filtered and annotated CNAs.
      --tumortype <str>                Tumor type identifier (default: HGSOC).
      --refgen <str>                   Reference genome version (default: GRCh38).
      --cores <int>                    Number of cores to use for processing (default: 0).
      --post_annotate                  Perform post-annotation.
      --ploidy_threshold <float>       Threshold for filtering by ploidy (default: 2.5).
      --pid <str>                      Patient ID.
      --homogeneity_threshold <float>  Homogeneity threshold for variant filtering (default: 0.05).
      --rf_score_threshold <float>     Random Forest score threshold for variant filtering (default: 0.95).
      --ada_score_threshold <float>    AdaBoost score threshold for variant filtering (default: 0.95).
    
    Examples:
      python local_annotator.py --output path/to/output --somatic_variants path/to/snvs.tsv --ascatestimates path/to/ascat.tsv
      python local_annotator.py --output path/to/output --copy_number_alterations path/to/cnas.tsv --ascatestimates path/to/ascat.tsv
      python local_annotator.py --output path/to/output --somatic_variants path/to/snvs.tsv --cn_annotations path/to/cnas --pid patient_id --ascatestimates path/to/ascat.tsv
      python local_annotator.py --output path/to/output --somatic_variants path/to/snvs.tsv --post_annotate --ascatestimates path/to/ascat.tsv
'''
def main(**kwargs):
    """
    Import genomic variants produced by sequencing analysis pipelines to populate corresponding models into Django database.

    Usage:
    -------
    `python manage.py import_genomic_variants --somatic_variants <filepath> --copy_number_alterations <filepath> --filter=<field> --equal=<value>`
    """
    help = "Import genomic alterations into the database."

    if kwargs.get("copy_number_alterations"):
        process_copy_number_alterations(kwargs)

    if kwargs.get("somatic_variants") and kwargs.get("cn_annotations") and kwargs.get("pid"):
        process_somatic_variants(kwargs)

    if kwargs.get("somatic_variants") and kwargs.get("post_annotate"):
        post_annotate_somatic_variants(kwargs)

def process_copy_number_alterations(kwargs):
    output = kwargs["output"]
    cores = int(kwargs.get("cores", 1))
    cna_data = pd.read_csv(kwargs["copy_number_alterations"], sep="\t")

    refgenome = kwargs.get("refgen", "GRCh38")
    ploidy_coeff = float(kwargs.get("ploidy_threshold", 2.5))
    tumortype = kwargs.get("tumortype", "HGSOC")

    ascats = pd.read_csv(kwargs["ascatestimates"], sep="\t", encoding='utf-8')
    annotator = CopyNumberAnnotator(refgenome=refgenome, tumortype=tumortype, ascats=ascats, ploidy_coeff=ploidy_coeff)

    cnas_filtered = df_apply(cna_data, annotator.filter_cnas_by_ploidy).dropna()
    cnadf = pd.DataFrame(dict(zip(cnas_filtered.index, cnas_filtered.values))).T
    cnadf.to_csv(output, sep='\t')
    print(cnadf)

def process_somatic_variants(kwargs):
    if kwargs['pid']:
        snv_file = f"{kwargs['somatic_variants']}/{kwargs['pid']}.csv"
    else:
        snv_file = kwargs['somatic_variants']

    output = kwargs["output"]
    cores = int(kwargs.get("cores", 0))
    refgenome = kwargs.get("refgen", "GRCh38")
    tumortype = kwargs.get("tumortype", "HGSOC")
    homogeneity_threshold = float(kwargs.get("homogeneity_threshold", 0.05))
    rf_score_threshold = float(kwargs.get("rf_score_threshold", 0.95))
    ada_score_threshold = float(kwargs.get("ada_score_threshold", 0.95))

    ascats = pd.read_csv(kwargs["ascatestimates"], sep="\t", encoding='utf-8')
    snvs = pd.read_csv(snv_file, sep="\t", encoding='utf-8')
    if kwargs['pid']:
        # Annotate by patient
        cnas = dd.read_csv(f"{kwargs['cn_annotations']}/{kwargs['pid']}*.csv", sep="\t").compute().dropna(subset=['Gene', 'nMinor', 'nMajor', 'LOHstatus'])
    else:
        # Annotate all from single file
        cnas = dd.read_csv(f"{kwargs['cn_annotations']}", sep="\t").compute().dropna(subset=['Gene', 'nMinor', 'nMajor', 'LOHstatus'])

    samples = cnas['sample'].drop_duplicates()

    # Instantiate the SomaticVariantAnnotator
    annotator = SomaticVariantAnnotator(refgenome=refgenome, tumortype=tumortype, cnas=cnas, ascats=ascats, samples=samples, homogeneity_threshold=homogeneity_threshold, rf_score_threshold=rf_score_threshold, ada_score_threshold=ada_score_threshold)
    snvs_filtered = df_apply(snvs, annotator.filter_and_classify_snvs, None, None, True, True, cores=cores)

    joined = [item for sublist in snvs_filtered if isinstance(sublist, list) for item in sublist]
    snvdf = pd.DataFrame(joined)
    snvdf.to_csv(output, sep='\t')
    print(snvdf)


def post_annotate_somatic_variants(kwargs):
    output = kwargs["output"]
    cores = int(kwargs.get("cores", 0))
    refgenome = kwargs.get("refgen", "GRCh38")
    tumortype = kwargs.get("tumortype", "HGSOC")

    ascats = pd.read_csv(kwargs["ascatestimates"], sep="\t", encoding='utf-8')
    snvs = pd.read_csv(kwargs["somatic_variants"], sep="\t", encoding='utf-8')
    cnas = pd.read_csv(kwargs["cn_annotations"], sep="\t", encoding='utf-8')[['sample', 'Gene', 'nMinor', 'nMajor', 'LOHstatus']].dropna(subset=['Gene', 'nMinor', 'nMajor', 'LOHstatus'])

    cna_grps = cnas.groupby('sample')
    snv_grps = snvs.groupby('sample_id')

    for gname, cna_grp in cna_grps:
        process_sample(gname, cna_grp, snv_grps, output, cores, ascats)

def process_sample(gname, cna_grp, snv_grps, output, cores, ascats):
    output_file = f"{output}/{gname}_snvs.csv"
    if os.path.exists(output_file):
        print(f"Output exists for {gname}! Skipping..")
        return

    try:
        annotator = SomaticVariantAnnotator(cnas=cna_grp, ascats=ascats)
        snv_grp = snv_grps.get_group(gname)
        snvs_filtered = df_apply(snv_grp, annotator.post_filter_and_classify_snvs_by_sample, None, None, True, True, cores=cores)
        snvdf = pd.DataFrame(snvs_filtered)
        snvdf.to_csv(output_file, sep='\t', header=False)
    except Exception as e:
        print(f"Error processing sample {gname}: {e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Annotate genomic alterations.")
    parser.add_argument("--output", type=str, required=True, help="Path to output files")
    parser.add_argument("--somatic_variants", type=str, help="Path to somatic variants files")
    parser.add_argument("--copy_number_alterations", type=str, help="Path to copy number alterations files")
    parser.add_argument("--ascatestimates", type=str, required=True, help="Path to ASCAT estimates file")
    parser.add_argument("--cn_annotations", type=str, help="Path to filtered and annotated CNAs")
    parser.add_argument("--tumortype", type=str, default="HGSOC", help="Tumor type identifier (default: HGSOC)")
    parser.add_argument("--refgen", type=str, default="GRCh38", help="Reference genome version (default: GRCh38)")
    parser.add_argument("--cores", type=int, default=0, help="Number of cores to use for processing")
    parser.add_argument("--post_annotate", action='store_true', help="Perform post-annotation")
    parser.add_argument("--ploidy_threshold", type=float, default=2.5, help="Threshold for filtering by ploidy")
    parser.add_argument("--pid", type=str, help="Patient ID")
    parser.add_argument("--homogeneity_threshold", type=float, default=0.05,
                        help="Homogeneity threshold for variant filtering")
    parser.add_argument("--rf_score_threshold", type=float, default=0.95,
                        help="Random Forest score threshold for variant filtering")
    parser.add_argument("--ada_score_threshold", type=float, default=0.95,
                        help="AdaBoost score threshold for variant filtering")

    args = parser.parse_args()

    main(**vars(args))