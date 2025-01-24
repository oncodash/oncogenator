import argparse

import pandas as pd
import dask.dataframe as dd
import somatic_variant_annotator as sva
import copy_number_annotator as cna
from utils import df_apply
import os

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
    cnas_filtered = df_apply(cna_data, cna.filter_cnas_by_ploidy, None, None, True, False, cores=cores, ascats=ascats, ploidy_coeff=ploidy_coeff).dropna()
    cnadf = pd.DataFrame(dict(zip(cnas_filtered.index, cnas_filtered.values))).T
    cnadf.to_csv(output, sep='\t')
    print(cnadf)

def process_somatic_variants(kwargs):
    snv_file = f"{kwargs['somatic_variants']}/{kwargs['pid']}.csv"
    output = kwargs["output"]
    cores = int(kwargs.get("cores", 0))
    refgenome = kwargs.get("refgen", "GRCh38")
    tumortype = kwargs.get("tumortype", "HGSOC")

    ascats = pd.read_csv(kwargs["ascatestimates"], sep="\t", encoding='utf-8')
    snvs = pd.read_csv(snv_file, sep="\t", encoding='utf-8')
    cnas = dd.read_csv(f"{kwargs['cn_annotations']}/{kwargs['pid']}*.csv", sep="\t").compute().dropna(subset=['Gene', 'nMinor', 'nMajor', 'LOHstatus'])

    samples = cnas['sample'].drop_duplicates()
    snvs_filtered = df_apply(snvs, sva.filter_and_classify_snvs, None, None, True, True, cores=cores, cnas=cnas, ascats=ascats, samples=samples)

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
        snv_grp = snv_grps.get_group(gname)
        snvs_filtered = df_apply(snv_grp, sva.post_filter_and_classify_snvs_by_sample, None, None, True, True, cores=cores, cnas=cna_grp, ascats=ascats)
        snvdf = pd.DataFrame(snvs_filtered)
        snvdf.to_csv(output_file, sep='\t', header=False)
    except Exception as e:
        print(f"Error processing sample {gname}: {e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=str,
        required=False,
        help="Path to CSV or TSV file",
    )
    parser.add_argument(
        "--somatic_variants",
        type=str,
        required=False,
        help="Path to CSV or TSV file to import data",
    )
    parser.add_argument(
        "--copy_number_alterations",
        type=str,
        required=False,
        help="Path to CSV or TSV file to import data",
    )
    parser.add_argument(
        "--clinical_data",
        type=str,
        required=False,
        help="Path to CSV or TSV file to import data",
    )
    parser.add_argument(
        "--oncokb_actionable_targets",
        type=str,
        required=False,
        help="Path to CSV or TSV file to import data",
    )
    parser.add_argument(
        "--ascatestimates",
        type=str,
        required=False,
        help="Path to CSV or TSV file to import data",
    )

    parser.add_argument(
        "--alphamissenses",
        type=str,
        required=False,
        help="Path to CSV or TSV file to import data",
    )

    parser.add_argument(
        "--filter",
        type=str,
        required=False,
        help="Field to filter data",
    )
    parser.add_argument(
        "--contains",
        type=str,
        required=False,
        help="Value to filter data",
    )
    parser.add_argument(
        "--equal",
        type=str,
        required=False,
        help="Value to filter data",
    )
    parser.add_argument(
        "--notequal",
        type=str,
        required=False,
        help="Value to filter data",
    )
    parser.add_argument(
        "--gt",
        type=str,
        required=False,
        help="Greater than value to filter data",
    )
    parser.add_argument(
        "--lt",
        type=str,
        required=False,
        help="Less than value to filter data",
    )
    parser.add_argument(
        "--deletesnvs",
        action='store_true',
        required=False,
        help="Remove variants from database",
    )
    parser.add_argument(
        "--deletecnas",
        action='store_true',
        required=False,
        help="Remove variants from database",
    )
    parser.add_argument(
        "--noheader",
        action='store_true',
        required=False,
        help="Use hardcoded header",
    )
    parser.add_argument(
        "--cnfilter",
        action='store_true',
        required=False,
        help="Filter by copynumber threshold",

    )
    parser.add_argument(
        "--all",
        action='store_true',
        required=False,
        help="Pass all variants",

    )
    parser.add_argument(
        "--tumortype",
        type=str,
        required=False,
        help="tumortype identifier OncoKB: HGSOC = CGI: OVE",
    )

    parser.add_argument(
        "--refgen",
        type=str,
        required=False,
        help="Reference genome version in format GRCh38",
    )
    parser.add_argument(
        "--cores",
        type=str,
        required=False,
        help="",
    )
    parser.add_argument(
        "--cn_annotations",
        type=str,
        required=False,
        help="Filtered and annotated CNAs",
    )
    parser.add_argument(
        "--post_annotate",
        action='store_true',
        required=False,
        help="",

    )
    parser.add_argument(
        "--ploidy_threshold",
        type=str,
        required=False,
        help="Threshold for filtering by ploidy",
    )

    parser.add_argument(
        "--pid",
        type=str,
        required=False,
        help="Patient ID",
    )

    args = parser.parse_args()

    main(**vars(args))