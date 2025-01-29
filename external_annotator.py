import pandas as pd
import argparse
import time

from oncokb_annotator import query_oncokb_cnas_to_csv, query_oncokb_somatic_mutations
from cgi_annotator import generate_temp_cgi_query_files, query_cgi_job, launch_cgi_job_with_mulitple_variant_types

'''
Usage: external_annotator.py [OPTIONS]

Options:
  --oncokbcna                  Query OncoKB for copy number alterations.
  --oncokbsnv                  Query OncoKB for somatic mutations.
  --cgiquery                   Query Cancer Genome Interpreter.
  --cgijobid <str>             Download results from CGI by jobid and apply annotations.
  --copy_number_alterations <str> Path to copy number alterations file.
  --somatic_variants <str>     Path to somatic variants file.
  --output <str>               Path to output file.


Examples:
  python external_annotator.py --oncokbcna --copy_number_alterations path/to/cnas.tsv --output path/to/output
  python external_annotator.py --oncokbsnv --somatic_variants path/to/snvs.tsv --output path/to/output
  python external_annotator.py --cgiquery --somatic_variants path/to/snvs.tsv --output path/to/output
  python external_annotator.py --cgiquery --copy_number_alterations path/to/cnas.tsv --output path/to/output
  python external_annotator.py --cgiquery --cgijobid <jobid> --somatic_variants path/to/snvs.tsv --output path/to/output
  python external_annotator.py --cgiquery --cgijobid <jobid> --copy_number_alterations path/to/cnas.tsv --output path/to/output
'''

def main(**kwargs):

    output = kwargs.get("output", ".")
    if kwargs["oncokbcna"] and kwargs["copy_number_alterations"]:

        cnas = pd.read_csv(kwargs["copy_number_alterations"], sep="\t")
        cnas['oncogenic'] = ""
        cnas['mutationEffectDescription'] = ""
        cnas['gene_role'] = ""
        cnas['citationPMids'] = ""
        cnas['level_of_evidence'] = ""
        cnas['cgi_level'] = ""
        cnas['geneSummary'] = ""
        cnas['variantSummary'] = ""
        cnas['tumorTypeSummary'] = ""

        # Query in chunks of 5000
        chunks = [cnas[x:x + 4999] for x in range(0, len(cnas), 5000)]
        i = 0
        for c in chunks:
            i += 1
            query_oncokb_cnas_to_csv(c, output, i)


    if kwargs["oncokbsnv"] and kwargs["somatic_variants"]:
        snvs = pd.read_csv(kwargs["somatic_variants"], sep="\t")

        snvs['consequence'] = ""
        snvs['oncogenic'] = ""
        snvs['mutationEffectDescription'] = ""
        snvs['gene_role'] = ""
        snvs['citationPMids'] = ""
        snvs['level_of_evidence'] = ""
        snvs['cgi_level'] = ""
        snvs['geneSummary'] = ""
        snvs['variantSummary'] = ""
        snvs['tumorTypeSummary'] = ""

        # Query in chunks of 5000
        chunks = [snvs[x:x + 4999] for x in range(0, len(snvs), 5000)]
        i = 0
        for c in chunks:
            i += 1
            query_oncokb_somatic_mutations(c, output, i)

    if kwargs["cgiquery"] and kwargs["somatic_variants"]:
        snvs = pd.read_csv(kwargs["somatic_variants"], sep="\t", dtype='string')

        if kwargs["cgijobid"]:
            jobid = kwargs["cgijobid"]
        else:
            generate_temp_cgi_query_files(snv_annotations=snvs)
            jobid = launch_cgi_job_with_mulitple_variant_types(mutations_file="./tmp/snvs.ext", cancer_type="OVSE", reference="hg38").replace('"', '')
        time.sleep(30)
        while query_cgi_job(jobid, output, snv_annotations=snvs) == 0:
            print("Waiting 30 seconds for the next try...")
            time.sleep(30)

    if kwargs["cgiquery"] and kwargs["copy_number_alterations"]:
        cnas = pd.read_csv(kwargs["copy_number_alterations"], sep="\t", dtype='string')

        if kwargs["cgijobid"]:
            jobid = kwargs["cgijobid"]
        else:
            generate_temp_cgi_query_files(cna_annotations=cnas)
            jobid = launch_cgi_job_with_mulitple_variant_types(cnas_file="./tmp/cnas.ext", cancer_type="OVSE", reference="hg38").replace('"', '')

        time.sleep(30)
        while query_cgi_job(jobid, output, cna_annotations=cnas) == 0:
            print("Waiting 30 seconds for the next try...")
            time.sleep(30)

if __name__ == "__main__":

    def add_arguments(parser):
        parser.add_argument('--oncokbcna', action='store_true', help='Query OncoKB for copy number alterations')
        parser.add_argument('--oncokbsnv', action='store_true', help='Query OncoKB for somatic mutations')
        parser.add_argument('--cgiquery', action='store_true', help='Query Cancer Genome Interpreter')
        parser.add_argument('--cgijobid', type=str, help='Download results from CGI by jobid')
        parser.add_argument('--copy_number_alterations', type=str, help='Path to copy number alterations file')
        parser.add_argument('--somatic_variants', type=str, help='Path to somatic variants file')
        parser.add_argument('--output', type=str, default=".", help='Path to output directory for annotated files')



    parser = argparse.ArgumentParser()
    add_arguments(parser)
    args = parser.parse_args()
    main(**vars(args))
