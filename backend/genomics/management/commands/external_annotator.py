import pandas as pd
import argparse
import time

from oncokb_annotator import query_oncokb_cnas_to_csv, query_oncokb_somatic_mutations
from cgi_annotator import generate_temp_cgi_query_files, query_cgi_job, launch_cgi_job_with_mulitple_variant_types

def main(**kwargs):

# OncoKB queries
    #USAGE: docker compose run --rm backend sh -c "python manage.py genomic_db_query_utils --oncokbcna --geneid=ENSG00000230280 --patientid=1"
    # if kwargs["import_to_django"]:
    #     if kwargs["copy_number_alterations"]:
    #         cna_data = pd.read_csv(kwargs.get("copy_number_alterations", ""), sep="\t")
    #         cnadf_to_django_model(cna_data)
    #     if kwargs["somatic_variants"]:
    #         snvs = pd.read_csv(kwargs["somatic_variants"], sep="\t")
    #         snvs_to_django_model(snvs)



    # if kwargs["oncokbsnv"] and kwargs["proteinchange"] and kwargs["cohortcode"]: # Query all exonic mutations for given patient
    #     pid = map_cohort_code_to_patient_id(kwargs["cohortcode"])
    #     snvs = get_actionable_snvs_by_aaChangeRefGene(pid)
    #     query_oncokb_somatic_mutations(snvs, kwargs["cancer"])
    # if kwargs["oncokbsnv"] and kwargs["exonic"] and kwargs["cohortcode"]: # Query all exonic mutations for given patient
    #     targets = ActionableTarget.objects.all()
    #     pid = map_cohort_code_to_patient_id(kwargs["cohortcode"])
    #     snvs = get_all_exonic_snvs_of_patient(pid, targets)
    #     query_oncokb_somatic_mutations(snvs, kwargs["cancer"])
    #     #chunks = [snvs[x:x+10] for x in range(0, len(snvs), 10)]
    #     #for c in chunks:
    #     #    query_oncokb_somatic_mutations(c, kwargs["cancer"])
    #     #    time.sleep(1)
    # if kwargs["oncokbcna"] and kwargs["targetall"] and kwargs["cohortcode"]:
    #     targets = ActionableTarget.objects.all()
    #     pid = map_cohort_code_to_patient_id(kwargs["cohortcode"])
    #     cnas = get_cnas_by_gene_list(pid, targets)
    #     query_oncokb_cnas(cnas, kwargs["cancer"])
    #
    # if kwargs["oncokbcna"] and kwargs["actionable"] and kwargs["cohortcode"]:
    #     pid = map_cohort_code_to_patient_id(kwargs["cohortcode"])
    #     targets = ActionableTarget.objects.all()
    #     cnas = get_cnas_by_cn_and_ploidy(pid, kwargs["cnthr"], targets)
    #     query_oncokb_cnas(cnas, kwargs["cancer"])
    #
    # if kwargs["oncokbsnv"] and kwargs["actionable"] and kwargs["cohortcode"] and not kwargs["single"]:
    #     targets = ActionableTarget.objects.all()
    #     pid = map_cohort_code_to_patient_id(kwargs["cohortcode"])
    #     snvs = get_snvs_by_gene_list(pid, targets)
    #     #chunks = [snvs[x:x + 30] for x in range(0, len(snvs), 30)]
    #     #for c in chunks:
    #     query_oncokb_somatic_mutations(snvs, kwargs["cancer"])
    #
    # if kwargs["oncokbsnv"] and kwargs["actionable"] and kwargs["cohortcode"] and kwargs["single"]:
    #     targets = ActionableTarget.objects.all()
    #     pid = map_cohort_code_to_patient_id(kwargs["cohortcode"])
    #     snvs = get_snvs_by_gene_list(pid, targets)
    #     for snv in snvs:
    #         query_oncokb_somatic_mutation(snv, kwargs["cancer"])
    #         #time.sleep(1)


    if kwargs["oncokbcna"] and kwargs["copy_number_alterations"] and kwargs["all"]:
        #for pid in ClinicalData.objects.order_by().values('patient_id').distinct():
        #    cnas = cna_annotation.objects.filter(patient_id=pid.get('patient_id'))
        #   if cnas:
        cnas = pd.read_csv(kwargs["copy_number_alterations"], sep="\t")
        cnas['tumorType'] = ""
        cnas['oncogenic'] = ""
        cnas['mutationEffectDescription'] = ""
        cnas['gene_role'] = ""
        cnas['citationPMids'] = ""
        cnas['level_of_evidence'] = ""
        cnas['cgi_level'] = ""
        cnas['geneSummary'] = ""
        cnas['variantSummary'] = ""
        cnas['tumorTypeSummary'] = ""
        cnas['treatments'] = ""

        chunks = [cnas[x:x + 4999] for x in range(0, len(cnas), 5000)]
        for c in chunks:
            query_oncokb_cnas_to_csv(c)


    if kwargs["oncokbsnv"] and kwargs["somatic_variants"] and kwargs["all"]:
        #for pid in ClinicalData.objects.order_by().values('patient_id').distinct():
        #    cnas = cna_annotation.objects.filter(patient_id=pid.get('patient_id'))
        #   if cnas:
        snvs = pd.read_csv(kwargs["somatic_variants"], sep="\t")

        snvs['tumorType'] = ""
        snvs['consequence_okb'] = ""
        snvs['oncogenic'] = ""
        snvs['mutationEffectDescription'] = ""
        snvs['gene_role'] = ""
        snvs['citationPMids'] = ""
        snvs['level_of_evidence'] = ""
        snvs['cgi_level'] = ""
        snvs['geneSummary'] = ""
        snvs['variantSummary'] = ""
        snvs['tumorTypeSummary'] = ""
        snvs['treatments'] = ""

        chunks = [snvs[x:x + 4999] for x in range(0, len(snvs), 5000)]
        for c in chunks:
            query_oncokb_somatic_mutations(c)

    if kwargs["cgiquery"] and kwargs["somatic_variants"] and kwargs["all"]:
        snvs = pd.read_csv(kwargs["somatic_variants"], sep="\t", dtype='string')
        if kwargs["cgijobid"]:
            jobid = kwargs["cgijobid"]
        else:
            generate_temp_cgi_query_files(snvs, None, None)
            jobid = launch_cgi_job_with_mulitple_variant_types("./tmp/snvs.ext",None, None, "OVSE", "hg38").replace('"', '')
        time.sleep(30)
        while query_cgi_job(jobid, snvs) == 0:
            print("Waiting 30 seconds for the next try...")
            time.sleep(30)

    if kwargs["cgiquery"] and kwargs["copy_number_alterations"] and kwargs["all"]:
        cnas = pd.read_csv(kwargs["copy_number_alterations"], sep="\t", dtype='string')
        if kwargs["cgijobid"]:
            jobid = kwargs["cgijobid"]
        else:
            generate_temp_cgi_query_files(None, cnas, None)
            jobid = launch_cgi_job_with_mulitple_variant_types(None, "./tmp/cnas.ext", None, "OVSE", "hg38").replace('"', '')

        time.sleep(30)
        while query_cgi_job(jobid, None, cnas) == 0:
            print("Waiting 30 seconds for the next try...")
            time.sleep(30)

# CGI queries

    #USAGE: docker compose run --rm backend sh -c "python manage.py genomic_db_query_utils --cgiquery --exonic --patientid=1"
    # if kwargs["cgiquery"] and kwargs["exonic"] and kwargs["snv"] and kwargs["cohortcode"]: # Query all exonic mutations for given patient
    #     pid = map_cohort_code_to_patient_id(kwargs["cohortcode"])
    #     snv = get_all_exonic_snvs_of_patient(pid)
    #     generate_temp_cgi_query_files(snv,[],[])
    #     jobid = launch_cgi_job_with_mulitple_variant_types("./tmp/snvs.ext", None, None, kwargs["cancer"], "hg38")
    #     if jobid != 0:
    #         while query_cgi_job(pid, jobid.replace('"', '')) == 0:
    #             print("Waiting 120 seconds for the next try...")
    #             time.sleep(30)
    #  #USAGE: docker compose run --rm backend sh -c "python manage.py genomic_db_query_utils --cgiquery --proteinchange --patientid=1"
    # if kwargs["cgiquery"] and kwargs["proteinchange"] and kwargs["snv"] and kwargs["cohortcode"]: # Query all protein affecting mutations for all patients
    #     pid = map_cohort_code_to_patient_id(kwargs["cohortcode"])
    #     snvs = get_actionable_snvs_by_aaChangeRefGene(kwargs["patientid"])
    #     generate_proteinchange_query_file(snvs)
    #     jobid = launch_cgi_job_with_mulitple_variant_types("./tmp/prot.ext", None, None, kwargs["cancer"], "hg38")
    #     while query_cgi_job(pid, jobid.replace('"', '')) == 0:
    #         print("Waiting 120 seconds for the next try...")
    #         time.sleep(120)

    # May be impossible to query every patient at once from cgi, could be done if distinct genes of every patient mapped to same query file
    # if kwargs["cgiquery"] and kwargs["cna"] and kwargs["targetall"] and kwargs["cohortcode"]:
    #     targets = ActionableTarget.objects.all()
    #     pid = map_cohort_code_to_patient_id(kwargs["cohortcode"])
    #     cnas = get_cnas_by_gene_list(pid, targets)
    #     if cnas:
    #         generate_temp_cgi_query_files([], cnas, [])
    #         jobid = launch_cgi_job_with_mulitple_variant_types(None, "./tmp/cnas.ext", None, kwargs["cancer"], "hg38")
    #         time.sleep(10)
    #         while query_cgi_job(pid, jobid.replace('"', '')) == 0:
    #             print("Waiting 30 seconds for the next try...")
    #             time.sleep(30)
    #
    # if kwargs["cgiquery"] and kwargs["cna"] and kwargs["actionable"] and kwargs["cohortcode"]:
    #     pid = map_cohort_code_to_patient_id(kwargs["cohortcode"])
    #     targets = ActionableTarget.objects.all()
    #     cnas = get_cnas_by_cn_and_ploidy(pid, kwargs["cnthr"], targets)
    #     if cnas:
    #         genfiles = generate_temp_cgi_query_files([], cnas, [])
    #         if genfiles:
    #             jobid = launch_cgi_job_with_mulitple_variant_types(None, "./tmp/cnas.ext", None, kwargs["cancer"], "hg38")
    #             time.sleep(10)
    #             if jobid != 0:
    #                 while query_cgi_job(pid, jobid.replace('"', '')) == 0:
    #                     print("Waiting 30 seconds for the next try...")
    #                     time.sleep(30)
    #         else:
    #             print("No cgi variant files generated!")
    #     else:
    #         print("No CNAs!")
    # if kwargs["cgiquery"] and kwargs["snv"] and kwargs["actionable"] and kwargs["cohortcode"]:
    #     targets = ActionableTarget.objects.all()
    #     pid = map_cohort_code_to_patient_id(kwargs["cohortcode"])
    #     snvs = get_snvs_by_gene_list(pid, targets)
    #     if snvs:
    #         generate_temp_cgi_query_files(snvs, [], [])
    #         jobid = launch_cgi_job_with_mulitple_variant_types("./tmp/snvs.ext", None, None, kwargs["cancer"], "hg38")
    #         time.sleep(10)
    #         while query_cgi_job(pid, jobid.replace('"', '')) == 0:
    #                 print("Waiting 30 seconds for the next try...")
    #                 time.sleep(30)

if __name__ == "__main__":

    def add_arguments(parser):
        parser.add_argument('--snv',  action='store_true', help='')
        parser.add_argument('--patientid', type=str, help='Give patient id (in internal DB)')
        parser.add_argument('--cohortcode', type=str, help='Give cohort code of patient')
        parser.add_argument('--geneid', type=str, help='Give gene id (in format eg. ENSG00000272262, NM_032264 or NBPF3)')
        parser.add_argument('--gene', type=str, help='Give gene name (in format eg. BRCA1)')
        parser.add_argument('--refid', type=str,  help='Give reference id (in dbSNP format eg. rs61769312)')
        parser.add_argument('--cna',  action='store_true', help='')
        parser.add_argument('--cancer', type=str, default='HGSOC', help='HGSOC, OV, CANCER')
        parser.add_argument('--cnatype', type=str, help='AMP,DEL')
        parser.add_argument('--targetall', action='store_true', help='Query all targetable variants from OncoKB actionable target list')
        parser.add_argument('--fusgenes', type=str,  help='Give a list of fusion genes eg. BCR__ABL1,PML__PARA')
        parser.add_argument('--cgijobid', type=str, help='Download results from CGI by jobid')
        parser.add_argument('--cgiquery',  action='store_true', help='Download results from CGI by jobid')
        parser.add_argument('--oncokbcna',  action='store_true',  help='Query OncoKB by gene id from given patient CNA. Input parameter: gene id')
        parser.add_argument('--oncokbsnv', action='store_true',  help='Query OncoKB by genomic location parsed patient SNV. Input parameter: ref id')
        parser.add_argument('--sqlsnvs',  action='store_true', help='')
        parser.add_argument('--exonic',  action='store_true', help='')
        parser.add_argument('--actionable', action='store_true', help='Query by significant copy number threshold. Default is minimum of 8 copies for AMP or minimum ploidy 2, and 0 for DEL. For somatic mutations OncoKB actionable target list is used as reference.')
        parser.add_argument('--proteinchange',  action='store_true', help='')
        parser.add_argument('--genelist', type=str, help='')
        parser.add_argument('--logrthr', type=float, default=1.5, help='logr threshold for filtering by ploidy')
        parser.add_argument('--cnthr', type=int, default=9, help='Copy number threshold')
        parser.add_argument('--single', action='store_true', help='')
        parser.add_argument('--all', action='store_true', help='')
        parser.add_argument('--direct', action='store_true', help='')
        parser.add_argument('--retrieve', action='store_true', help='')
        parser.add_argument('--import_to_django', action='store_true', help='')
        parser.add_argument('--copy_number_alterations', type=str, help='')
        parser.add_argument('--somatic_variants', type=str, help='')


    parser = argparse.ArgumentParser()
    add_arguments(parser)
    args = parser.parse_args()
    main(**vars(args))
