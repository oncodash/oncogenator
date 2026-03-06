import zipfile
from enum import Enum

import pandas as pd
from utils import *
import httpx
import io
import urllib3
from config import CGI_LOGIN, CGI_TOKEN, CGI_API_URL, CGI_DEFAULT_CANCER_TYPE, CGI_DEFAULT_REFERENCE

class cna_alt_to_cgi(Enum):
    AMPLIFICATION = "AMP"
    DELETION = "DEL"
    def __str__(self):
        return str(self.value)

class cgi2oncokb_level(Enum):
    A = "LEVEL_1"
    B = "LEVEL_2"
    C = "LEVEL_3A"
    D = "LEVEL_3B"
    E = "LEVEL_4"
    R1 = "LEVEL_R1"
    R2 = "LEVEL_R2"
    def __str__(self):
        return str(self.value)

def map_cgi_evidence(biomarker):
    """
        Map CGI evidence to OncoKB levels.

        Parameters:
        biomarker (Series): A Series containing biomarker data.

        Returns:
        str: Mapped OncoKB level.
    """
    evidence = biomarker['Evidence']
    response = biomarker['Response']
    if pd.isna(evidence):
        return None
    if response == "Responsive":
        return cgi2oncokb_level[evidence].value
    if response == "Resistant":
        if cgi2oncokb_level[evidence] in ["LEVEL_1", "LEVEL_2"]:
            return cgi2oncokb_level["R1"].value
        if cgi2oncokb_level[evidence] in ["LEVEL_3A", "LEVEL_3B", "LEVEL_4"]:
            return cgi2oncokb_level["R2"].value
    return None

def handle_treatments_cgi(row, alt_type, alteration):
    """
        Handle treatments from CGI data.

        Parameters:
        row (Series): A Series containing treatment data.
        alt_type (str): Alteration type.
        alteration (str): Alteration description.

        Returns:
        Series: A Series containing treatment information.
    """

    drugs = row['Drugs']
    pmids = row['Source']
    approvedIndications = row['Biomarker']
    tumortype = row['Tumor type']
    level = map_cgi_evidence(row)
    description = ""
    return pd.Series({
        'alteration_type': alt_type,
        'alteration': alteration,
        'approvedIndications': approvedIndications,
        'description': description,
        'treatment': drugs,
        'level_of_evidence': level,
        'cgi_level': handle_string_field(row['Evidence'])+"("+handle_string_field(row['Response'])+")",
        'citations': pmids,
        'tumorType': tumortype
    })

def generate_cgi_cna_file_from_list(genelist):
    """
        Launch a CGI job with multiple variant types.

        Parameters:
        mutations_file (str): Path to the mutation file.
        cnas_file (str): Path to the CNAs file.
        transloc_file (str): Path to the translocation file.
        cancer_type (str): Type of cancer.
        reference (str): Reference genome.

        Returns:
        str: Job ID if the request is successful, otherwise 0.
    """
    header = "gene\tcna\n"
    with open("./tmp/cnas.ext", "w") as file2:
        file2.write(header)
        genes = genelist
        for gene in genes:
            row = gene + '\tAMP\n'
            print(row)
            file2.write(row)
        file2.close()

def launch_cgi_job_with_mulitple_variant_types(mutations_file=None, cnas_file=None, transloc_file=None, cancer_type=None, reference=None):
    """
        This function launches a CGI (Cancer Genome Interpreter) job with multiple variant types,
        using the CGI API. It takes in mutation, cnas, and translocation files, cancer type, and
        reference as input, and returns a job ID if the request is successful.

        Args:
        mutations_file (str): The path to the mutation file.
        cnas_file (str): The path to the cnas file.
        transloc_file (str): The path to the translocation file.
        cancer_type (str): The type of cancer.
        reference (str): The reference genome.

        Returns:
        jobid (str): The job ID if the request is successful.

        Raises:
        None.
        """
    
    # Use defaults from config if not provided
    if cancer_type is None:
        cancer_type = CGI_DEFAULT_CANCER_TYPE
    if reference is None:
        reference = CGI_DEFAULT_REFERENCE

    request_url = CGI_API_URL
    login = CGI_LOGIN
    token = CGI_TOKEN

    print("Request CGI")
    # CGI api requires every type mutation files to be provided
    headers = {
        'Authorization': login+' '+token
    }

    if cnas_file:
        payload = {
            'cancer_type': cancer_type,
            'title': 'Title',
            'reference': reference,
            'cnas': ('cnas.ext', open(cnas_file, 'rb').read(), 'application/octet-stream')
        }
    if mutations_file:
        payload = {
            'cancer_type': cancer_type,
            'title': 'Title',
            'reference': reference,
            'mutations': ('somatic_mutations.ext', open(mutations_file, 'rb').read(), 'application/octet-stream'),
        }

    # Make the POST request using multipart/form-data with the files parameter
    http = urllib3.PoolManager()

    # Make the POST request using multipart/form-data with the files parameter
    response = http.request(
        'POST',
        'https://www.cancergenomeinterpreter.org/api/v1',
        fields=payload,
        headers=headers,
        multipart_boundary="----WebKitFormBoundary7MA4YWxkTrZu0gW",
        preload_content=False  # Set preload_content to False to allow streaming the files
    )

    if (response.status == 200):

        jobid = response.data.decode("utf-8")
        print(jobid)
        return jobid

    else:
        print("[ERROR] Unable to request. Response: ", print(response.data))
        return 0


def query_cgi_job(jobid, output, somatic_mutation_annotations: pd.DataFrame = None, cna_annotations: pd.DataFrame = None, mode="x"):
    """
    Query the CGI API with a job ID and save the results to the database.

    Parameters:
    jobid (str): The job ID for the CGI job to query.
    somatic_mutation_annotations (DataFrame): DataFrame containing somatic_mutation annotations.
    cna_annotations (DataFrame): DataFrame containing CNA annotations.

    Returns:
    int: 1 if successful, otherwise 0.
    """
    request_url = CGI_API_URL + "/"
    print("Request CGI job by id")

    cgilogin = CGI_LOGIN
    cgitoken = CGI_TOKEN

    headers = {
        'Authorization': cgilogin + ' ' + cgitoken
    }
    payload = {'action': 'download'}
    # response = httpx.request("GET",request_url+jobid, headers=headers, fields=payload)
    response = httpx.get(request_url + jobid, params=payload, headers=headers, timeout=None)
    print("CGI response status code: ", response.status_code)
    if response.status_code == 200:
        z = zipfile.ZipFile(io.BytesIO(response.content))
        treatmentsdf = None
        cgi_somatic_mutationdf = None
        cgi_cnadf = None
        treatments = []

        for fn in z.namelist():
            base_name = fn.rsplit("/", 1)[-1]
            if base_name not in {"alterations.tsv", "cna_analysis.tsv", "biomarkers.tsv"}:
                continue
            with z.open(fn) as f:
                df = pd.read_csv(f, sep="\t")

            if base_name == "alterations.tsv":
                cgi_somatic_mutationdf = df
            if base_name == "cna_analysis.tsv":
                cgi_cnadf = df
            if base_name == "biomarkers.tsv":
                treatmentsdf = df

        if treatmentsdf is None:
            print("No CGI biomarkers.tsv found for job id: " + str(jobid))
            return 0

        cgi_cna_lookup = {}
        if isinstance(cgi_cnadf, pd.DataFrame) and {'sample', 'driver', 'gene_role'}.issubset(cgi_cnadf.columns):
            cgi_cna_lookup = cgi_cnadf.drop_duplicates(subset=['sample'], keep='first').set_index('sample')[['driver', 'gene_role']].to_dict(orient='index')

        cgi_snv_lookup = {}
        if isinstance(cgi_somatic_mutationdf, pd.DataFrame) and {'CGI-Sample ID', 'CGI-Consequence', 'CGI-Oncogenic Summary', 'CGI-Oncogenic Prediction'}.issubset(cgi_somatic_mutationdf.columns):
            cgi_snv_lookup = cgi_somatic_mutationdf.drop_duplicates(subset=['CGI-Sample ID'], keep='first').set_index('CGI-Sample ID')[['CGI-Consequence', 'CGI-Oncogenic Summary', 'CGI-Oncogenic Prediction']].to_dict(orient='index')

        bioms = treatmentsdf.loc[treatmentsdf['Match'] == 'YES'].copy()
        bioms['sample_id'] = bioms['Sample ID'].apply(handle_string_field)

        cna_bioms = bioms.loc[bioms['sample_id'].str.startswith('CNA:', na=False)].copy()
        if not cna_bioms.empty:
            cna_split = cna_bioms['sample_id'].str.split(':', n=2, expand=True)
            cna_bioms['hugoSymbol'] = cna_split[1]
            cna_bioms['cna_alteration'] = cna_split[2]
            cna_bioms['alteration_key'] = cna_bioms['hugoSymbol'] + ':' + cna_bioms['cna_alteration']
            
            for _, row in cna_bioms.iterrows():
                print("Appending treatment for CNA biomarker: ", row['alteration_key'])
                treatments.append(handle_treatments_cgi(row, 'CNA', row['alteration_key']))

            if isinstance(cna_annotations, pd.DataFrame) and {'hugoSymbol', 'alteration', 'oncogenic'}.issubset(cna_annotations.columns):
                cna_updates = cna_annotations.loc[
                    (cna_annotations['oncogenic'] == "Unknown") | (cna_annotations['oncogenic'].isna()),
                    ['hugoSymbol', 'alteration']
                ].reset_index()

                if not cna_updates.empty:
                    cna_updates = cna_updates.merge(
                        cna_bioms[['hugoSymbol', 'cna_alteration', 'sample_id']].drop_duplicates(subset=['hugoSymbol', 'cna_alteration'], keep='first'),
                        left_on=['hugoSymbol', 'alteration'],
                        right_on=['hugoSymbol', 'cna_alteration'],
                        how='left'
                    )
                    cna_lookup_df = pd.DataFrame.from_dict(cgi_cna_lookup, orient='index').reset_index().rename(columns={'index': 'sample_id'})
                    cna_updates = cna_updates.merge(cna_lookup_df, on='sample_id', how='left')

                    cna_valid = cna_updates.loc[cna_updates['driver'].notna()]
                    if not cna_valid.empty:
                        cna_annotations.loc[cna_valid['index'], 'oncogenic'] = cna_valid['driver'].apply(handle_string_field).values
                        cna_annotations.loc[cna_valid['index'], 'gene_role'] = cna_valid['gene_role'].apply(handle_string_field).values

        snv_bioms = bioms.loc[bioms['sample_id'].str.startswith('SNV:', na=False)].copy()
        if not snv_bioms.empty:
            snv_split = snv_bioms['sample_id'].str.split(':', n=5, expand=True)
            snv_bioms['alteration_key'] = (
                snv_split[1].astype(str) + ':' +
                snv_split[2].astype(str) + ':' +
                snv_split[3].astype(str) + ':' +
                snv_split[4].astype(str) + ':' +
                snv_split[5].astype(str)
            )
            for _, row in snv_bioms.iterrows():
                treatments.append(handle_treatments_cgi(row, 'SNV', row['alteration_key']))

            if isinstance(somatic_mutation_annotations, pd.DataFrame) and {'alteration', 'oncogenic'}.issubset(somatic_mutation_annotations.columns):
                snv_updates = somatic_mutation_annotations.loc[
                    (somatic_mutation_annotations['oncogenic'] == "Unknown") | (somatic_mutation_annotations['oncogenic'].isna()),
                    ['alteration']
                ].reset_index()

                if not snv_updates.empty:
                    snv_updates = snv_updates.merge(
                        snv_bioms[['alteration_key', 'sample_id']].drop_duplicates(subset=['alteration_key'], keep='first'),
                        left_on='alteration',
                        right_on='alteration_key',
                        how='left'
                    )
                    snv_lookup_df = pd.DataFrame.from_dict(cgi_snv_lookup, orient='index').reset_index().rename(columns={'index': 'sample_id'})
                    snv_updates = snv_updates.merge(snv_lookup_df, on='sample_id', how='left')

                    snv_valid = snv_updates.loc[snv_updates['CGI-Oncogenic Summary'].notna()]
                    if not snv_valid.empty:
                        somatic_mutation_annotations.loc[snv_valid['index'], 'consequence'] = snv_valid['CGI-Consequence'].apply(handle_string_field).values
                        somatic_mutation_annotations.loc[snv_valid['index'], 'oncogenic'] = snv_valid['CGI-Oncogenic Summary'].apply(handle_string_field).values
                        somatic_mutation_annotations.loc[snv_valid['index'], 'gene_role'] = snv_valid['CGI-Oncogenic Prediction'].apply(handle_string_field).values

        if isinstance(somatic_mutation_annotations, pd.DataFrame):
            somatic_mutation_annotations.to_csv(output, mode=mode, index=False, sep="\t", columns=['patient_id', 'sample_id', 'alteration', 'hugoSymbol', 'ensembl_id', 'tumorType', 'consequence', 'oncogenic', 'mutationEffectDescription', 'gene_role', 'citationPMids', 'geneSummary', 'variantSummary', 'tumorTypeSummary'])
            trdf = pd.DataFrame(treatments)
            trdf.to_csv("treatments.csv", mode="a", index=False, sep="\t")

        if isinstance(cna_annotations, pd.DataFrame):
            cna_annotations.to_csv(output, mode=mode, index=False, sep="\t", columns=['patient_id', 'sample_id', 'alteration', 'hugoSymbol', 'ensembl_id', 'tumorType', 'oncogenic', 'mutationEffectDescription', 'gene_role', 'citationPMids', 'geneSummary', 'variantSummary', 'tumorTypeSummary'])
            trdf = pd.DataFrame(treatments)
            trdf.to_csv("treatments.csv", mode="a", index=False, sep="\t")

        return 1
    else:
    #print(response.status_code)
        print("No CGI results available for job id: "+str(jobid))
        return 0
    
def query_cgi_job_old(jobid, output, somatic_mutation_annotations: pd.DataFrame = None, cna_annotations: pd.DataFrame = None, mode="x"):
    """
    Query the CGI API with a job ID and save the results to the database.

    Parameters:
    jobid (str): The job ID for the CGI job to query.
    somatic_mutation_annotations (DataFrame): DataFrame containing somatic_mutation annotations.
    cna_annotations (DataFrame): DataFrame containing CNA annotations.

    Returns:
    int: 1 if successful, otherwise 0.
    """
    request_url = CGI_API_URL + "/"
    print("Request CGI job by id")

    cgilogin = CGI_LOGIN
    cgitoken = CGI_TOKEN

    headers = {
        'Authorization': cgilogin + ' ' + cgitoken
    }
    payload = {'action': 'download'}
    # response = httpx.request("GET",request_url+jobid, headers=headers, fields=payload)
    response = httpx.get(request_url + jobid, params=payload, headers=headers, timeout=None)

    if response.status_code == 200:
        z = zipfile.ZipFile(io.BytesIO(response.content))
        fnames = z.namelist()
        treatmentsdf = None
        cgi_somatic_mutationdf = None
        cgi_cnadf = None
        treatments = []

        for fn in fnames:
            z.extract(fn)
            df = pd.read_csv(fn, sep="\t")
            print(fn)
            print(df)

            # Mutation response
            # ['Input ID', 'CHROMOSOME', 'POSITION', 'REF', 'ALT', 'chr', 'pos', 'ref','alt', 'ALT_TYPE', 'STRAND', 'CGI-Sample ID', 'CGI-Gene', 'CGI-Protein Change', 'CGI-Oncogenic Summary', 'CGI-Oncogenic Prediction', 'CGI-External oncogenic annotation','CGI-Mutation', 'CGI-Consequence', 'CGI-Transcript', 'CGI-STRAND', 'CGI-Type', 'CGI-HGVS', 'CGI-HGVSc', 'CGI-HGVSp']

            if fn == "alterations.tsv":
                cgi_somatic_mutationdf = df
            if fn == "cna_analysis.tsv":
                cgi_cnadf = df
            if fn == "biomarkers.tsv":
                treatmentsdf = df

        bioms = treatmentsdf.loc[treatmentsdf['Match'] == 'YES']
        i = 0
        for index, biom in bioms.iterrows():
            # TODO: identify CNA and somatic_mutations from ID and handle separately
            id = handle_string_field(biom["Sample ID"])
            idsplit = id.split(":")
            print(id)
            if idsplit[0] == "CNA":
                alteration = idsplit[1]+":"+idsplit[2]
                treatment = handle_treatments_cgi(biom, 'CNA', alteration)
                print(treatment)
                treatments.append(treatment)
                updatedf = cna_annotations.loc[
                    (((cna_annotations['oncogenic'] == "Unknown") |
                      (cna_annotations['oncogenic'].isna() == True)) & (
                             cna_annotations['hugoSymbol'] == idsplit[1]) & (
                             cna_annotations['alteration'] == idsplit[2]))]
                print(len(updatedf))

                for indxs, row in updatedf.iterrows():
                    i += 1
                    cgi_cna = cgi_cnadf.loc[cgi_cnadf['sample'] == id].iloc[0]
                    cna_annotations.at[indxs, 'oncogenic'] = handle_string_field(cgi_cna["driver"])
                    cna_annotations.at[indxs, 'gene_role'] = handle_string_field(cgi_cna["gene_role"]),

            if idsplit[0] == "SNV":
                hugoSymbol = idsplit[1]
                chromosome = str(idsplit[2])
                position = int(idsplit[3])
                reference_allele = str(idsplit[4])
                sample_allele = str(idsplit[5])
                alteration = hugoSymbol + ":" + chromosome + ":" + str(
                    position) + ":" + reference_allele + ":" + sample_allele

                treatment = handle_treatments_cgi(biom, 'SNV', alteration)
                #print(treatment)
                treatments.append(treatment)
                
                # TODO: try update only if oncokb oncogenic result is None e.g. not known by oncokb
                print(alteration, somatic_mutation_annotations['alteration'], somatic_mutation_annotations['oncogenic'])
                #updatedf = somatic_mutation_annotations.loc[(((somatic_mutation_annotations['oncogenic'] == "Unknown") | (somatic_mutation_annotations['oncogenic'].isna() == True)) & somatic_mutation_annotations['alteration'] == alteration)]
                updatedf = somatic_mutation_annotations.loc[((somatic_mutation_annotations['oncogenic'] == "Unknown") | (somatic_mutation_annotations['oncogenic'].isna() == True)) & (somatic_mutation_annotations['alteration'] == alteration)]
                print("somatic_mutation updatedf:"+str(len(updatedf)))

                for indxs, row in updatedf.iterrows():
                    cgi_somatic_mutation = cgi_somatic_mutationdf.loc[cgi_somatic_mutationdf['CGI-Sample ID'] == id].iloc[0]
                    somatic_mutation_annotations.at[indxs, 'consequence'] = handle_string_field(cgi_somatic_mutation["CGI-Consequence"]),
                    somatic_mutation_annotations.at[indxs, 'oncogenic'] = handle_string_field(cgi_somatic_mutation["CGI-Oncogenic Summary"])
                    somatic_mutation_annotations.at[indxs, 'gene_role'] = handle_string_field(cgi_somatic_mutation["CGI-Oncogenic Prediction"]),

        if isinstance(somatic_mutation_annotations, pd.DataFrame):
            somatic_mutation_annotations.to_csv(output, mode=mode, index=False, sep="\t", columns=['patient_id', 'sample_id', 'alteration', 'hugoSymbol', 'ensembl_id', 'tumorType', 'consequence', 'oncogenic', 'mutationEffectDescription', 'gene_role', 'citationPMids', 'geneSummary', 'variantSummary', 'tumorTypeSummary'])
            trdf = pd.DataFrame(treatments)
            trdf.to_csv("treatments.csv", mode="a", index=False, sep="\t")

        if isinstance(cna_annotations, pd.DataFrame):
            cna_annotations.to_csv(output, mode=mode, index=False, sep="\t", columns=['patient_id', 'sample_id', 'alteration', 'hugoSymbol', 'ensembl_id', 'tumorType', 'oncogenic', 'mutationEffectDescription', 'gene_role', 'citationPMids', 'geneSummary', 'variantSummary', 'tumorTypeSummary'])
            trdf = pd.DataFrame(treatments)
            trdf.to_csv("treatments.csv", mode="a", index=False, sep="\t")

        return 1
    else:
    #print(response.status_code)
        print("No CGI results available for job id: "+str(jobid))
        return 0

def generate_cgi_cna_file_from_list(genelist):
    header = "gene\tcna\n"
    with open("./tmp/cnas.ext", "w") as file2:
        file2.write(header)
        genes = genelist
        for gene in genes:
            row = gene + '\tAMP\n'
            print(row)
            file2.write(row)
        file2.close()

def generate_temp_cgi_query_files(somatic_mutation_annotations: pd.DataFrame = None, cna_annotations: pd.DataFrame = None, translocs: pd.DataFrame = None, append_to_annotations: bool = True):
    """
        Generate temporary CGI query files from annotations.

        Parameters:
        somatic_mutation_annotations (DataFrame): DataFrame containing somatic_mutation annotations.
        cna_annotations (DataFrame): DataFrame containing CNA annotations.
        translocs (DataFrame): DataFrame containing translocation data.
    """
    header = "chr\tpos\tref\talt\tsample\n"
    try:
        if isinstance(somatic_mutation_annotations, pd.DataFrame):
            if append_to_annotations:
                with open("./tmp/somatic_mutations.ext", "w") as file1:
                    file1.write(header)

                    uniques = somatic_mutation_annotations[['alteration']].drop_duplicates()
                    for indx, somatic_mutation in uniques.iterrows():
                        id = "SNV:"+somatic_mutation['alteration']
                        alt_split = somatic_mutation['alteration'].split(':')
                        #print(alt_split)
                        row = alt_split[1]+'\t'+alt_split[2]+'\t'+alt_split[3]+'\t'+alt_split[4]+'\t'+id+'\n'
                        #row = somatic_mutation['chromosome']+'\t'+str(somatic_mutation['position'])+'\t'+somatic_mutation['reference_allele']+'\t'+somatic_mutation['sample_allele']+'\t'+id+'\n' #+'\t'+cryptocode.encrypt(somatic_mutation.samples, settings.CRYPTOCODE)+'\n'

                        file1.write(row)
                    file1.close()
            else:
                with open("./tmp/somatic_mutations.ext", "w") as file1:
                    file1.write(header)

                    uniques = somatic_mutation_annotations[['hugoSymbol', 'chromosome', 'position', 'reference_allele', 'sample_allele', 'tumorType', 'referenceGenome']].drop_duplicates()
                    for indx, somatic_mutation in uniques.iterrows():
                        id = "SNV:"+somatic_mutation['hugoSymbol']+':'+somatic_mutation['chromosome']+':'+str(somatic_mutation['position'])+':'+somatic_mutation['reference_allele']+':'+somatic_mutation['sample_allele']
                        row = somatic_mutation['chromosome']+'\t'+str(somatic_mutation['position'])+'\t'+somatic_mutation['reference_allele']+'\t'+somatic_mutation['sample_allele']+'\t'+id+'\n' #+'\t'+cryptocode.encrypt(somatic_mutation.samples, settings.CRYPTOCODE)+'\n'
                        file1.write(row)
                    file1.close()

        if isinstance(cna_annotations, pd.DataFrame):
            header = "gene\tcna\tsample\n"
            with open("./tmp/cnas.ext", "w") as file2:
                file2.write(header)
                print("CNA annotations:")
                print(cna_annotations)
                uniques = cna_annotations[['hugoSymbol', 'alteration', 'tumorType']].drop_duplicates()
                print(type(uniques))
                for indx, cna in uniques.iterrows():
                    print(cna)
                    id = "CNA:"+str(cna['hugoSymbol']) + ':' + str(cna['alteration'])
                    row = cna['hugoSymbol']+'\t'+cna_alt_to_cgi[cna['alteration']].value+'\t'+id+'\n'
                    file2.write(row)
                file2.close()

    except Exception as e:
        print(f"Unexpected {e=}, {type(e)=}")
        raise
    return 1
