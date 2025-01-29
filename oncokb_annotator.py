from utils import *
import pandas as pd

import json
import httpx

ONCOKB_TOKEN = ""

def handle_treatments_oncokb(jsondata, alt_type, alteration):
    """
        Handle treatments from OncoKB data.

        Parameters:
        jsondata (list): List of dictionaries containing treatment data.
        alt_type (str): Alteration type.
        alteration (str): Alteration description.

        Returns:
        list: A list of Series containing treatment information.
    """
    treatments = []
    for row in jsondata:
        drugs = ""
        print(row)
        if row["drugs"]:
            print(row["drugs"])
            for drug in row["drugs"]:
                print(drug)
                drugs += drug["drugName"]+";"
        pmids = ";".join(row['pmids'])
        approvedIndications = ";".join(row['approvedIndications'])
        tumortype = row['levelAssociatedCancerType']['mainType']['name']
        level = row['level']
        description = row['description']
        treatments.append(pd.Series({
            'alteration_type': alt_type,
            'alteration': alteration,
            'approvedIndications': approvedIndications,
            'description': description,
            'treatment': drugs,
            'level_of_evidence': level,
            'cgi_level':"",
            'citations': pmids,
            'tumorType': tumortype
        }))
    return treatments

def handle_drugs_field(jsondata):
    """
        Handle the drugs field from OncoKB data.

        Parameters:
        jsondata (list): List of dictionaries containing drug data.

        Returns:
        str: A semicolon-separated string of drug names.
    """
    if jsondata:
        drugs = ""
        for rec in jsondata:
            darr = drugs.split(";")
            if rec["drugs"][0]["drugName"] not in darr:
                drugs += rec["drugs"][0]["drugName"]+";"
        return drugs[0:len(drugs)-1]
    else:
        return None


def query_oncokb_cnas_to_csv(cna_annotations: pd.DataFrame, output, i):

    """
    Query OncoKB API to get annotations for copy number alterations (CNAs) and save the results to a CSV file.

    Parameters:
    cna_annotations (DataFrame): DataFrame containing CNA annotations.

    Returns:
    Response: The HTTP response from the OncoKB API.
    """

    api_url = "https://www.oncokb.org/api/v1/annotate/copyNumberAlterations"
    #request_url = api_url + 'copyNameAlterationType='+AlterationType[cna.CNstatus].value+'&hugoSymbol='+hugosymbol+'&tumorType='+tumorType
    header = {'accept':'application/json', 'Content-Type': 'application/json', 'Authorization':'Bearer '+ONCOKB_TOKEN}

    print("Request OncoKB API "+api_url)

    # TODO: No need to query same alteration for every patient and sample, get unique by cnas[i].hugoSymbol cnas[i].alteration

    cnas = cna_annotations.groupby(
        ['hugoSymbol', 'alteration', 'referenceGenome', 'tumorType'])
    uniques = []
    for keys, group in cnas:
        uniques.append(dict(
            {'hugoSymbol': keys[0], 'alteration': keys[1], 'referenceGenome': keys[2], 'tumorType': keys[3]}))

    data = [
        {
            "copyNameAlterationType": f"{str.upper(cna['alteration'])}",
            "referenceGenome": f"{cna['referenceGenome']}",
            "gene": {
                "hugoSymbol": f"{str.upper(cna['hugoSymbol'])}",
            },
            "tumorType": f"{cna['tumorType']}",
        }
        for cna in uniques
    ]

    #header = str(header).replace("'",'"')
    #data = str(data).replace("'",'"')
    print("Querying " +str(len(uniques))+ " CNAs....")

    # Sending a POST request and getting back response as HTTPResponse object.
    #response = urllib3.PoolManager().request("POST", api_url, body=data, headers={'accept':'application/json','Content-Type':'application/json','Authorization':'Bearer '})
    response = httpx.post(api_url, json=data, headers={'Authorization':'Bearer {ONCOKB_TOKEN}'}, timeout=None)


    if (response.status_code == 200):
        treatments = []
        respjson = json.loads(response.text)

        for rjson in respjson:
            hugosymbol = handle_string_field(rjson["query"]["hugoSymbol"])
            alteration = str.upper(handle_string_field(rjson["query"]["alteration"]))

            updatedf = cna_annotations.loc[(cna_annotations['hugoSymbol']==hugosymbol) & (cna_annotations['alteration']==alteration)]
            for indxs, row in updatedf.iterrows():

                cna_annotations.at[indxs,'hugoSymbol'] = handle_string_field(rjson["query"]["hugoSymbol"])
                cna_annotations.at[indxs,'referenceGenome'] = handle_string_field(rjson["query"]["referenceGenome"])
                cna_annotations.at[indxs,'tumorType'] = handle_string_field(rjson["query"]["tumorType"])
                cna_annotations.at[indxs,'consequence'] = handle_string_field(rjson["query"]["consequence"])
                cna_annotations.at[indxs,'oncogenic'] = handle_string_field(rjson["oncogenic"])
                cna_annotations.at[indxs,'mutationEffectDescription'] = handle_string_field(rjson["mutationEffect"]["description"])
                cna_annotations.at[indxs,'gene_role'] = handle_string_field(rjson["mutationEffect"]["knownEffect"])
                cna_annotations.at[indxs,'citationPMids'] = handle_string_field(",".join(rjson["mutationEffect"]["citations"]["pmids"]))
                cna_annotations.at[indxs,'level_of_evidence'] = handle_string_field(rjson["highestSensitiveLevel"]) if handle_string_field(rjson["highestSensitiveLevel"]) else handle_string_field(rjson["highestResistanceLevel"])
                # Hematologic malignancies only
                #updatedf['prognosticSummary'] = handle_string_field(rjson["prognosticSummary"])
                #updatedf['diagnosticSummary'] = handle_string_field(rjson["diagnosticSummary"])
                #updatedf['diagnosticImplications'] = handle_string_field(rjson["diagnosticImplications"])
                #updatedf['prognosticImplications'] = handle_string_field(rjson["prognosticImplications"])
                cna_annotations.at[indxs,'geneSummary'] = handle_string_field(rjson["geneSummary"])
                cna_annotations.at[indxs,'variantSummary'] = handle_string_field(rjson["variantSummary"])
                cna_annotations.at[indxs,'tumorTypeSummary'] = handle_string_field(rjson["tumorTypeSummary"])
                treatments.extend(handle_treatments_oncokb(rjson["treatments"], 'CNA', hugosymbol + ':' + alteration))

            #print("Updated "+str(updatedf.count())+" CNAs")
        #cna_annotations.drop(columns=cna_annotations.columns[0], axis=1, inplace=True)
        header = False if i > 1 else True
        cna_annotations.to_csv(output, mode="a", index=False, header=header, sep="\t", columns=['patient_id', 'sample_id', 'alteration', 'hugoSymbol', 'tumorType', 'consequence', 'oncogenic', 'mutationEffectDescription', 'gene_role', 'citationPMids', 'level_of_evidence', 'cgi_level', 'geneSummary', 'variantSummary', 'tumorTypeSummary'])
        trdf = pd.DataFrame(treatments)
        trdf.to_csv("treatments.csv", mode="a", header=header, index=False, sep="\t")
    else:
        print("Unable to request. Response: ", response.text)

    return response


def query_oncokb_somatic_mutations(snv_annotations: pd.DataFrame, output, i):
    """
    Query OncoKB API to get annotations for somatic mutations and save the results to a CSV file.

    Parameters:
    snv_annotations (DataFrame): DataFrame containing SNV annotations.

    Returns:
    None
    """

    header = {"accept":"application/json", 'Content-Type': 'application/json', "Authorization":'Bearer '+ONCOKB_TOKEN}
    request_url = "https://www.oncokb.org/api/v1/annotate/mutations/byGenomicChange"
    #request_url = "https://www.oncokb.org/api/v1/annotate/mutations/byHGVSg"

    snvs = snv_annotations.groupby(['chromosome', 'position', 'reference_allele', 'sample_allele', 'tumorType', 'referenceGenome'])
    uniques = []
    for keys, group in snvs:
        uniques.append(dict({'chromosome':keys[0], 'position':keys[1], 'reference_allele':keys[2], 'sample_allele':keys[3], 'tumorType':keys[4], 'referenceGenome':keys[5]}))

    data = [
        {
            "id": f"{row['chromosome']+':'+str(row['position'])+':'+row['reference_allele']+':'+row['sample_allele']}",
            "genomicLocation": f"{row['chromosome']+','+str(row['position'])+','+str(int(row['position'])+len(row['sample_allele']))+','+row['reference_allele']+','+row['sample_allele']}",
            "tumorType": f"{row['tumorType']}",
            "referenceGenome": f"{row['referenceGenome']}",
        }
        for row in uniques
    ]

    print("Request OncoKB API "+request_url)
    print("Querying " + str(len(uniques)) + " CNAs....")

    #response = urllib3.PoolManager().request("POST", request_url, body=data, headers={'accept':'application/json','Content-Type':'application/json','Authorization':'Bearer'})
    response = httpx.post(request_url, json=data, headers={'Authorization':'Bearer {ONCOKB_TOKEN}'}, timeout=None)
    print(response.status_code)

    #TODO: check why EGFR chr7,55181426,55181427,A,C  is not found but is found from web api (and also from CGI)
    if (response.status_code == 200):
        treatments = []

        respjson = json.loads(response.text)
        for rjson in respjson:

            id = str(rjson["query"]["id"])
            idsplit = id.split(":")
            chromosome = str(idsplit[0])
            position = int(idsplit[1])
            reference_allele = str(idsplit[2])
            sample_allele = str(idsplit[3])
            updatedf = snv_annotations.loc[(snv_annotations['chromosome']==chromosome) & (snv_annotations['position']==position) & (snv_annotations['reference_allele']==reference_allele) & (snv_annotations['sample_allele']==sample_allele)]

            for indxs, row in updatedf.iterrows():
                alteration = snv_annotations.at[indxs,'hugoSymbol']+":"+chromosome+":"+str(position)+":"+reference_allele+":"+sample_allele
                snv_annotations.at[indxs, 'alteration'] = alteration
                snv_annotations.at[indxs, 'referenceGenome'] = handle_string_field(rjson["query"]["referenceGenome"])
                snv_annotations.at[indxs,'tumorType'] = handle_string_field(rjson["query"]["tumorType"])
                snv_annotations.at[indxs,'consequence'] = handle_string_field(rjson["query"]["consequence"])
                snv_annotations.at[indxs,'oncogenic'] = handle_string_field(rjson["oncogenic"])
                snv_annotations.at[indxs,'mutationEffectDescription'] = handle_string_field(rjson["mutationEffect"]["description"])
                snv_annotations.at[indxs,'gene_role'] = handle_string_field(rjson["mutationEffect"]["knownEffect"])
                snv_annotations.at[indxs,'citationPMids'] = handle_string_field(",".join(rjson["mutationEffect"]["citations"]["pmids"]))
                snv_annotations.at[indxs,'level_of_evidence'] = handle_string_field(rjson["highestSensitiveLevel"]) if handle_string_field(rjson["highestSensitiveLevel"]) else handle_string_field(rjson["highestResistanceLevel"])
                snv_annotations.at[indxs,'geneSummary'] = handle_string_field(rjson["geneSummary"])
                snv_annotations.at[indxs,'variantSummary'] = handle_string_field(rjson["variantSummary"])
                snv_annotations.at[indxs,'tumorTypeSummary'] = handle_string_field(rjson["tumorTypeSummary"])
                treatments.extend(handle_treatments_oncokb(rjson["treatments"], 'SNV', alteration))

        print(snv_annotations)
        header = False if i > 1 else True
        snv_annotations.to_csv(output, mode="a", header=header, index=False, sep="\t", columns=['patient_id', 'sample_id', 'alteration', 'hugoSymbol', 'tumorType', 'consequence', 'oncogenic', 'mutationEffectDescription', 'gene_role', 'citationPMids', 'level_of_evidence', 'cgi_level', 'geneSummary', 'variantSummary', 'tumorTypeSummary'])
        trdf = pd.DataFrame(treatments)
        trdf.to_csv("treatments.csv", header=header, mode="a", index=False, sep="\t")
        #print("Updated " + str(len(snvdf)) + " CNAs")
    else:
        print("[ERROR] Unable to request. Response: ", print(response.text))
        exit()