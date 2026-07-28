from utils import *
import pandas as pd

import json
import httpx
from config import ONCOKB_TOKEN, ONCOKB_CNA_ENDPOINT, ONCOKB_MUTATION_ENDPOINT

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
        #print(row)
        if row["drugs"]:
            #print(row["drugs"])
            for drug in row["drugs"]:
                #print(drug)
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

    api_url = ONCOKB_CNA_ENDPOINT
    #request_url = api_url + 'copyNameAlterationType='+AlterationType[cna.CNstatus].value+'&hugoSymbol='+hugosymbol+'&tumorType='+tumorType
    header = {'Authorization':'Bearer '+ONCOKB_TOKEN}

    print("Request OncoKB API "+api_url)

    # TODO: No need to query same alteration for every patient and sample, get unique by cnas[i].hugoSymbol cnas[i].alteration

    #cnas = cna_annotations.groupby(['hugoSymbol', 'alteration', 'referenceGenome', 'tumorType'])
    #uniques = []
    #for keys, row in cna_annotations.iterrows():
    #    uniques.append(dict(
    #        {"hugoSymbol": row["hugoSymbol"], "alteration": row["alteration"], "referenceGenome": row["referenceGenome"], "tumorType": row["tumorType"]}))

    data = [
        {
            "copyNameAlterationType": f"{str.upper(cna['alteration'])}",
            "referenceGenome": "GRCh38",#f"{cna['referenceGenome']}",
            "gene": {
                "hugoSymbol": f"{str.upper(cna['hugoSymbol'])}"
            }
            #"tumorType": f"{cna['tumorType']}",
        }
        for keys, cna in cna_annotations.iterrows()
    ]

    #with open(output.split(".")[0]+"_cna_payload.json", "w") as payload_file:
    #    json.dump(data, payload_file, indent=2)
    
    #header = str(header).replace("'",'"')
    #data = str(data).replace("'",'"')
    print("Querying " +str(len(cna_annotations))+ " CNAs....")
    
    # Sending a POST request and getting back response as HTTPResponse object.
    #response = urllib3.PoolManager().request("POST", api_url, body=data, headers={'accept':'application/json','Content-Type':'application/json','Authorization':'Bearer '})
    response = httpx.post(api_url, json=data, headers=header, timeout=None)
    
    if (response.status_code == 200):
        treatments = []
        #print(response.text)

        respjson = json.loads(response.text)
        #with open(output.split(".")[0]+"_response.json", "w") as payload_file:
        #    json.dump(respjson, payload_file, indent=2)
        
        for rjson in respjson:
            hugosymbol = handle_string_field(rjson["query"]["hugoSymbol"])
            alteration = str.upper(handle_string_field(rjson["query"]["alteration"]))
            #TODO: Do not update but create new dataf
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
                #cna_annotations.at[indxs,'level_of_evidence'] = handle_string_field(rjson["highestSensitiveLevel"]) if handle_string_field(rjson["highestSensitiveLevel"]) else handle_string_field(rjson["highestResistanceLevel"])

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
        header = True if i == 0 else False
        cna_annotations.to_csv(output, index=False, header=True, sep="\t", columns=['patient_id', 'sample_id', 'alteration', 'hugoSymbol', 'ensembl_id', 'referenceGenome', 'tumorType', 'consequence', 'oncogenic', 'mutationEffectDescription', 'gene_role', 'citationPMids', 'geneSummary', 'variantSummary', 'tumorTypeSummary'])
        trdf = pd.DataFrame(treatments)
        trdf.to_csv(output.split(".")[0]+"_treatments.csv", mode="a", header=True,index=False, sep="\t")
    else:
        print("Unable to request. Response: ", response.text)

    return response


def query_oncokb_somatic_mutations(somatic_mutation_annotations: pd.DataFrame, output, i):
    """
    Query OncoKB API to get annotations for somatic mutations and save the results to a CSV file.

    Parameters:
    somatic_mutation_annotations (DataFrame): DataFrame containing somatic_mutation annotations.

    Returns:
    None
    """


    header = {'Authorization':'Bearer '+ONCOKB_TOKEN}
    request_url = ONCOKB_MUTATION_ENDPOINT
    #request_url = "https://www.oncokb.org/api/v1/annotate/mutations/byHGVSg"

    #somatic_mutations = somatic_mutation_annotations.groupby(['chromosome', 'position', 'ensembl_id', 'reference_allele', 'sample_allele', 'tumorType', 'referenceGenome'])
    #uniques = []
    #for keys, group in somatic_mutations:
    #    uniques.append(dict({'chromosome':keys[0], 'position':keys[1], 'reference_allele':keys[2], 'sample_allele':keys[3], 'tumorType':keys[4], 'referenceGenome':keys[5]}))

    data = [
        {
            "id": f"{row['chromosome']+':'+str(row['position'])+':'+row['reference_allele']+':'+row['sample_allele']}",
            "genomicLocation": f"{row['chromosome']+','+str(row['position'])+','+str(int(row['position'])+(len(row['sample_allele'])-len(row['reference_allele'])))+','+row['reference_allele']+','+row['sample_allele']}",
            #"tumorType": f"{row['tumorType']}",
            "referenceGenome": f"{row['referenceGenome']}"
        }
        for keys, row in somatic_mutation_annotations.iterrows()
    ]

    print("Request OncoKB API "+request_url)
    print("Querying " + str(len(somatic_mutation_annotations)) + " somatic_mutations....")
    print(str(data))
    #response = urllib3.PoolManager().request("POST", request_url, body=data, headers={'accept':'application/json','Content-Type':'application/json','Authorization':'Bearer'})
    response = httpx.post(request_url, json=data, headers=header, timeout=None)
    #print(response.status_code)

    #TODO: check why EGFR chr7,55181426,55181427,A,C  is not found but is found from web api (and also from CGI)
    if (response.status_code == 200):
        treatments = []

        respjson = json.loads(response.text)
        #with open(output.split(".")[0]+"_response.json", "w") as payload_file:
        #    json.dump(respjson, payload_file, indent=2)
        #print(respjson)
        for rjson in respjson:

            id = str(rjson["query"]["id"])
            idsplit = id.split(":")
            chromosome = str(idsplit[0])
            position = int(idsplit[1])
            reference_allele = str(idsplit[2])
            sample_allele = str(idsplit[3])
            updatedf = somatic_mutation_annotations.loc[(somatic_mutation_annotations['chromosome']==chromosome) & (somatic_mutation_annotations['position']==position) & (somatic_mutation_annotations['reference_allele']==reference_allele) & (somatic_mutation_annotations['sample_allele']==sample_allele)]

            for indxs, row in updatedf.iterrows():
                alteration = somatic_mutation_annotations.at[indxs,'hugoSymbol']+":"+chromosome+":"+str(position)+":"+reference_allele+":"+sample_allele
                somatic_mutation_annotations.at[indxs, 'alteration'] = alteration
                somatic_mutation_annotations.at[indxs, 'referenceGenome'] = handle_string_field(rjson["query"]["referenceGenome"])
                somatic_mutation_annotations.at[indxs,'tumorType'] = handle_string_field(rjson["query"]["tumorType"])
                somatic_mutation_annotations.at[indxs,'oncokb_consequence'] = handle_string_field(rjson["query"]["consequence"])
                somatic_mutation_annotations.at[indxs,'oncokb_oncogenic'] = handle_string_field(rjson["oncogenic"])
                somatic_mutation_annotations.at[indxs,'mutationEffectDescription'] = handle_string_field(rjson["mutationEffect"]["description"])
                somatic_mutation_annotations.at[indxs,'gene_role'] = handle_string_field(rjson["mutationEffect"]["knownEffect"])
                somatic_mutation_annotations.at[indxs,'citationPMids'] = handle_string_field(",".join(rjson["mutationEffect"]["citations"]["pmids"]))
                #somatic_mutation_annotations.at[indxs,'level_of_evidence'] = handle_string_field(rjson["highestSensitiveLevel"]) if handle_string_field(rjson["highestSensitiveLevel"]) else handle_string_field(rjson["highestResistanceLevel"])
                somatic_mutation_annotations.at[indxs,'geneSummary'] = handle_string_field(rjson["geneSummary"])
                somatic_mutation_annotations.at[indxs,'variantSummary'] = handle_string_field(rjson["variantSummary"])
                somatic_mutation_annotations.at[indxs,'tumorTypeSummary'] = handle_string_field(rjson["tumorTypeSummary"])
                print(row)
                if somatic_mutation_annotations.at[indxs, 'consensus_pathogenecity_source'] is not None and somatic_mutation_annotations.at[indxs, 'consensus_pathogenecity_source'] == "ClinVar":
                    consensus_pathogenecity = somatic_mutation_annotations.at[indxs, 'consensus_pathogenecity']
                    consensus_prediction_source = somatic_mutation_annotations.at[indxs, 'consensus_pathogenecity_source'] 
                else:   
                    consensus_pathogenecity, consensus_prediction_source = get_consensus_pathogenecity_prediction(
                        clinvar_pathogenecity=handle_string_field(row['clinvar_pathogenecity']) if 'clinvar_pathogenecity' in row else None,
                        consensus_pathogenecity_source=handle_string_field(row['consensus_pathogenecity_source']) if 'consensus_pathogenecity_source' in row else None,
                        gene_role=handle_string_field(rjson["mutationEffect"]["knownEffect"]),
                        oncokb_oncogenic=handle_string_field(rjson["oncogenic"])
                    )
                somatic_mutation_annotations.at[indxs, 'consensus_pathogenecity'] = consensus_pathogenecity
                somatic_mutation_annotations.at[indxs, 'consensus_pathogenecity_source'] = consensus_prediction_source
                # FIXME: for some reason treatments are not being handled properly for SNVs
            
                treatments.extend(handle_treatments_oncokb(rjson["treatments"], 'SNV', alteration))

        #print(somatic_mutation_annotations.columns)
        header = True if i == 0 else False
        print(somatic_mutation_annotations.columns)
        somatic_mutation_annotations.to_csv(output, header=True, index=False, sep="\t", columns=['patient_id', 'sample_id', 'alteration', 'hugoSymbol', 'ensembl_id', 'tumorType', 'consequence', 'annovar_consequence', 'oncokb_consequence', 'oncokb_oncogenic', 'mutationEffectDescription', 'gene_role', 'citationPMids', 'geneSummary', 'variantSummary', 'tumorTypeSummary',  'expressed', 'refCount', 'altCount', 'consensus_pathogenecity', 'consensus_pathogenecity_source'])
        trdf = pd.DataFrame(treatments)
        trdf.to_csv(output.split(".")[0]+"_treatments.csv", header=True, index=False, sep="\t")
        #print("Updated " + str(len(somatic_mutationdf)) + " CNAs")
    else:
        print("[ERROR] Unable to request. Response: ", print(response.text))
        exit()
