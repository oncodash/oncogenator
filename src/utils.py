import pandas as pd
import psutil
from pandarallel import pandarallel
from typing import Callable, Optional, Tuple
import json
import httpx


def handle_boolean_field(value, field=None, default=False):
    if str(value).lower() in ["yes", "t", "true"]:
        return True
    elif str(value).lower() in ["no", "f", "false"]:
        return False
    else:
        return None if pd.isna(value) else value


def handle_float_field(value):
    return None if pd.isna(value) or str(value) == "." or str(value).upper() == "NA"else float(value.replace(",", "."))


def handle_cn_type_field(value):
    if str(value).lower() in ["amp", "amplification"]:
        return "AMPLIFICATION"
    elif str(value).lower() in ["del", "deletion"]:
        return "DELETION"
    else:
        return None if pd.isna(value) else "UNKNOWN"


def handle_string_field(value):
    return None if pd.isna(value) or str(value) == "." or str(value).upper() == "NA" else value


def handle_cn_field(value):
    return None if pd.isna(value) or str(value) == "." or str(value).upper() == "NA"else value


def handle_int_field(value):
    return None if pd.isna(value) or str(value) == "." or str(value).upper() == "NA" else int(value)


def handle_date_field(value):
    return None if pd.isna(value) or str(value) == "." else value


def handle_decimal_field(value):
    return None if pd.isna(value) or str(value) == "." or str(value).upper() == "NA" else float(value)


def get_consensus_pathogenecity_prediction(
    clinvar_pathogenecity=None,
    consensus_pathogenecity_source=None,
    gene_role=None,
    consequence=None,
    oncokb_oncogenic=None,
    cgi_oncogenic=None,
    am_score=None,
    polyphen_score=None,
    sift_score=None,
):
    clinvar_sig = str(clinvar_pathogenecity).strip().lower() if clinvar_pathogenecity is not None else None
    if clinvar_sig in ["pathogenic", "likely_pathogenic"]:
        return True, "ClinVar"
    if clinvar_sig in ["benign", "likely_benign"]:
        return False, "ClinVar"

    if consensus_pathogenecity_source != "ClinVar":
        if gene_role == "Loss-of-function" and consequence in ["truncating", "splicing"]: # Likely Loss-of-function
            return True, "truncation_rule"
        
        oncokb_sig = str(oncokb_oncogenic).strip().lower() if oncokb_oncogenic is not None else None
        if oncokb_sig in ["oncogenic", "likely oncogenic", "resistance"]:
            return True, "OncoKB"
        if oncokb_sig in ["likely neutral"]:
            return False, "OncoKB"

    if consensus_pathogenecity_source != "ClinVar" and consensus_pathogenecity_source != "OncoKB" and consensus_pathogenecity_source != "truncation_rule":
            cgi_sig = str(cgi_oncogenic).strip().lower() if cgi_oncogenic is not None else None
            if cgi_sig in ["oncogenic","oncogenic (predicted)","likely oncogenic"]:
                return True, "CGI"
            if cgi_sig in ["non-oncogenic"]:
                return False, "CGI"
    if consensus_pathogenecity_source != "ClinVar" and consensus_pathogenecity_source != "OncoKB" and consensus_pathogenecity_source != "truncation_rule" and consensus_pathogenecity_source != "CGI":  
        if handle_decimal_field(am_score) is not None and handle_decimal_field(am_score) >= 0.5:
            return True, "missense_prediction_algorithms"

        prediction_scores = []
        if handle_decimal_field(polyphen_score) is not None:
            prediction_scores.append(handle_decimal_field(polyphen_score))
        if handle_decimal_field(sift_score) is not None:
            prediction_scores.append(handle_decimal_field(sift_score))
        if len(prediction_scores) > 0:
            consensus_pathogenecity = sum(prediction_scores) >= 0.5
            return consensus_pathogenecity, "missense_prediction_algorithms"

    return None, None


def df_apply(
        df: pd.DataFrame,
        func: Callable,
        col: Optional[str] = None,
        extra_col: Optional[str] = None,
        parallel: bool = True,
        pbar: bool = False,
        cores: Optional[int] = 1,
        **kwargs,
) -> pd.Series:
    """Apply or parallel apply a function to any col of a DataFrame.

    Parameters
    ----------
        df : pd.DataFrame
            Input DataFrame.
        func : Callable
            A callable function.
        col : str, optional,
            The name of the column of the df that is used as the input
            to apply operation.
        extra_col : str, optional
            An extra column that can be used in the apply operation.
        parallel : bool, default=False
            Flag, whether to parallelize the operation with pandarallel
        pbar : bool, default=False
            Show progress bar when executing in parallel mode. Ignored if
            `parallel=False`
        **kwargs:
            Arbitrary keyword args for the `func` callable.

    Returns
    -------
        Output of applied function
    """

    if not parallel:
        if col is None:
            res = df.apply(func, **kwargs)
        else:
            if extra_col is None:
                res = df[col].apply(func, **kwargs)
            else:
                res = df[[col, extra_col]].apply(lambda x: func(*x, **kwargs), axis=1)
    else:
        pandarallel.initialize(nb_workers=cores, verbose=1, progress_bar=pbar)
        if col is None:
            res = df.parallel_apply(func, **kwargs, axis=1)
        else:
            if extra_col is None:
                res = df[col].parallel_apply(func, **kwargs)
            else:
                res = df[[col, extra_col]].parallel_apply(
                    lambda x: func(*x, **kwargs), axis=1
                )

    return res

def gene_id_convert(geneids, target):
    # geneids given as list with whitespace separator, target can be one of the target namespaces in https://biit.cs.ut.ee/gprofiler/convert
    request_url = "https://biit.cs.ut.ee/gprofiler/api/convert/convert/"
    print("Request gProfiler API "+request_url)
    data = '{"organism":"hsapiens", "target":"'+target+'", "query":"'+geneids+'"}'
    headers = {"Content-Type": "application/json"}
    body = json.dumps(data).encode('utf-8')
    response = httpx.post(request_url, json=body, headers=headers, timeout=None)

    print("response.status", response.status)
    print(data)
    print(response.json)
    if (response.status == 200):
        rjson = json.loads(response.data.decode('utf-8'))
        print(dict(rjson['result'][0]).get('converted'))
        return dict(rjson['result'][0]).get('converted')

    else:
        print("[ERROR] Unable to request. Response: ", print(response.data))
        exit()
