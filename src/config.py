"""
Configuration file for oncogenator.
Contains API credentials and settings.
"""

# CGI (Cancer Genome Interpreter) API credentials
CGI_LOGIN = "ilari.maarala@helsinki.fi"
CGI_TOKEN = "f7ab7f3911b1194f4629"

# CGI API settings
CGI_API_URL = "https://www.cancergenomeinterpreter.org/api/v1"
CGI_DEFAULT_CANCER_TYPE = "CANCER" #HGSOC for high serous ovarian, LUNG for lung, etc.
CGI_DEFAULT_REFERENCE = "hg38"

# OncoKB API credentials
ONCOKB_TOKEN = "a35f6d69-8895-4742-8b3d-16bd7d765793"

# OncoKB API settings
ONCOKB_API_URL = "https://www.oncokb.org/api/v1"
ONCOKB_CNA_ENDPOINT = "https://www.oncokb.org/api/v1/annotate/copyNumberAlterations"
ONCOKB_MUTATION_ENDPOINT = "https://www.oncokb.org/api/v1/annotate/mutations/byGenomicChange"
