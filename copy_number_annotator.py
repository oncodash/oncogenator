from utils import *

class CopyNumberAnnotator:

    def __init__(self, refgenome="GRCh38", tumortype="HGSOC", ascats=None, ploidy_coeff=2.5):
        """
            Initialize the CopyNumberAnnotator class.

            Parameters:
            refgenome (str): Reference genome version, default is "GRCh38".
            tumortype (str): Type of tumor, default is "HGSOC".
            ascats (DataFrame): DataFrame containing ASCAT results.
            ploidy_coeff (float): Coefficient for ploidy filtering, default is 2.5.
        """
        self.refgenome = refgenome
        self.tumortype = tumortype
        self.ascats = ascats
        self.ploidy_coeff = ploidy_coeff

    def filter_cnas_by_ploidy(self, row):
        """
                Filter copy number alterations (CNAs) based on ploidy.

                Parameters:
                row (Series): A row from a DataFrame containing CNA data.

                Returns:
                Series: A Series containing filtered CNA data if conditions are met, otherwise None.
        """

        ploidy = self.ascats.loc[self.ascats['sample'] == row['sample']]['ploidy']
        nminor = handle_int_field(row['nMinor'])
        nmajor = handle_int_field(row['nMajor'])
        if nminor and nmajor:
            cn = int(nminor) + int(nmajor)
            ploidy = ploidy.iloc[0]
            if ploidy > 0 and float(ploidy) > 0:
                if cn < 1 or cn > self.ploidy_coeff * float(ploidy):
                    return pd.Series({
                        'patient_id': row["sample"].split("_")[0],
                        # Use cohort code here, map to pid later to reduce queries sample name includes cohort code which is mapped to patient id
                        'sample_id': handle_string_field(row["sample"]),
                        'referenceGenome': self.refgenome,
                        'hugoSymbol': handle_string_field(row["Gene"]),
                        'alteration': handle_cn_type_field(row["CNstatus"]),
                        'tumorType': handle_string_field(self.tumortype),
                        'nMajor': handle_int_field(row["nMajor"]),
                        'nMinor': handle_int_field(row["nMinor"]),
                        'lohstatus': handle_string_field(row["LOHstatus"]),
                        'ploidy': handle_decimal_field(float(ploidy)),
                        'start': handle_string_field(row["start"]),
                        'end': handle_string_field(row["end"]),
                        'strand': handle_string_field(row["strand"]),
                        'band': handle_string_field(row["band"]),
                        'nProbesCr': handle_string_field(row["nProbesCr"]),
                        'nProbesAf': handle_string_field(row["nProbesAf"]),
                        'logR': handle_decimal_field(row["logR"]),
                        'baf': handle_decimal_field(row["baf"]),
                        'nAraw': handle_decimal_field(row["nAraw"]),
                        'nBraw': handle_decimal_field(row["nBraw"]),
                        'purifiedLogR': handle_decimal_field(row["purifiedLogR"]),
                        'purifiedBaf': handle_decimal_field(row["purifiedBaf"]),
                        'purifiedLoh': handle_decimal_field(row["purifiedLoh"]),
                        'minPurifiedLogR': handle_decimal_field(row["minPurifiedLogR"]),
                        'maxPurifiedLogR': handle_decimal_field(row["maxPurifiedLogR"]),
                        'breaksInGene': handle_string_field(row["breaksInGene"]),

                    })
        return None
