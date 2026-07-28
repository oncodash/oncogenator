# TODO get ploidy from Purple output /mnt/storageBig8/web-pub/projects/HERCULES/WGS/StructuralVariations_all/v2.0/purity_ploidy_estimates.tsv
from utils import *

class CopyNumberAnnotator:

    def __init__(self, refgenome="GRCh38", tumortype="HGSOC", ascats=None, ploidy_coeff=2.5, sample_info=None):
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
        self.sample_info = sample_info

    @staticmethod
    def _to_bool(series):
        return series.astype(str).str.strip().str.lower().map({
            "true": True,
            "false": False,
    })

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
        passes_sample_info_filters = True
        if self.sample_info is not None:
            sinfo = self.sample_info.loc[self.sample_info['sample'] == row['sample']]
            usable = self._to_bool(sinfo["usable"]) == True
            contam_filter = self._to_bool(sinfo["contamFilter"]) == False
            duplicate = self._to_bool(sinfo["duplicate"]) == False
            cell_line = self._to_bool(sinfo["cellLine"]) == False
            passes_sample_info_filters = (usable & contam_filter & duplicate & cell_line).all()

        if not passes_sample_info_filters:
            print(f"Sample {row['sample']} failed sample info filters. Skipping...")
            return None

        if not (nminor and nmajor):
            print(f"DEBUG: Sample {row['sample']} gene {row.get('Gene', '?')} skipped — nMinor={nminor!r} or nMajor={nmajor!r} is missing/zero.")
            return None

        if nminor and nmajor:
            cn = int(nminor) + int(nmajor)
            ploidy_val = ploidy.iloc[0] if len(ploidy) > 0 else None
            if ploidy_val is None or float(ploidy_val) <= 0:
                print(f"DEBUG: Sample {row['sample']} gene {row.get('Gene', '?')} skipped — ploidy={ploidy_val!r} is missing or non-positive. nminor={nminor!r}, nmajor={nmajor!r}, cn={cn}. Ploidy coeff={self.ploidy_coeff * float(ploidy_val) if ploidy_val is not None else 'N/A'}.")
                return None
            ploidy = ploidy_val
            if ploidy > 0 and float(ploidy) > 0:
                if cn < 1 or cn > self.ploidy_coeff * float(ploidy):
                    return pd.Series({
                        'patient_id': row["sample"].split("_")[0],
                        # Use cohort code here, map to pid later to reduce queries sample name includes cohort code which is mapped to patient id
                        'sample_id': handle_string_field(row["sample"]),
                        'referenceGenome': self.refgenome,
                        'ensembl_id': handle_string_field(row["ID"]),
                        'hugoSymbol': handle_string_field(row["Gene"]),
                        'alteration': handle_cn_type_field(row["CNstatus"]),
                        'tumorType': handle_string_field(self.tumortype),
                        'nMajor': handle_int_field(row["nMajor"]),
                        'nMinor': handle_int_field(row["nMinor"]),
                        'lohstatus': handle_string_field(row["LOHstatus"]),
                        'ploidy': handle_decimal_field(float(ploidy)),
                        'chromosome': handle_string_field(row["chr"]),
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
            else:
                print(f"DEBUG: Sample {row['sample']} gene {row.get('Gene', '?')} skipped — cn={cn} does not meet threshold (must be <1 or >{self.ploidy_coeff} * ploidy={float(ploidy):.2f}).")
        return None
