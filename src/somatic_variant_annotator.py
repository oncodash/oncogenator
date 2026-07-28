# TODO: get purity from Purple output /mnt/storageBig8/web-pub/projects/HERCULES/WGS/StructuralVariations_all/v2.0/purity_ploidy_estimates.tsv
# TODO integrate expression data fields refCount and altCount from /mnt/storageBig8/work/joikkone/cohort/mutations_in_RNA/result_SNV_calling_RNA/D327_p2Asc1_RNA1-ASEcall-bash/out1.table by matching patientID chromosome(contig) position referenceAllele and sampleAllele

from utils import *
from scipy import stats
#import dask.dataframe as dd
import pandas as pd
import re

# Default thresholds
expression_threshold = 5
homogeneity_threshold=0.05
ada_score_threshold=0.95
rf_score_threshold=0.95
tumor_type="CANCER" # or eg. HGSOC

class SomaticVariantAnnotator:
    def __init__(self, refgenome="GRCh38", tumortype=tumor_type, cnas=None, ascats=None, samples=None, rna_path=None, sample_info=None, homogeneity_threshold=homogeneity_threshold, ada_score_threshold=ada_score_threshold, rf_score_threshold=rf_score_threshold):
        """
                Initialize the SomaticVariantAnnotator class.

                Parameters:
                refgenome (str): Reference genome version, default is "GRCh38".
                tumortype (str): Type of tumor, default is "CANCER".
                cnas (DataFrame): DataFrame containing CNA data.
                ascats (DataFrame): DataFrame containing ASCAT results.
                samples (list): List of sample identifiers.
                homogeneity_threshold (float): Threshold for homogeneity filtering, default is 0.05.
                ada_score_threshold (float): Threshold for ADA score filtering, default is 0.95.
                rf_score_threshold (float): Threshold for RF score filtering, default is 0.95.
        """
        self.homogeneity_threshold = homogeneity_threshold
        self.rf_score_threshold = rf_score_threshold
        self.ada_score_threshold = ada_score_threshold
        self.samples = samples
        self.ascats = ascats
        self.cnas = cnas
        self.refgenome = refgenome
        self.tumortype = tumortype
        self.rna_path = rna_path
        self.sample_info = sample_info
    @staticmethod
    def get_variant_assoc_cnas(cnas, sid, gene):
        """
                Get associated CNAs for a given variant.

                Parameters:
                cnas (DataFrame): DataFrame containing CNA data.
                sid (str): Sample identifier.
                gene (str): Gene symbol.

                Returns:
                Series: A Series containing CNA data for the given sample and gene, or an empty list if no data is found.
        """
        cnar = cnas.loc[(cnas['sample'] == sid) & (cnas['Gene'] == gene)]
        return cnar.iloc[0] if not cnar.empty else []

    @staticmethod
    def expectedAF(N_t, CN_t, TF):
        """
                Calculate the expected allele frequency.

                Parameters:
                N_t (int): Tumor copy number.
                CN_t (int): Total copy number.
                TF (float): Tumor fraction.

                Returns:
                float: Expected allele frequency.
        """
        return (N_t * TF) / (CN_t * TF + 2 * (1 - TF))

    @staticmethod
    def parse_isoforms(aaChangeRefGene):
        """
                Parse isoforms from the AAChange.refGene field.

                Parameters:
                aaChangeRefGene (str): AAChange.refGene field value.

                Returns:
                list: List of unique isoforms.
        """
        records = aaChangeRefGene.split(",")
        isoforms = []
        for rec in records:
            fields = rec.split(":")
            if len(fields) > 1:
                gene = fields[0]
                protein = fields[len(fields) - 1].split(".")[1]
                isoform = gene + ":" + protein
                isoforms.append(isoform)
        return list(dict.fromkeys(isoforms))

    def create_somatic_mutation_annotation(self, row, pid, sample_id, gene, alteration, consequence, nMinor, nMajor, lohstatus, expHomAF, expHomCI_lo, expHomCI_hi, expHom_pbinom_lower, homogenous, ad0, ad1, ensemble_id="", depth=0, polyphen_score="", polyphen_category="", AM_class="", AM_score="", sift_score="", sift_category="", classification="", consensus_pathogenecity=None, consensus_pathogenecity_source=None, expressed=False, refCount=0, altCount=0):
        """
                Create an somatic_mutation annotation.

                Parameters:
                row (Series): A row from a DataFrame containing somatic_mutation data.
                pid (str): Patient identifier.
                sample_id (str): Sample identifier.
                gene (str): Gene symbol.
                alteration (str): Alteration description.
                consequence (str): Consequence of the alteration.
                nMinor (int): Number of minor alleles.
                nMajor (int): Number of major alleles.
                lohstatus (str): LOH status.
                expHomAF (float): Expected homogenous allele frequency.
                expHomCI_lo (float): Lower bound of the expected homogenous confidence interval.
                expHomCI_hi (float): Upper bound of the expected homogenous confidence interval.
                expHom_pbinom_lower (float): Lower bound of the binomial probability.
                homogenous (bool): Homogeneity status.

                ad0 (int): Allele depth for reference allele.
                ad1 (int): Allele depth for alternate allele.
                depth (int): Total depth.
                AM_class (str): AMIS category.
                AM_score (float): AMIS score.
                polyphen_score (float): PolyPhen score.
                polyphen_category (str): PolyPhen category.
                sift_score (float): SIFT score.
                sift_category (str): SIFT category.

                Returns:
                Series: A Series containing the somatic_mutation annotation.
        """
        # v4.10 fields
        # patient CHROM   POS     REF     ALT     ID      FILTER  cytoBand        Func.MANE       Gene.MANE       GeneDetail.MANE ExonicFunc.MANE AAChange.MANE   Func.refGene    
        # Gene.refGene    GeneDetail.refGeneExonicFunc.refGene       AAChange.refGene        genomicSuperDups        dbscsomatic_mutation_ADA_SCORE       dbscsomatic_mutation_RF_SCORE        COSMIC_ID       
        # COSMIC_OCCURRENCE       COSMIC_TOTAL_OCC        COSMIC_CONF_SOMA  CLNSIG   CLNSIGCONF      CLNDN   CLNREVSTAT      CLNALLELEID     CLNDISDB        ONC     ONCCONF ONCDN  
        #  ONCDISDB        ONCREVSTAT      SCI     SCIDN   SCIDISDB        SCIREVSTAT      AM_variant      AM_score   AM_class        PolyPhenVal     PolyPhenCat     SIFTval SIFTcat 
        # Interpro_transcript     Interpro_domain regulomeDB      CADD_raw        CADD_phred      1000G_ALL       1000G_EUR       gnomAD_joint_ALL   gnomAD_joint_NFE        
        # gnomAD_joint_FIN        gnomAD_joint_max        Truncal readCounts      VAFs    samples

        return pd.Series({
            'patient_id': pid,
            'sample_id': sample_id,
            'ref_id': handle_string_field(row["ID"]),
            'chromosome': handle_string_field(row["CHROM"]),
            'position': handle_int_field(row["POS"]),
            'reference_allele': handle_string_field(row["REF"]),
            'sample_allele': handle_string_field(row["ALT"]),
            'referenceGenome': "GRCh38",
            'ensembl_id': ensemble_id,
            'hugoSymbol': gene,
            'alteration': alteration,
            'tumorType': tumor_type,
            'consequence': handle_string_field(classification).lower() if handle_string_field(classification) else "",
            'annovar_consequence': consequence,
            'cytoBand': handle_string_field(row['cytoBand']),
            'exonicFuncMane': handle_string_field(row["ExonicFunc.MANE"]),
            'funcMane': handle_string_field(row['Func.MANE']),
            'AAChangeMANE': handle_string_field(row['AAChange.MANE']),
            'funcRefGene': handle_string_field(row['Func.refGene']),
            'exonicFuncRefGene': handle_string_field(row['ExonicFunc.refGene']),
            'AAChangerefGene': handle_string_field(row['AAChange.refGene']),
            'nMinor': nMinor,
            'nMajor': nMajor,
            'ad0': ad0,
            'ad1': ad1,
            'af': expHomAF,
            'depth': depth,
            'lohstatus': lohstatus,
            'hom_lo': "{:.9f}".format(expHomCI_lo),
            'hom_hi': "{:.9f}".format(expHomCI_hi),
            'hom_pbinom_lo': "{:.9f}".format(expHom_pbinom_lower),
            'homogenous': homogenous,
            'cadd_score': handle_decimal_field(row["CADD_phred"]),
            'ada_score': handle_decimal_field(row["dbscSNV_ADA_SCORE"]),
            'rf_score': handle_decimal_field(row["dbscSNV_RF_SCORE"]),
            'sift_category': sift_category,
            'sift_score': sift_score,
            'polyphen_category': polyphen_category,
            'polyphen_score': polyphen_score,
            'AM_class': AM_class,
            'AM_score': AM_score,
            'cosmic_id': handle_string_field(row["COSMIC_ID"]),
            'clinvar_id': handle_string_field(row["CLNALLELEID"]),
            'clinvar_sig': handle_string_field(row["CLNSIG"]), 
            'clinvar_assoc': handle_string_field(row["CLNDN"]), 
            'clinvar_status': handle_string_field(row["CLNREVSTAT"]),
            #'classification': handle_string_field(classification),
            'consensus_pathogenecity': consensus_pathogenecity,
            'consensus_pathogenecity_source': consensus_pathogenecity_source,
            'refCount': refCount,
            'altCount': altCount,
            'expressed': expressed
        })

    @staticmethod
    def _to_bool(series):
        return series.astype(str).str.strip().str.lower().map({
            "true": True,
            "false": False,
    })

    def filter_and_classify_somatic_mutations(self, row):
        """
                Filter and classify somatic_mutations based on various criteria.

                Parameters:
                row (Series): A row from a DataFrame containing somatic_mutation data.

                Returns:
                list: A list of Series containing somatic_mutation annotations.
        """
        sv_class = None
        somatic_mutation_annotations = []
		
        print(row)
        
        exonicFuncMane = handle_string_field(row["ExonicFunc.MANE"])
        funcMane = handle_string_field(row["Func.MANE"])
        funcRefgene = handle_string_field(row["Func.refGene"])
        ada_score = handle_decimal_field(row["dbscSNV_ADA_SCORE"])
        rf_score = handle_decimal_field(row["dbscSNV_RF_SCORE"])

        for sample_id in self.samples:

            # Check sample info filters if provided
            passes_sample_info_filters = True
            if self.sample_info is not None:
                sinfo = self.sample_info.loc[self.sample_info['sample'] == sample_id]
                usable = self._to_bool(sinfo["usable"]) == True
                contam_filter = self._to_bool(sinfo["contamFilter"]) == False
                duplicate = self._to_bool(sinfo["duplicate"]) == False
                cell_line = self._to_bool(sinfo["cellLine"]) == False
                passes_sample_info_filters = (usable & contam_filter & duplicate & cell_line).all()

            if not passes_sample_info_filters:
                print(f"Sample {sample_id} failed sample info filters. Skipping...")
                continue
            
            sample_name_split = sample_id.split("_")
            pid = sample_name_split[0]
            siteid = sample_name_split[1]

            rna_expression = None
            # Add rna expression data if available
            if self.rna_path:
                for rna_num in range(1, 6):  # Try RNA1 through RNA5
                    try:
                        rna_path = self.rna_path + f"/{pid}_{siteid}_RNA{rna_num}-ASEcall-bash/out1.table"
                        rna_expression = pd.read_csv(rna_path, sep="\t")
                        if len(rna_expression) > 0:
                            print(f"Found RNA expression data for {sample_id} in RNA{rna_num}: {len(rna_expression)} records")
                            break
                    except Exception as e:
                        print(f"No RNA expression data found {rna_path}. Error: {e}")
                        pass
                
                if rna_expression is None or len(rna_expression) == 0:
                    print(f"No RNA expression data found for sample {sample_id}")
                    

            purities = self.ascats.loc[self.ascats['sample'] == sample_id]['purity']
            purity = purities.iloc[0] if len(purities) > 0 else 0.0

            try:
                depth = int(row[str(sample_id)+".DP"])
                ad0 = int(row[str(sample_id)+".AD"].split(',')[0])
                ad1 = int(row[str(sample_id)+".AD"].split(',')[1])
            except Exception as e:
                print(f"Sample not found from somatic variants. Error processing sample {sample_id}: {e}")
                continue

            # Get genes associated with the variant from both MANE and refGene annotations    
            geneMANE = re.split(r'[;,\s]+', handle_string_field(row["Gene.MANE"]))
            genes = set(geneMANE)
            geneRefGene = re.split(r'[;,\s]+', handle_string_field(row["Gene.refGene"]))
            for g in geneRefGene:
                genes.add(g)
            i = 0
            for gene in genes:
                AM_score = row['AM_score'].split(',')[i] if ',' in str(row['AM_score']) and len(row['AM_score'].split(',')) > i else row['AM_score']
                AM_class = row['AM_class'].split(',')[i] if ',' in str(row['AM_class']) and len(row['AM_class'].split(',')) > i else row['AM_class']
                pathogenecity = handle_string_field(row["CLNSIG"].split(',')[i] if ',' in str(row["CLNSIG"]) and len(row["CLNSIG"].split(',')) > i else row["CLNSIG"])
                polyphen_score = handle_decimal_field(row["PolyPhenVal"].split(',')[i] if ',' in str(row["PolyPhenVal"]) and len(row["PolyPhenVal"].split(',')) > i else row["PolyPhenVal"])
                sift_score = handle_decimal_field(row["SIFTval"].split(',')[i] if ',' in str(row["SIFTval"]) and len(row["SIFTval"].split(',')) > i else row["SIFTval"])
                sift_category = handle_string_field(row["SIFTcat"].split(',')[i] if ',' in str(row["SIFTcat"]) and len(row["SIFTcat"].split(',')) > i else row["SIFTcat"])
                poylphen_category = handle_string_field(row["PolyPhenCat"].split(',')[i] if ',' in str(row["PolyPhenCat"]) and len(row["PolyPhenCat"].split(',')) > i else row["PolyPhenCat"])
                consensus_pathogenecity, consensus_prediction_source = get_consensus_pathogenecity_prediction(
                    clinvar_pathogenecity=pathogenecity,
                    am_score=AM_score,
                    polyphen_score=polyphen_score,
                    sift_score=sift_score,
                )
                # Get CNA data for the sample and gene
                vcnas = self.get_variant_assoc_cnas(self.cnas, sample_id, gene)
                ensemble_id = handle_string_field(vcnas["ID"]) if len(vcnas) > 0 else None
                nMajor = handle_cn_field(vcnas['nMajor']) if len(vcnas) > 0 else None
                nMinor = handle_cn_field(vcnas['nMinor']) if len(vcnas) > 0 else None
                lohstatus = vcnas['LOHstatus'] if len(vcnas) > 0 else None

                expHomAF = 0.0
                expHomCI_lo = 0.0
                expHomCI_hi = 0.0
                expHom_pbinom_lower = 0.0
                homogenous = None

                # Calculate homogeneity estimate if CNA data is available
                # TODO: DOn't filter by homogeneity just show in frontend
                if nMajor and nMinor:
                    cn = float(nMinor) + float(nMajor)
                    expHomAF = float(self.expectedAF(cn, cn, purity))
                    expHomCI_lo = float(stats.binom.ppf(0.025, depth, expHomAF))
                    expHomCI_hi = float(stats.binom.ppf(0.975, depth, expHomAF))
                    expHomCI_cover = expHomCI_lo <= ad1
                    expHom_pbinom_lower = float(stats.binom.cdf(ad1, depth, expHomAF))
                    homogenous = expHom_pbinom_lower > self.homogeneity_threshold

                # Classify variant based on gene function and homogeneity
                if exonicFuncMane == "nonsynonymous_SNV":
                    sv_class = "Missense"
                if exonicFuncMane in ["frameshift_insertion", "frameshift_deletion", "stopgain"]:
                    sv_class = "Truncating"
                if exonicFuncMane in ["nonframeshift_deletion", "nonframeshift_substitution", "nonframeshift_insertion"]:
                    sv_class = "Other"
                consequence = exonicFuncMane
                if not sv_class:
                    if funcMane in ["splicing", "splicesite", "intron", "intronic"] or funcRefgene in ["splicing", "splicesite", "intron", "intronic"]:
                        consequence = funcMane
                        if (ada_score and float(ada_score) > self.ada_score_threshold) or (rf_score and float(rf_score) > self.rf_score_threshold):
                            sv_class = "Splicing"
                        else:
                            continue

                alteration = f"{gene}:{row['CHROM']}:{row['POS']}:{row['REF']}:{row['ALT']}"
                expressed = ""
                refCount = ""
                altCount = ""
                if rna_expression is not None:
                    refCount = rna_expression.loc[(rna_expression['contig'] == row['CHROM']) & (rna_expression['position'] == row['POS']) & (rna_expression['refAllele'] == row['REF']) & (rna_expression['altAllele'] == row['ALT']), 'refCount'].values[0] if len(rna_expression.loc[(rna_expression['contig'] == row['CHROM']) & (rna_expression['position'] == row['POS']) & (rna_expression['refAllele'] == row['REF']) & (rna_expression['altAllele'] == row['ALT'])]) > 0 else 0
                    altCount = rna_expression.loc[(rna_expression['contig'] == row['CHROM']) & (rna_expression['position'] == row['POS']) & (rna_expression['refAllele'] == row['REF']) & (rna_expression['altAllele'] == row['ALT']), 'altCount'].values[0] if len(rna_expression.loc[(rna_expression['contig'] == row['CHROM']) & (rna_expression['position'] == row['POS']) & (rna_expression['refAllele'] == row['REF']) & (rna_expression['altAllele'] == row['ALT'])]) > 0 else 0
                
                # Expression threshold: altCount > 5
                    expressed = True if (altCount) > expression_threshold else False       
                if sv_class:
                    somatic_mutation_annotations.append(self.create_somatic_mutation_annotation(row, pid, sample_id, gene, alteration, consequence, nMinor, nMajor, lohstatus, expHomAF, expHomCI_lo, expHomCI_hi, expHom_pbinom_lower, homogenous, ad0, ad1, ensemble_id=ensemble_id,depth=depth, AM_class=AM_class, AM_score=AM_score, polyphen_score=polyphen_score, polyphen_category=poylphen_category, sift_score=sift_score, sift_category=sift_category, classification=sv_class, consensus_pathogenecity=consensus_pathogenecity, consensus_pathogenecity_source=consensus_prediction_source, expressed=expressed, refCount=refCount, altCount=altCount))  

        return somatic_mutation_annotations



    def post_filter_and_classify_somatic_mutations(self, row):
        """
                Post-filter and classify somatic_mutations based on various criteria.

                Parameters:
                row (Series): A row from a DataFrame containing somatic_mutation data.

                Returns:
                Series: A Series containing the updated somatic_mutation data.
        """
        try:
            exonicFuncMane = handle_string_field(row["exonicFuncMane"])
            sample_id = row['sample_id']

            # Calculate homogeneity estimate
            tfs = self.ascats.loc[self.ascats['sample'] == sample_id]['purity']
            tf = tfs.iloc[0] if len(tfs) > 0 else 0.0  # loc[ascats['sample'] == sample_id]['purity'].values[0]
            ad0 = int(row['ad0'])
            ad1 = int(row['ad1'])
            depth = ad0 + ad1
            gene = handle_string_field(row["hugoSymbol"])

            nMajor = None
            nMinor = None
            lohstatus = None
            vcnas = self.get_variant_assoc_cnas(self.cnas, sample_id, gene)
            nMajor = handle_cn_field(vcnas['nMajor']) if len(vcnas) > 0 else None
            nMinor = handle_cn_field(vcnas['nMinor']) if len(vcnas) > 0 else None
            lohstatus = vcnas['LOHstatus'] if len(vcnas) > 0 else None
            expHomAF = 0.0
            expHomCI_lo = 0.0
            expHomCI_hi = 0.0
            expHom_pbinom_lower = 0.0
            homogenous = None

            if nMajor and nMinor:
                cn = int(nMinor) + int(nMajor)
                expHomAF = float(self.expectedAF(cn, cn, tf))
                expHomCI_lo = float(stats.binom.ppf(0.025, depth, expHomAF))
                expHomCI_hi = float(stats.binom.ppf(0.975, depth, expHomAF))
                expHomCI_cover = expHomCI_lo <= ad1
                expHom_pbinom_lower = float(stats.binom.cdf(ad1, depth, expHomAF))
                homogenous = expHom_pbinom_lower > self.homogeneity_threshold


            row['nMinor'] = nMinor
            row['nMajor'] = nMajor
            row['lohstatus'] = lohstatus
            row['hom_lo'] = "{:.9f}".format(expHomCI_lo),
            row['hom_hi'] = "{:.9f}".format(expHomCI_hi),
            row['hom_pbinom_lo'] = "{:.9f}".format(expHom_pbinom_lower),
            row['homogenous'] = bool(homogenous)
            row['af'] = expHomAF
            if homogenous and exonicFuncMane == "nonsynonymous_SNV":
                row['classification'] = "Missense"
            row['hom_lo'] = float(row['hom_lo'][0] if isinstance(row['hom_lo'], tuple) else row['hom_lo'])
            row['hom_hi'] = float(row['hom_hi'][0] if isinstance(row['hom_hi'], tuple) else row['hom_hi'])
            row['hom_pbinom_lo'] = float(row['hom_pbinom_lo'][0] if isinstance(row['hom_pbinom_lo'], tuple) else row[
                'hom_pbinom_lo'])

        except Exception as e:
            print(e)
            pass

        return row

