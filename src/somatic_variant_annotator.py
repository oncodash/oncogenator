from utils import *
from scipy import stats


class SomaticVariantAnnotator:
    def __init__(self, refgenome="GRCh38", tumortype="HGSOC", cnas=None, ascats=None, samples=None, homogeneity_threshold=0.05, ada_score_threshold=0.95, rf_score_threshold=0.95):
        """
                Initialize the SomaticVariantAnnotator class.

                Parameters:
                refgenome (str): Reference genome version, default is "GRCh38".
                tumortype (str): Type of tumor, default is "HGSOC".
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

    def create_somatic_mutation_annotation(self, row, pid, sample_id, gene, alteration, consequence, nMinor, nMajor, lohstatus, expHomAF, expHomCI_lo, expHomCI_hi, expHom_pbinom_lower, homogenous, ad0, ad1, depth, AM_category, amisscore, classification):
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
                AM_category (str): AMIS category.
                amisscore (float): AMIS score.

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
            'hugoSymbol': gene,
            'alteration': alteration,
            'tumorType': "HGSOC",
            'consequence': consequence,
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
            'ada_score': handle_decimal_field(row["dbscsomatic_mutation_ADA_SCORE"]),
            'rf_score': handle_decimal_field(row["dbscsomatic_mutation_RF_SCORE"]),
            'sift_category': handle_string_field(row["SIFTval"]),
            'sift_score': handle_string_field(row["SIFTcat"]),
            'polyphen_category': handle_string_field(row["PolyPhenCat"]),
            'polyphen_score': handle_string_field(row["PolyPhenVal"]),
            'AM_category': handle_string_field(AM_category),
            'AM_score': handle_decimal_field(amisscore),
            'cosmic_id': handle_string_field(row["COSMIC_ID"]),
            'clinvar_id': handle_string_field(row["CLNALLELEID"]),
            'classification': handle_string_field(classification)
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
        AM_variant = handle_string_field(row["AM_variant"])
        amisscore = row['AM_score']
        AM_category = row['AM_class']
        pathogenecity = AM_category

        exonicFuncMane = handle_string_field(row["ExonicFunc.MANE"])
        funcMane = handle_string_field(row["Func.MANE"])
        funcRefgene = handle_string_field(row["Func.refGene"])
        ada_score = handle_decimal_field(row["dbscsomatic_mutation_ADA_SCORE"])
        rf_score = handle_decimal_field(row["dbscsomatic_mutation_RF_SCORE"])

        for sample_id in self.samples:
            pid = sample_id.split("_")[0]

            tfs = self.ascats.loc[self.ascats['sample'] == sample_id]['purity']
            tf = tfs.iloc[0] if len(tfs) > 0 else 0.0
            depth = int(row[str(sample_id)+".DP"])
            ad0 = int(row[str(sample_id)+".AD"].split(',')[0])
            ad1 = int(row[str(sample_id)+".AD"].split(',')[1])

            geneMANE = handle_string_field(row["Gene.MANE"]).split(';')
            genes = set(geneMANE)
            geneRefGene = handle_string_field(row["Gene.refGene"]).split(';')
            for g in geneRefGene:
                genes.add(g)

            for gene in genes:
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
                    cn = float(nMinor) + float(nMajor)
                    expHomAF = float(self.expectedAF(cn, cn, tf))
                    expHomCI_lo = float(stats.binom.ppf(0.025, depth, expHomAF))
                    expHomCI_hi = float(stats.binom.ppf(0.975, depth, expHomAF))
                    expHomCI_cover = expHomCI_lo <= ad1
                    expHom_pbinom_lower = float(stats.binom.cdf(ad1, depth, expHomAF))
                    homogenous = expHom_pbinom_lower > self.homogeneity_threshold

                if homogenous and exonicFuncMane == "nonsynonymous_somatic_mutation":
                    sv_class = "Missense"
                if exonicFuncMane in ["frameshift_insertion", "frameshift_deletion", "stopgain"]:
                    sv_class = "Truncating"
                if exonicFuncMane == ["nonframeshift_deletion", "nonframeshift_substitution", "nonframeshift_insertion"]:
                    sv_class = "Other"
                consequence = exonicFuncMane
                if not sv_class:
                    if funcMane in ["splicing", "splicesite", "intron", "intronic"] or funcRefgene in ["splicing", "splicesite", "intron", "intronic"]:
                        consequence = funcMane
                        if (ada_score and float(ada_score) > self.ada_score_threshold) or (rf_score and float(rf_score) > self.rf_score_threshold):
                            sv_class = "Splicing"
                        else:
                            continue

                alteration = f"{gene}:{row['CHROM']}:{row['POS']}:{row['REF']}>{row['ALT']}"
                if sv_class:
                    somatic_mutation_annotations.append(self.create_somatic_mutation_annotation(row, pid, sample_id, gene, alteration, consequence, nMinor, nMajor, lohstatus, expHomAF, expHomCI_lo, expHomCI_hi, expHom_pbinom_lower, homogenous, ad0, ad1, depth, AM_category, amisscore, sv_class))

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
            if homogenous and exonicFuncMane == "nonsynonymous_somatic_mutation":
                row['classification'] = "Missense"
            row['hom_lo'] = float(row['hom_lo'][0] if isinstance(row['hom_lo'], tuple) else row['hom_lo'])
            row['hom_hi'] = float(row['hom_hi'][0] if isinstance(row['hom_hi'], tuple) else row['hom_hi'])
            row['hom_pbinom_lo'] = float(row['hom_pbinom_lo'][0] if isinstance(row['hom_pbinom_lo'], tuple) else row[
                'hom_pbinom_lo'])

        except Exception as e:
            print(e)
            pass

        return row


    def post_filter_and_classify_somatic_mutations_by_sample(self, row):
        """
                Post-filter and classify somatic_mutations by sample based on various criteria.

                Parameters:
                row (Series): A row from a DataFrame containing somatic_mutation data.

                Returns:
                Series: A Series containing the updated somatic_mutation data.
        """
        try:

            exonicFuncMane = handle_string_field(row["exonicFuncMane"])
            sample_id = row['sample_id']

            # Calculate homogeneity estimate
            tfs = self.ascats.loc[self.ascats['sample'] == sample_id]['purity'] # Tumor fraction = purity of the tumor
            tf = tfs.iloc[0] if len(tfs) > 0 else 0.0  # loc[ascats['sample'] == sample_id]['purity'].values[0]
            ad0 = int(row['ad0'])
            ad1 = int(row['ad1'])
            depth = ad0 + ad1
            gene = handle_string_field(row["hugoSymbol"])

            nMajor = None
            nMinor = None
            lohstatus = None
            pathogenecity = handle_string_field(row["CLINSIG"])
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

            if homogenous == True and row['classification'] == "Truncating":
                pathogenecity = "Pathogenic"

            row['nMinor'] = nMinor
            row['nMajor'] = nMajor
            row['lohstatus'] = lohstatus
            row['hom_lo'] = "{:.9f}".format(expHomCI_lo),
            row['hom_hi'] = "{:.9f}".format(expHomCI_hi),
            row['hom_pbinom_lo'] = "{:.9f}".format(expHom_pbinom_lower),
            row['homogenous'] = homogenous
            row['pathogenecity'] = pathogenecity
            row['af'] = expHomAF
            if homogenous and exonicFuncMane == "nonsynonymous_somatic_mutation":
                row['classification'] = "Missense"
            row['hom_lo'] = row['hom_lo'][0] if isinstance(row['hom_lo'], tuple) else row['hom_lo']
            row['hom_hi'] = row['hom_hi'][0] if isinstance(row['hom_hi'], tuple) else row['hom_hi']
            row['hom_pbinom_lo'] = row['hom_pbinom_lo'][0] if isinstance(row['hom_pbinom_lo'], tuple) else row[
                'hom_pbinom_lo']

        except Exception as e:
            print(e)
            pass

        return row
