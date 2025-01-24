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

    def create_snv_annotation(self, row, pid, sample_id, gene, alteration, consequence, nMinor, nMajor, lohstatus, expHomAF, expHomCI_lo, expHomCI_hi, expHom_pbinom_lower, homogenous, ad0, ad1, depth, amis_category, amisscore):
        """
                Create an SNV annotation.

                Parameters:
                row (Series): A row from a DataFrame containing SNV data.
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
                amis_category (str): AMIS category.
                amisscore (float): AMIS score.

                Returns:
                Series: A Series containing the SNV annotation.
        """
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
            'ada_score': handle_decimal_field(row["dbscSNV_ADA_SCORE"]),
            'rf_score': handle_decimal_field(row["dbscSNV_RF_SCORE"]),
            'sift_category': None,
            'sift_score': None,
            'polyphen_category': None,
            'polyphen_score': None,
            'amis_category': amis_category,
            'amis_score': handle_decimal_field(amisscore),
            'cosmic_id': handle_string_field(row["COSMIC_ID"]),
            'clinvar_id': handle_string_field(row["CLNALLELEID"]),
        })

    def filter_and_classify_snvs(self, row):
        """
                Filter and classify SNVs based on various criteria.

                Parameters:
                row (Series): A row from a DataFrame containing SNV data.

                Returns:
                list: A list of Series containing SNV annotations.
        """
        sv_class = None
        snv_annotations = []

        amisscore = row['AM_score']
        amis_category = row['AM_class']
        pathogenecity = amis_category

        exonicFuncMane = handle_string_field(row["ExonicFunc.MANE"])
        funcMane = handle_string_field(row["Func.MANE"])
        funcRefgene = handle_string_field(row["Func.refGene"])
        ada_score = handle_decimal_field(row["dbscSNV_ADA_SCORE"])
        rf_score = handle_decimal_field(row["dbscSNV_RF_SCORE"])

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
                    cn = int(nMinor) + int(nMajor)
                    expHomAF = float(self.expectedAF(cn, cn, tf))
                    expHomCI_lo = float(stats.binom.ppf(0.025, depth, expHomAF))
                    expHomCI_hi = float(stats.binom.ppf(0.975, depth, expHomAF))
                    expHomCI_cover = expHomCI_lo <= ad1
                    expHom_pbinom_lower = float(stats.binom.cdf(ad1, depth, expHomAF))
                    homogenous = expHom_pbinom_lower > self.homogeneity_threshold

                if homogenous and exonicFuncMane == "nonsynonymous_SNV":
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
                    snv_annotations.append(self.create_snv_annotation(row, pid, sample_id, gene, alteration, consequence, nMinor, nMajor, lohstatus, expHomAF, expHomCI_lo, expHomCI_hi, expHom_pbinom_lower, homogenous, ad0, ad1, depth, amis_category, amisscore))

        return snv_annotations



    def post_filter_and_classify_snvs(self, row):
        """
                Post-filter and classify SNVs based on various criteria.

                Parameters:
                row (Series): A row from a DataFrame containing SNV data.

                Returns:
                Series: A Series containing the updated SNV data.
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
            row['homogenous'] = homogenous
            row['af'] = expHomAF
            if homogenous and exonicFuncMane == "nonsynonymous_SNV":
                row['classification'] = "Missense"
            row['hom_lo'] = row['hom_lo'][0] if isinstance(row['hom_lo'], tuple) else row['hom_lo']
            row['hom_hi'] = row['hom_hi'][0] if isinstance(row['hom_hi'], tuple) else row['hom_hi']
            row['hom_pbinom_lo'] = row['hom_pbinom_lo'][0] if isinstance(row['hom_pbinom_lo'], tuple) else row[
                'hom_pbinom_lo']

        except Exception as e:
            print(e)
            pass

        return row


    def post_filter_and_classify_snvs_by_sample(self, row):
        """
                Post-filter and classify SNVs by sample based on various criteria.

                Parameters:
                row (Series): A row from a DataFrame containing SNV data.

                Returns:
                Series: A Series containing the updated SNV data.
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
            pathogenecity = handle_string_field(row["clinvar_sig"])
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
            if homogenous and exonicFuncMane == "nonsynonymous_SNV":
                row['classification'] = "Missense"
            row['hom_lo'] = row['hom_lo'][0] if isinstance(row['hom_lo'], tuple) else row['hom_lo']
            row['hom_hi'] = row['hom_hi'][0] if isinstance(row['hom_hi'], tuple) else row['hom_hi']
            row['hom_pbinom_lo'] = row['hom_pbinom_lo'][0] if isinstance(row['hom_pbinom_lo'], tuple) else row[
                'hom_pbinom_lo']

        except Exception as e:
            print(e)
            pass

        return row
