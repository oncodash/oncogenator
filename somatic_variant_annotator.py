from utils import *
from scipy import stats

class somatic_variant_annotator:
    def __init__(self, models, files):
        self.models = models
        self.files = files


def get_variant_assoc_cnas(cnas, sid, gene):
    cnar = cnas.loc[(cnas['sample'] == sid) & (cnas['Gene'] == gene)]
    return cnar.iloc[0] if not cnar.empty else []


def expectedAF(N_t, CN_t, TF):
    return (N_t * TF) / (CN_t * TF + 2 * (1 - TF))


def parse_isoforms(aaChangeRefGene):
    records = aaChangeRefGene.split(",")
    isoforms = []
    for rec in records:
        fields = rec.split(":")
        if len(fields) > 1:
            gene = fields[0]
            protein = fields[len(fields) - 1].split(".")[1]
            isoform = gene + ":" + protein
            isoforms.append(isoform)

            # else aaChangeRefGene can have UNKNOWN status, what to do with it?
    return list(dict.fromkeys(isoforms))


def filter_and_classify_snvs(row, cnas=None, ascats=None, samples=None):
    sv_class = None

    snv_annotations = []

    sift_category = None
    sift_score = None  # Will be dded to WGS output in near future https://sift.bii.a-star.edu.sg/
    polyphen_category = None
    polyphen_score = None  # Will be dded to WGS output in near future http://genetics.bwh.harvard.edu/pph2/
    pathogenecity = handle_string_field(row["CLNSIG"])  # Currently Clinvar annotation is included in WGS pipeline output
    # TODO: Change pathogenecity annotation to voting after all scorings are included in output
    # if len(str(row['AAChange.refGene'])) > 2:
    #     isoforms = parse_isoforms(str(row['AAChange.refGene']))
    #     for isof in isoforms:
    #         ps = isof.split(':')
    #         protein = ps[1]
    #         amis = alphamissenses.loc[alphamissenses['protein_variant'] == protein]
    #         amisscore = amis.iloc[0]['am_pathogenicity'] if len(amis) > 0 else 0.0
    #         # TODO: vote pathogenecity by multiple estimates from sift,polyphen,oncokb,cgi
    #         amis_category = amis.iloc[0]['am_class'] if len(amis) > 0 else None
    #         pathogenecity = amis_category
    amisscore = row['AM_score']
    amis_category = row['AM_class']
    pathogenecity = amis_category


    #samples = handle_string_field(row["samples"]).split(';') if handle_string_field(row["samples"]) else []
    #readcounts = handle_string_field(row["readCounts"]).split(';') if handle_string_field(row["readCounts"]) else []
    i = 0

    exonicFuncMane = handle_string_field(row["ExonicFunc.MANE"])
    funcMane = handle_string_field(row["Func.MANE"])
    funcRefgene = handle_string_field(row["Func.refGene"])
    ada_score = handle_decimal_field(row["dbscSNV_ADA_SCORE"])
    rf_score = handle_decimal_field(row["dbscSNV_RF_SCORE"])
    for sample_id in samples:
        pid = sample_id.split("_")[0]

        # Splicing mutations
        #  funcMane | funRefgene = splicing / splicesite / intronic = > dbscSNV_ADA / RF_SCORE > 0.95
        # funcMane = handle_string_field(row["Func.MANE"])
        # funcRefgene = handle_string_field(row["Func.refGene"])

        # Truncating and Missenses => test homogeneity
        # truncating

        # if exonicFuncMane == ("frameshift_insertion" or "frameshift_deletion" or "stopgain" or "nonsynonymous_SNV"):
        # Calculate homogeneity estimate
        tfs = ascats.loc[ascats['sample'] == sample_id]['purity']
        tf = tfs.iloc[0] if len(tfs) > 0 else 0.0  # loc[ascats['sample'] == sample_id]['purity'].values[0]
        depth = int(row[str(sample_id)+".DP"])
        ad0 = int(row[str(sample_id)+".AD"].split(',')[0])
        ad1 = int(row[str(sample_id)+".AD"].split(',')[1])

        geneMANE = handle_string_field(row["Gene.MANE"]).split(';')
        genes = set(geneMANE)
        geneRefGene = handle_string_field(row["Gene.refGene"]).split(';')
        for g in geneRefGene:
            genes.add(g)
        for gene in genes:
            nMajor = None
            nMinor = None
            lohstatus = None
            vcnas = get_variant_assoc_cnas(cnas, sample_id, gene)
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
                expHomAF = float(expectedAF(cn, cn, tf))
                expHomCI_lo = float(stats.binom.ppf(0.025, depth, expHomAF))
                expHomCI_hi = float(stats.binom.ppf(0.975, depth, expHomAF))
                expHomCI_cover = expHomCI_lo <= ad1
                expHom_pbinom_lower = float(stats.binom.cdf(ad1, depth, expHomAF))
                homogenous = expHom_pbinom_lower > 0.05
            # TODO: Check FGFR2 in H266, why not estimated as homogenous
            if homogenous and exonicFuncMane == "nonsynonymous_SNV":
                sv_class = "Missense"
            if exonicFuncMane in ["frameshift_insertion", "frameshift_deletion", "stopgain"]:
                sv_class = "Truncating"
            if exonicFuncMane == ["nonframeshift_deletion", "nonframeshift_substitution", "nonframeshift_insertion"]:
                sv_class = "Other"
            consequence = exonicFuncMane
            if not sv_class:
                if funcMane in ["splicing", "splicesite", "intron", "intronic"] or funcRefgene in ["splicing",
                                                                                                   "splicesite",
                                                                                                   "intron",
                                                                                                   "intronic"]:
                    consequence = funcMane
                    if (ada_score and float(ada_score) > 0.95) or (rf_score and float(rf_score) > 0.95):
                        sv_class = "Splicing"
                    else:
                        continue

            alteration = gene + ":" + handle_string_field(row["CHROM"]) + ":" + str(
                handle_int_field(row["POS"])) + ":" + handle_string_field(row["REF"]) + ">" + handle_string_field(
                row["ALT"])
            if sv_class:
                snv_annotations.append(pd.Series({
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
                    # ensemblGeneId ': "",
                    # alteration: string;
                    'tumorType': "HGSOC",
                    'consequence': consequence,
                    'cytoBand': handle_string_field(row['cytoBand']),
                    'exonicFuncMane': exonicFuncMane,
                    'funcMane':handle_string_field(row['Func.MANE']),
                    'AAChangeMANE':handle_string_field(row['AAChange.MANE']),
                    'funcRefGene':handle_string_field(row['Func.refGene']),
                    'exonicFuncRefGene':handle_string_field(row['ExonicFunc.refGene']),
                    'AAChangerefGene':handle_string_field(row['AAChange.refGene']),
                    # proteinStart: string;
                    # proteinEnd: string;
                    # oncogenic: string;
                    # mutationEffectDescription: string;
                    # gene_role: string;
                    # citationPMids: string;
                    # geneSummary: string;
                    # variantSummary: string;
                    # tumorTypeSummary: string;
                    # diagnosticSummary: string;
                    # diagnosticImplications: string;
                    # prognosticImplications: string;
                    # treatments: string;
                    'nMinor': nMinor,
                    'nMajor': nMajor,
                    # oncokb_level: string;
                    # cgi_level: string;
                    # rank: number;
                    'ad0': ad0,
                    'ad1': ad1,
                    'af': expHomAF,
                    'depth': depth,
                    'lohstatus': lohstatus,
                    'hom_lo': "{:.9f}".format(expHomCI_lo),
                    'hom_hi': "{:.9f}".format(expHomCI_hi),
                    'hom_pbinom_lo': "{:.9f}".format(expHom_pbinom_lower),
                    'homogenous': homogenous,
                    # 'funcMane':handle_string_field(row["Func.MANE"]),
                    # 'funcRefgene':handle_string_field(row["Func.refGene"]),
                    # 'exonicFuncMane':handle_string_field(row["ExonicFunc.MANE"]),
                    'cadd_score': handle_decimal_field(row["CADD_phred"]),
                    'ada_score': ada_score,
                    'rf_score': rf_score,
                    'sift_category': sift_category,
                    'sift_score': sift_score,  # Will be dded to WGS output in near future https://sift.bii.a-star.edu.sg/
                    'polyphen_category': polyphen_category,
                    'polyphen_score': polyphen_score,
                    # Will be dded to WGS output in near future http://genetics.bwh.harvard.edu/pph2/
                    'amis_category': amis_category,
                    'amis_score': handle_decimal_field(amisscore),
                    'cosmic_id': handle_string_field(row["COSMIC_ID"]),
                    'clinvar_id': handle_string_field(row["CLNALLELEID"]),
                    'clinvar_sig': handle_string_field(row["CLNSIG"]),
                    'clinvar_status': handle_string_field(row["CLNREVSTAT"]),
                    'clinvar_assoc': handle_string_field(row["CLNDN"]),
                    'pathogenecity': handle_string_field(pathogenecity),
                    'classification': sv_class,
                }))
        i += 1
    if len(snv_annotations) > 0:
        return snv_annotations

def post_filter_and_classify_snvs(row, cnas=None, ascats=None):
    try:
        exonicFuncMane = handle_string_field(row["exonicFuncMane"])
        sample_id = row['sample_id']

        # Calculate homogeneity estimate
        tfs = ascats.loc[ascats['sample'] == sample_id]['purity']
        tf = tfs.iloc[0] if len(tfs) > 0 else 0.0  # loc[ascats['sample'] == sample_id]['purity'].values[0]
        ad0 = int(row['ad0'])
        ad1 = int(row['ad1'])
        depth = ad0 + ad1
        gene = handle_string_field(row["hugoSymbol"])

        nMajor = None
        nMinor = None
        lohstatus = None
        vcnas = get_variant_assoc_cnas(cnas, sample_id, gene)
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
            expHomAF = float(expectedAF(cn, cn, tf))
            expHomCI_lo = float(stats.binom.ppf(0.025, depth, expHomAF))
            expHomCI_hi = float(stats.binom.ppf(0.975, depth, expHomAF))
            expHomCI_cover = expHomCI_lo <= ad1
            expHom_pbinom_lower = float(stats.binom.cdf(ad1, depth, expHomAF))
            homogenous = expHom_pbinom_lower > 0.05

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


def post_filter_and_classify_snvs_by_sample(row, cnas=None, ascats=None):
    try:

        exonicFuncMane = handle_string_field(row["exonicFuncMane"])
        sample_id = row['sample_id']

        # Calculate homogeneity estimate
        tfs = ascats.loc[ascats['sample'] == sample_id]['purity']
        tf = tfs.iloc[0] if len(tfs) > 0 else 0.0  # loc[ascats['sample'] == sample_id]['purity'].values[0]
        ad0 = int(row['ad0'])
        ad1 = int(row['ad1'])
        depth = ad0 + ad1
        gene = handle_string_field(row["hugoSymbol"])

        nMajor = None
        nMinor = None
        lohstatus = None
        pathogenecity = handle_string_field(row["clinvar_sig"])
        vcnas = get_variant_assoc_cnas(cnas, sample_id, gene)
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
            expHomAF = float(expectedAF(cn, cn, tf))
            expHomCI_lo = float(stats.binom.ppf(0.025, depth, expHomAF))
            expHomCI_hi = float(stats.binom.ppf(0.975, depth, expHomAF))
            expHomCI_cover = expHomCI_lo <= ad1
            expHom_pbinom_lower = float(stats.binom.cdf(ad1, depth, expHomAF))
            homogenous = expHom_pbinom_lower > 0.05
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
