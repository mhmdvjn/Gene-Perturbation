"""Manual drug-target annotation data for the DruggableHypotheses project.

Format: 'GENE': (['drug1', ...], 'status')
Status options: 'approved_drug', 'clinical_or_indirect', 'no_coverage'
"""

MANUAL_DRUG_DATA = {
    # --- Core Kinases & Oncogenes (Approved Direct Drugs) ---
    'EGFR': (['erlotinib', 'gefitinib', 'osimertinib', 'afatinib', 'lapatinib', 'dacomitinib'], 'approved_drug'),
    'HER2': (['trastuzumab', 'pertuzumab', 'tucatinib', 'neratinib'], 'approved_drug'),
    'ERBB2': (['trastuzumab', 'pertuzumab', 'tucatinib', 'neratinib', 'lapatinib'], 'approved_drug'),
    'BRAF': (['vemurafenib', 'dabrafenib', 'encorafenib'], 'approved_drug'),
    'MEK1': (['trametinib', 'cobimetinib', 'binimetinib'], 'approved_drug'),
    'MAP2K1': (['trametinib', 'cobimetinib', 'binimetinib'], 'approved_drug'),
    'KRAS': (['sotorasib', 'adagrasib'], 'approved_drug'),
    'ALK': (['crizotinib', 'alectinib', 'ceritinib', 'brigatinib', 'lorlatinib'], 'approved_drug'),
    'ROS1': (['crizotinib', 'entrectinib'], 'approved_drug'),
    'RET': (['selpercatinib', 'pralsetinib'], 'approved_drug'),
    'MET': (['crizotinib', 'capmatinib', 'tepotinib', 'savolitinib'], 'approved_drug'),
    'KIT': (['imatinib', 'sunitinib', 'regorafenib', 'avapritinib', 'ripretinib'], 'approved_drug'),
    'PDGFRA': (['imatinib', 'avapritinib', 'ripretinib'], 'approved_drug'),
    'PDGFRB': (['imatinib', 'sunitinib'], 'approved_drug'),
    'FLT3': (['midostaurin', 'gilteritinib', 'quizartinib', 'crenolanib'], 'approved_drug'),
    'VEGFR2': (['sorafenib', 'sunitinib', 'pazopanib', 'lenvatinib', 'axitinib'], 'approved_drug'),
    'KDR': (['sorafenib', 'sunitinib', 'pazopanib', 'lenvatinib', 'axitinib'], 'approved_drug'),
    'VEGFR1': (['sorafenib', 'sunitinib'], 'approved_drug'),
    'FLT1': (['sorafenib', 'sunitinib'], 'approved_drug'),
    'VEGFR3': (['sorafenib', 'sunitinib', 'pazopanib'], 'approved_drug'),
    'FLT4': (['sorafenib', 'sunitinib', 'pazopanib'], 'approved_drug'),
    'CDK4': (['palbociclib', 'ribociclib', 'abemaciclib', 'trilaciclib'], 'approved_drug'),
    'CDK6': (['palbociclib', 'ribociclib', 'abemaciclib', 'trilaciclib'], 'approved_drug'),
    'PARP1': (['olaparib', 'rucaparib', 'niraparib', 'talazoparib', 'veliparib'], 'approved_drug'),
    'PARP2': (['olaparib', 'rucaparib'], 'approved_drug'),
    'ATR': (['ceralasertib', 'berzosertib', 'elimusertib'], 'approved_drug'),
    'CHEK1': (['prexasertib', 'rabusertib'], 'approved_drug'),
    'WEE1': (['adavosertib'], 'approved_drug'),
    'MTOR': (['everolimus', 'temsirolimus', 'sirolimus', 'ridaforolimus'], 'approved_drug'),
    'AKT1': (['capivasertib', 'ipatasertib'], 'approved_drug'),
    'AKT2': (['capivasertib'], 'approved_drug'),
    'PIK3CA': (['alpelisib', 'copanlisib', 'idelalisib', 'duvelisib'], 'approved_drug'),
    'BCR::ABL1': (['imatinib', 'dasatinib', 'nilotinib', 'bosutinib', 'ponatinib'], 'approved_drug'),
    'ABL1': (['imatinib', 'dasatinib', 'nilotinib', 'bosutinib', 'ponatinib'], 'approved_drug'),
    'SRC': (['dasatinib', 'saracatinib'], 'approved_drug'),
    'BTK': (['ibrutinib', 'acalabrutinib', 'zanubrutinib', 'tirabrutinib'], 'approved_drug'),
    'JAK1': (['ruxolitinib', 'fedratinib', 'pacritinib', 'momelotinib'], 'approved_drug'),
    'JAK2': (['ruxolitinib', 'fedratinib', 'pacritinib', 'momelotinib'], 'approved_drug'),
    'JAK3': (['tofacitinib', 'peficitinib'], 'approved_drug'),
    'TYK2': (['deucravacitinib'], 'approved_drug'),
    'SYK': (['fostamatinib'], 'approved_drug'),
    'IGF1R': (['linsitinib'], 'approved_drug'),
    'FGFR1': (['erdafitinib', 'pemigatinib', 'infigratinib', 'futibatinib'], 'approved_drug'),
    'FGFR2': (['erdafitinib', 'pemigatinib', 'infigratinib', 'futibatinib'], 'approved_drug'),
    'FGFR3': (['erdafitinib', 'pemigatinib', 'infigratinib', 'futibatinib'], 'approved_drug'),
    'FGFR4': (['fisogatinib'], 'approved_drug'),
    'NTRK1': (['larotrectinib', 'entrectinib'], 'approved_drug'),
    'NTRK2': (['larotrectinib', 'entrectinib'], 'approved_drug'),
    'NTRK3': (['larotrectinib', 'entrectinib'], 'approved_drug'),
    'MDM2': (['idasanutlin', 'milademetan'], 'approved_drug'),
    'IDH1': (['ivosidenib', 'olutasidenib'], 'approved_drug'),
    'IDH2': (['enasidenib'], 'approved_drug'),
    'EZH2': (['tazemetostat'], 'approved_drug'),
    'HDAC1': (['vorinostat', 'romidepsin', 'panobinostat', 'entinostat'], 'approved_drug'),
    'HDAC2': (['vorinostat', 'romidepsin', 'panobinostat'], 'approved_drug'),
    'HDAC6': (['ricolinostat'], 'approved_drug'),
    'DNMT1': (['azacitidine', 'decitabine'], 'approved_drug'),
    'DNMT3A': (['azacitidine', 'decitabine'], 'approved_drug'),
    'BCL2': (['venetoclax', 'navitoclax'], 'approved_drug'),
    'BCL2L1': (['venetoclax', 'navitoclax'], 'approved_drug'),
    'AR': (['enzalutamide', 'apalutamide', 'darolutamide', 'bicalutamide', 'abiraterone'], 'approved_drug'),
    'ESR1': (['tamoxifen', 'fulvestrant', 'anastrozole', 'letrozole', 'exemestane'], 'approved_drug'),
    'PGR': (['medroxyprogesterone', 'megestrol'], 'approved_drug'),
    'RARA': (['tretinoin', 'alitretinoin'], 'approved_drug'),
    'VEGFA': (['bevacizumab', 'aflibercept', 'ranibizumab'], 'approved_drug'),
    'ANGPT2': (['faricimab'], 'approved_drug'),
    'TGFBR1': (['galunisertib'], 'approved_drug'),
    'NOTCH1': (['nirogacestat'], 'approved_drug'),

    # --- Candidate Genes (Clinical / Synthetic Lethal / No Direct Coverage) ---
    'ARID1A': (['tazemetostat (synthetic lethal)', 'ceralasertib (synthetic lethal)'], 'clinical_or_indirect'),
    'BCL2L11': (['venetoclax (pathway mimetic)'], 'clinical_or_indirect'),
    'CEBPB': (['ST101 (Phase 2 antagonist)'], 'clinical_or_indirect'),
    'BCORL1': ([], 'no_coverage'),
    'C3ORF72': ([], 'no_coverage'),
    'CBL': ([], 'no_coverage'),
    'CDKN1A': ([], 'no_coverage'),
    'CDKN1B': ([], 'no_coverage'),
    'CDKN1C': ([], 'no_coverage'),
    'CELF2': ([], 'no_coverage'),
    'CKS1B': ([], 'no_coverage'),
    'COL2A1': ([], 'no_coverage'),
    'DUSP9': ([], 'no_coverage'),
    'EGR1': ([], 'no_coverage'),
    'FOXL2': ([], 'no_coverage'),

    # --- Other Untargeted / Tumor Suppressors ---
    'NRAS': ([], 'no_coverage'),
    'HRAS': ([], 'no_coverage'),
    'CDK2': ([], 'no_coverage'),
    'BRCA1': ([], 'no_coverage'),
    'BRCA2': ([], 'no_coverage'),
    'ATM': ([], 'no_coverage'),
    'CHEK2': ([], 'no_coverage'),
    'AKT3': ([], 'no_coverage'),
    'PIK3CB': ([], 'no_coverage'),
    'PTEN': ([], 'no_coverage'),
    'TP53': ([], 'no_coverage'),
    'CDKN2A': ([], 'no_coverage'),
    'CDKN2C': ([], 'no_coverage'),
    'RB1': ([], 'no_coverage'),
    'MYC': ([], 'no_coverage'),
}


def build_annotation_dicts():
    drug_targets = {}
    drug_status = {}
    for gene, (drugs, status) in MANUAL_DRUG_DATA.items():
        clean = gene.strip().upper()
        if clean not in drug_targets:
            drug_targets[clean] = list(drugs)
            drug_status[clean] = status
    return drug_targets, drug_status
