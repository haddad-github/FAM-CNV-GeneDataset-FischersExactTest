import pandas as pd
import scipy.stats as stats

###Étape 1: Lecture des échantillons

def lireFAM(fichierFAM):
        #1. lire fichier FAM
        dataFAM = pd.read_csv(fichierFAM, header=None, sep=" ")

        #2. garder fondateurs en assigant
        dataFAM.columns = ["fid", "iid", "father", "mother", "sex", "status"]
        dataFAM = dataFAM[ (dataFAM["father"] == "0") & (dataFAM["mother"] == "0") ]
        
        #3. creer deux ensembles, controles et cas
        dataFAM = dataFAM.set_index("iid")
        controles = dataFAM.index[dataFAM["status"] == 1]
        cases = dataFAM.index[dataFAM["status"] == 2]

        return controles, cases

dataFAM = lireFAM("dataset.fam")

#controles et cas
con = dataFAM[0]
cas = dataFAM[1]

###Étape 2: Lecture des CNV

def lireCNV(fichierCNV):
    #1. lire fichier CNV
    dataCNV = pd.read_csv(fichierCNV, sep="\s+")

    #2. exclure CNV indice de qualitee inferieur ou egal a 30
    dataCNV.columns = ["iid", "chrom", "start", "end", "cn", "bf"]
    dataCNV = dataCNV[dataCNV.bf > 30]

    return dataCNV

dataCNV = lireCNV("dataset.cnv")


#Étape 3: Lecture des gènes

def lireGenes(fichierGenes):
    genes = pd.read_csv(fichierGenes, sep="\t")

    return genes

genes = lireGenes("ensembl_grch38_genes_chr18.cleaned.txt")


#Étape 4: Cas et contrôles

def get_gene_counts(cnvs, cases, controls, genes):
    """Counts the number of CNVs for cases and controls overlapping each gene.

    Args:
        cnvs (pandas.DataFrame): the CNV calls.
        cases (set): the set of cases.
        controls (set): the set of controls.
        genes (pandas.DataFrame): the genes.

    Returns:
        dict: a Counter object which counts, for each genes, the number of CNV
              calls that it overlaps. The keys are genes, and the values are
              tuples where the first element is the number of CASES, and the
              second element is the number of CONTROLS.

    Note
    ----
        You can hypothesize that, for a given sample, there are no overlapping
        CNVs.

    """
    # The final counter
    gene_counter = {}

    # The cases and controls
    is_case = cnvs.iid.isin(cases)
    is_controls = cnvs.iid.isin(controls)

    # Cycling through all the genes
    for _, gene in genes.iterrows():
        # Getting the gene position
        gene_start = gene["Gene start (bp)"]
        gene_end = gene["Gene end (bp)"]

        # Overlapping with the CNVs
        overlapping_cnvs = find_overlapping_cnvs(
            gene_start=gene_start,
            gene_end=gene_end,
            cnvs_start=cnvs.start,
            cnvs_end=cnvs.end,
        )

        # Counting the number of cases and controls
        nb_cases = (overlapping_cnvs & is_case).sum()
        nb_controls = (overlapping_cnvs & is_controls).sum()

        # Saving the data
        gene_counter[gene["Gene name"]] = (nb_cases, nb_controls)

    return gene_counter


def find_overlapping_cnvs(gene_start, gene_end, cnvs_start, cnvs_end):
    """Finds overlapping CNV calls.

    Args:
        gene_start (int): the starting position of the gene.
        gene_end (int): the ending position of the gene.
        cnvs_start (pandas.Series): the starting positions of the CNV calls.
        cnvs_end (pandas.Series): the ending positions of the CNV calls.

    Returns:
        pandas.Series: a Series (pandas) containing True/False values, where
                       True if a CNV call overlaps the gene.

    There are different ways a CNV call can overlap a gene:

        1. The starting position of the CNV call (C) is inside the gene (G):
            GGGGGGGGGGGGGG
                    CCCCCCCCC
                CCCCC

        2. The ending position of the CNV call (C) is inside the gene (G):
                GGGGGGGGGGGGGG
            CCCCCCCCC
                        CCCC

        3. The gene is located inside the CNV call.
                GGGGGGGGGGGGGG
            CCCCCCCCCCCCCCCCCCCCCCCC


    Note
    ----
        We assume the gene and the CNV are located on the same chromosome.

    """
    # Start inside the gene?
    good_start = (cnvs_start >= gene_start) & (cnvs_start <= gene_end)

    # End inside the gene?
    good_end = (cnvs_end >= gene_start) & (cnvs_end <= gene_end)

    # The gene is located inside
    inside = (cnvs_start <= gene_start) & (cnvs_end >= gene_end)

    return good_start | good_end | inside

dataGenes = get_gene_counts(dataCNV, cas, con, genes)


###Étape 5

#liste a stoquer
geneList = []


#calculer de Fisher pour chaque gene
for gene in dataGenes:
    #chercher les cas
    casCNV = dataGenes[gene][0]
    conCNV = dataGenes[gene][1]
    
    #nombre de cas/con - CNV cas/con pour total de cas/non sans CNV 
    casNoCNV = len(cas) - casCNV
    conNoCNV = len(con) - conCNV
   
    oddsratio, pvalue = stats.fisher_exact([[casCNV, conCNV], [casNoCNV, conNoCNV]])
    
    geneList += [(gene, oddsratio, pvalue)]
    
df = pd.DataFrame(geneList, columns=["symbol", "odds", "p"])

###Étape 6

#isoler les dataframe ayant un pvalue plus petit que 1e-4
for gene in df:
    genesInt = df[df["p"] < 1e-4]
    

###Étape 7

#ecire les fichiers
with open("burden_results.txt", "w") as fichier:
    fichier.write(df.to_string())
    
with open("burden_results.interesting.txt", "w") as fichier2:
    fichier2.write(genesInt.to_string())
   
