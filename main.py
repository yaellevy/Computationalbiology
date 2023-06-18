import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from Bio.Seq import Seq
from Bio import SeqIO, pairwise2


class TransmitterStatistics_statistics:
    pass


def geneLength(flist, genet):
    gensLength = {}
    for i in range(len(flist)):
        f = flist[i]
        try:
            if f.type:
                if 'gene' in flist[i].qualifiers.keys() and f.type == genet:
                    start_s = flist[i].location.start.position
                    end_s = flist[i].location.end.position
                    genName = flist[i].qualifiers['gene'][0]
                    gensLength[genName] = end_s - start_s
        except TransmitterStatistics_statistics as e:
            d = f[0]
            if 'gene' in d.qualifiers.keys() and d.type == genet:
                start_s = d.location.start.position
                end_s = d.location.end.position
                genName = d.qualifiers['gene'][0]
                gensLength[genName] = end_s - start_s
    return gensLength


def transmemranStatistics(atDict):
    TransmemranStatistics_statistics = {}
    total = 0
    countlen = 0
    TransmemranStatistics_statistics['max_value'] = max(atDict.values())
    TransmemranStatistics_statistics['min_value'] = min(atDict.values())
    for v in atDict.values():
        countlen += len(v)
        total += sum(v)
    TransmemranStatistics_statistics['mean_value'] = total / countlen
    return TransmemranStatistics_statistics


def At_statistics(atDict):
    sensationalistic = {'max_value': max(atDict.values()), 'min_value': min(atDict.values()),
                        'mean_value': sum(atDict.values()) / len(atDict)}
    listr = []
    std = np.std
    # appending all the values in the list
    for value in atDict.values():
        listr.append(value)
        # # # calculating standard deviation using np.std
        std = np.std(listr)
    sensationalistic['deviation'] = std
    return sensationalistic


def gene_cds_statistics(listGenes, genet):
    genecdsstatistics = {}
    gennocdsdictionary = geneLength(listGenes, genet)
    genecdsstatistics['max_value'] = max(gennocdsdictionary.values())
    genecdsstatistics['min_value'] = min(gennocdsdictionary.values())
    genecdsstatistics['mean_value'] = sum(gennocdsdictionary.values()) / len(gennocdsdictionary)
    listr = []
    std = np.std
    # appending all the values in the list
    for value in gennocdsdictionary.values():
        listr.append(value)
        # # # calculating standard deviation using np.std
        std = np.std(listr)
    genecdsstatistics['deviation'] = std
    return genecdsstatistics


def histogram(genesDict, lbx, lby, bins):
    #  genesDict    = {}
    # genesDict    = self.geneLength(list,'gene')
    plt.hist(genesDict.values(), bins=bins, edgecolor='black')
    plt.xlabel(lbx)
    plt.ylabel(lby)
    plt.show()


def alignment(seq1, seq2):
    aligmentdict = {}
    aligmentdictscore = {}
    matchesdict = {}
    match = []
    matchesdict = {'match': 0, 'mismatch': 0, 'gaps_1': 0, 'gaps_2': 0, 'transitions': 0, 'transversions': 0}
    for a, b in zip(seq1, seq2):
        if a == b:
            match.append('|')
            matchesdict['match'] += 1
        else:
            match.append(' ')
    print(seq1)
    print("".join(match))
    print(seq2)


def countAT(seq):
    seq.upper()
    countat = seq.count('A')
    countat += seq.count('T')
    if countat == 0:
        return 0
    else:
        return (countat / len(seq)) * 100


def getMeanValue(geneDict):
    return sum(geneDict.values()) / len(geneDict)


def getDictOfGensCDS(lst):
    from collections import defaultdict
    dictvalue = defaultdict(int)

    for i in range(len(lst)):
        try:
            dictvalue[lst[i].qualifiers['gene'][0]] += 1
        except Exception as e:
            continue
    return dictvalue


def ATPercenrHistogram(dica, dicb, dicc):
    fig, axs = plt.subplots(nrows=2, ncols=2)
    axs[0, 0].hist(dica.values(), bins=10, edgecolor='black')
    axs[0, 0].set_title('Group A')
    axs[0, 0].set_xlabel('%AT')
    axs[0, 0].set_ylabel('frequency')

    axs[0, 1].hist(dicb.values(), bins=10, edgecolor='black')
    axs[0, 1].set_title('Group B')
    axs[0, 1].set_xlabel('%AT')
    axs[0, 1].set_ylabel('frequency')
    # plt.subplot(1, 2, 2)
    axs[1, 0].hist(dicc.values(), bins=10, edgecolor='black')
    axs[1, 0].set_title('Group A-B')
    axs[1, 0].set_xlabel('%AT')
    axs[1, 0].set_ylabel('frequency')
    axs[1, 1].hist(dicc.values(), bins=10, edgecolor='red')
    axs[1, 1].hist(dicb.values(), bins=10, edgecolor='blue')
    axs[1, 1].set_xlabel('%AT')
    axs[1, 1].set_ylabel('frequency')
    axs[1, 1].set_title('Group B and Group A-B')
    fig.tight_layout()
    plt.show()


class Genes:
    def __init__(self, table, file):
        self.hidrofobicacid = ['A', 'V', 'L', 'I', 'P', 'F', 'C']
        self.table = table
        self.overAllDict = {'gene': [], 'start': [], 'end': [], 'geneSeq': [], 'cdsSeq': [], 'strand': [], '%AT': []}
        self.list_of_gene_cds = []
        self.list_of_gene_no_cds = []
        self.list_of_gene_regulator = []
        self.list_of_gene_rna = []
        self.cdsGenes = []
        self.featuresCount = {}
        self.dfUniprot = pd
        geneList = []
        startList = []
        endList = []
        seqList = []
        cdsSeq = []
        strand = []
        ATPercent = []

        assert (os.path.exists(file))  # making sure that the path is valid
        with open(file, "rU") as input_handle:
            for i, record_gb in enumerate(SeqIO.parse(input_handle, "genbank")): ''
        # 👇️ this attribute hides the method
        self.sequence = record_gb.seq
        features_list = record_gb.features
        self.features_list = features_list
        for i in range(len(self.features_list)):
            j = i + 1
            if j < len(self.features_list):
                f = self.features_list[i]
                r = self.features_list[j]
                if f.type == 'gene' and r.type == 'CDS' and 'gene' in f.qualifiers.keys():
                    geneList.append(f.qualifiers['gene'][0])
                    startList.append(f.location.start.position)
                    endList.append(f.location.end.position)
                    strand.append(f.location.strand)
                    ATPercent.append(countAT(self.sequence[f.location.start.position:f.location.end.position]))
                    seqList.append(self.sequence[f.location.start.position:f.location.end.position])
                    cdsSeq.append(r.qualifiers['translation'][0])
                    self.list_of_gene_cds.append(f)
                    self.cdsGenes.append(r)
                elif f.type == 'gene':
                    if r.type == 'rRNA':
                        self.list_of_gene_rna.append(f)
                    else:
                        self.list_of_gene_regulator.append(f)
        self.list_of_gene_no_cds.extend(self.list_of_gene_rna)
        self.list_of_gene_no_cds.extend(self.list_of_gene_regulator)
        self.overAllDict['gene'] = geneList
        self.overAllDict['start'] = startList
        self.overAllDict['end'] = endList
        self.overAllDict['geneSeq'] = seqList
        self.overAllDict['cdsSeq'] = cdsSeq
        self.overAllDict['%AT'] = ATPercent
        self.overAllDict['strand'] = strand
        self.df = pd.DataFrame.from_dict(self.overAllDict)
        self.df.set_index('gene', drop=True, inplace=True)
        self.dfUniprot = pd.read_excel('Uniprot.xlsx')
        if not os.path.exists('part_a.csv'):
            self.df.to_csv('part_a.csv', header=True)


    def countFeatures(self):
        # featuresCount={}
        from collections import defaultdict
        featuresCount = defaultdict(int)
        for i in range(len(self.features_list)):
            if featuresCount[self.features_list[i].type] in featuresCount.keys():
                featuresCount[self.features_list[i].type] += 1
            else:
                featuresCount[self.features_list[i].type] = 1
        print(featuresCount)

    # 👇️ same name as class variable
    def get_features(self):
        return self.features_list

    def sortGene(self):
        from collections import defaultdict
        genbank = defaultdict(int)
        for i in range(len(self.features_list)):
            f = self.features_list[i]
            if f.type in genbank.keys():
                genbank[f.type] += 1
            else:
                genbank[f.type] = 1
        return genbank

    def gene_cds(self, cds):
        if cds:
            return self.list_of_gene_cds
        else:
            return self.list_of_gene_no_cds

    def compareGenesCds(self):
        gene = []
        compare = []
        coding_dna = ''
        compareDict = {'gene': [], 'error': []}

        for index in self.df.itertuples():
            gene_name = index.Index
            if index.strand == -1:
                try:
                    coding_dna = Seq(self.sequence[index.start:index.end]).reverse_complement().translate(
                        table=self.table, cds=True)
                except Exception as e:
                    compare.append(e)
                    gene.append(gene_name)

            else:
                try:
                    coding_dna = Seq(self.sequence[index.start:index.end]).translate(
                        table=self.table, cds=True)
                except Exception as e:
                    compare.append(e)
                    gene.append(gene_name)
            alignments1 = pairwise2.align.globalxx(coding_dna, index.cdsSeq)
        compareDict['gene'] = gene
        compareDict['error'] = compare
        df = pd.DataFrame.from_dict(compareDict)
        df.set_index('gene', drop=True, inplace=True)
        if not os.path.exists('gene_exceptions.csv'):
          df.to_csv('gene_exceptions.csv', header=True)

    def getSeq(self):
        return self.sequence

    def cellWall(self):  # ex4
        counting = 0
        listGenes = []
        genesDict = {}
        for i in range(len(self.cdsGenes)):
            for j in range(len(self.cdsGenes[i].qualifiers.values())):
                # if list(self.cdsGenes[i].qualifiers.values())[j][0] == 'yabE':
                result = list(self.cdsGenes[i].qualifiers.values())[j][0].find('cell wall')
                if result > 0:
                    counting += 1
                    listGenes.append(self.cdsGenes[i])

        return counting, listGenes

    def countGenesAT(self, listGenes, type):
        genName = ''
        seq = ''
        from collections import defaultdict
        countATGen = defaultdict(float)
        for i in range(len(listGenes)):
            f = listGenes[i]
            try:
                if f.type:
                    if 'gene' in f.qualifiers.keys() and f.type == type:
                        start_s = f.location.start.position
                        end_s = f.location.end.position
                        genName = f.qualifiers['gene'][0]
                        seq = self.sequence[start_s:end_s]

            except Exception as e:
                d = f[0]
                if 'gene' in d.qualifiers.keys() and d.type == type:
                    start_s = d.location.start.position
                    end_s = d.location.end.position
                    genName = d.qualifiers['gene'][0]
                    seq = self.sequence[start_s:end_s]
            countATGen[genName] = countAT(seq)
        return countATGen

    def getDictFiveTopFiveLess(self):
        listdict = []
        listSortdict = []
        for i in range(len(self.features_list)):
            gensDict = {}
            f = self.features_list[i]
            if 'gene' in self.features_list[i].qualifiers.keys() and f.type == 'gene':
                start_s = self.features_list[i].location.start.position
                end_s = self.features_list[i].location.end.position
                genName = self.features_list[i].qualifiers['gene'][0]
                seq = self.sequence[start_s:end_s]
                gensDict['name'] = genName
                gensDict['start'] = start_s
                gensDict['end'] = end_s
                gensDict['%AT'] = countAT(seq)
                gensDict['strand'] = self.features_list[i].location.strand
                listdict.append(gensDict)
        listSortdict.append((sorted(listdict, key=lambda k: k['%AT'], reverse=True)[:5]))
        listSortdict.append((sorted(listdict, key=lambda k: k['%AT'])[:5]))
        return listSortdict

    def uniProtCompare(self):
        listinuniprotandgb = []
        listinuniprot = []
        listingb = []
        uniprotGenes = {'geneBoth': 0, 'inUniport': 0, 'inGenbank': 0}
        for index in self.dfUniprot.itertuples():
            genesName = index._5.split()
            flagExist = ''
            for j in range(len(genesName)):
                if genesName[j] in self.df.index:
                    flagExist = 'X'
                    listinuniprotandgb.append(genesName[j])
                    break
            if flagExist != 'X':
                listinuniprot.append(index._9)
        for index in self.df.itertuples():
            Exist = self.dfUniprot.loc[self.dfUniprot['Gene Names'].str.contains(index.Index, case=False)]
            if Exist.empty == True:
                listingb.append(index.Index)
        uniprotGenes['geneBoth'] = len(listinuniprotandgb)
        uniprotGenes['inUniport'] = len(listinuniprot)
        uniprotGenes['inGenbank'] = len(listingb)
        print(uniprotGenes)
        return listinuniprotandgb

    def aminoAcidPercent(self, seq):
        sumAcid = 0
        for i in range(len(self.hidrofobicacid)):
            sumAcid += seq.count(self.hidrofobicacid[i])
        return sumAcid / len(seq)

    def getTransmambernal(self):
        from collections import defaultdict
        dicB = defaultdict(float)
        from collections import defaultdict
        dicA = defaultdict(float)
        from collections import defaultdict
        dicC = defaultdict(float)
        dictMamber = {}
        aminoacid = {}
        dictSeq = {}
        length = []
        seq = []
        aminoAcidList = []

        DicAll = getDictOfGensCDS(self.list_of_gene_no_cds)
        DictCds = getDictOfGensCDS(self.list_of_gene_cds)
        DicAll.update(DictCds)

        for index in self.dfUniprot.itertuples():
            if index._9 in DictCds.keys():
                dicA[index._9] = countAT(index.Sequence)
            if not pd.isna(index.Transmembrane):

                doOnce = 'X'
                mamberList = index.Transmembrane.split('TRANSMEM')
                positions = []
                for j in range(len(mamberList)):
                    positions.append(mamberList[j].split(';')[0])
                    length = []
                    seq = []
                    aminoAcidList = []
                for k in range(len(positions)):
                    start = positions[k].split('..')[0]

                    if start != '':
                        length.append(int(positions[k].split('..')[1]) - int(positions[k].split('..')[0]))
                        sequence = index.Sequence[int(positions[k].split('..')[0]):int(positions[k].split('..')[1])]
                        seq.append(sequence)
                        aminoAcidList.append(self.aminoAcidPercent(sequence))
                        if index._9 in DicAll.keys():
                            if doOnce == 'X':
                                doOnce = ''
                                dicB[index._9] = countAT(index.Sequence)
                dictMamber[index._9] = length
                dictSeq[index._9] = seq
                aminoacid[index._9] = aminoAcidList
                dicC = dict(set(dicA.items()) - set(dicB.items()))
        return dictMamber, dictSeq, aminoacid, dicA, dicB, dicC


Gene = Genes(11, 'BS168.gb')
# Part A
# 1.print(Gene.sortGene())
# 2.a print(geneLength(Gene.features_list,'gene'))
# 2.b print(geneLength(Gene.list_of_gene_cds,'gene'))
#     print(geneLength(Gene.list_of_gene_rna,'gene'))
#     print(geneLength(Gene.list_of_gene_no_cds,'gene'))
# 2.c
# print(gene_cds_statistics(Gene.list_of_gene_no_cds,'gene'))
# print(gene_cds_statistics(Gene.list_of_gene_cds,'gene'))
# 2.d
# histogram(geneLength(Gene.features_list, 'gene'), 'Length', 'Frequency', 200)
# histogram(geneLength(Gene.list_of_gene_cds, 'gene'), 'Length', 'Frequency', 200)
# histogram(geneLength(Gene.list_of_gene_no_cds, 'gene'), 'Length', 'Frequency', 200)
# 3.a
# print('mean of %AT in all Genom:',countAT(Gene.sequence))#mean %AT in all Genom
# 3.b
# print('mean of %AT in Genes that coding to CDS:',
#      getMeanValue(Gene.countGenesAT(Gene.list_of_gene_cds, 'gene')))  # mean of %AT in Genes that coding to CDS
# print('mean of %AT in Genes that Are Not Coding to CDS:', getMeanValue(
#     Gene.countGenesAT(Gene.list_of_gene_no_cds, 'gene')))  # mean of %AT in Genes that Not coding to CDS
# 3.d
# histogram(Gene.countGenesAT(Gene.list_of_gene_cds,'gene'),'Length','Frequency',200)
# 3.e
#print(Gene.getDictFiveTopFiveLess())
# 4.a
count, llistGen = Gene.cellWall()
cellWallDict = geneLength(llistGen, 'CDS')
print('num of cellWall =', count)
# 4.b
# count,llistGen = Gene.cellWall()
# cellWallDict   = geneLength(llistGen,'CDS')
# print(cellWallDict)
# print(gene_cds_statistics(llistGen,'CDS'))
# histogram(cellWallDict,'Length Cell Wall Genes','Frequency',200)
# 4.c
# count, llistGen = Gene.cellWall()
# print(Gene.countGenesAT(llistGen, 'CDS'))
# histogram(Gene.countGenesAT(llistGen, 'CDS'), '%AT For Cell Wall Genes', 'Frequency', 200)
# print(At_statistics(Gene.countGenesAT(llistGen, 'CDS')))
# 5
#Gene.compareGenesCds()

# ---------------------------------------------------------------------------
# Part B
# B.a
# Gene.uniProtCompare()
# B.b.1
# dictMamber,dictSeq,aminoacid,dicA,dicB,dicC = Gene.getTransmambernal()
# print(transmemranStatistics(dictMamber))
# histogram(dictMamber,'Length','Frequency',10)
# B.b.2
# dictMamber,dictSeq,aminoacid,dicA,dicB,dicC = Gene.getTransmambernal()
# print(transmemranStatistics(aminoacid))
# histogram(aminoacid,'Length','Frequency',10)
# B.c.1
# dictMamber,dictSeq,aminoacid,dicA,dicB,dicC = Gene.getTransmambernal()
# histogram(dicB,'%AT','Frequency',10)
# B.c.2
dictMamber,dictSeq,aminoacid,dicA,dicB,dicC = Gene.getTransmambernal()
ATPercenrHistogram(dicA,dicB,dicC)
print(At_statistics(dicB))
print(At_statistics(dicC))
