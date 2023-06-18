import os
from typing import List, Any
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from Bio.Seq import Seq

from Bio import Entrez, SeqIO, pairwise2
class Genes():
    def __init__(self):

        self.overAllDict = {'gene':[],'start':[],'end':[],'geneSeq':[],'cdsSeq':[],'geneCds':[]}
        self.list_of_gene_cds=[]
        self.list_of_gene_no_cds=[]
        self.cdsGenes = []
        geneList =[]
        startList=[]
        endList  =[]
        seqList  =[]
        cdsSeq  = []
        geneCds = []
        genecdsdict = {}
        from Bio import SeqIO
        gene_bank_file = 'BS168.gb'
        assert (os.path.exists(gene_bank_file))  # making sure that the path is valid
        with open(gene_bank_file, "rU") as input_handle:
            for i, record_gb in enumerate(SeqIO.parse(input_handle, "genbank")):''
        # 👇️ this attribute hides the method
        self.sequence = record_gb.seq
        features_list = record_gb.features

        self.features_list = features_list

        for i in range(len(self.features_list)):
            j = i + 1
            if (j < len(self.features_list)):
                f = self.features_list[i]
                r = self.features_list[j]
                if (f.type == 'gene' and r.type == 'CDS' and 'gene' in f.qualifiers.keys()):
                    geneList.append(f.qualifiers['gene'][0])
                    startList.append(f.location.start.position)
                    endList.append(f.location.end.position)
                    seqList.append(self.sequence[f.location.start.position:f.location.end.position])
                    cdsSeq.append(r.qualifiers['translation'][0])
                    self.list_of_gene_cds.append(f)
                    self.cdsGenes.append(r)
                    if f.location.strand == -1:
                        coding_dna = Seq(self.sequence[f.location.start.position:f.location.end.position]).reverse_complement().translate(table=11)
                    else:
                        coding_dna = Seq(self.sequence[f.location.start.position:f.location.end.position]).translate(table=11)
                    geneCds.append(str(coding_dna))
                elif (f.type == 'gene'):
                    self.list_of_gene_no_cds.append(f)
        self.overAllDict['gene']  =geneList
        self.overAllDict['start'] =startList
        self.overAllDict['end']   =endList
        self.overAllDict['geneSeq']=seqList
        self.overAllDict['cdsSeq'] =cdsSeq
        self.overAllDict['geneCds']=geneCds
        self.df  = pd.DataFrame.from_dict(self.overAllDict)

    # 👇️ same name as class variable
    def get_features(self):
        return self.features_list

    def geneLength(self,flist,type):
        count = 0
        gensLength={}
        for i in range(len(flist)):
            f = flist[i]
            if 'gene' in flist[i].qualifiers.keys() and f.type == type:
                start_s = flist[i].location.start.position
                end_s = flist[i].location.end.position
                genName = flist[i].qualifiers['gene'][0]
                gensLength[genName] = end_s - start_s
        return gensLength

    def sortGene(self):
        from collections import defaultdict
        genbank = defaultdict(int)
        for i in range(len(self.features_list)):
            f = self.features_list[i]
            if f.type in genbank.keys():
                genbank[f.type] += 1
        return genbank

    def gene_cds(self,cds):
        if cds == True:
            return self.list_of_gene_cds
        else:
            return self.list_of_gene_no_cds

    def gene_cds_statistics(self,cds):
        gene_cds_statistics = {}
        gennocdsdictionary  = {}
        listGenes=[]
        if cds == True:
            listGenes = self.list_of_gene_cds
        else:
            listGenes = self.list_of_gene_no_cds
        gennocdsdictionary = self.geneLength(listGenes,'gene')
        gene_cds_statistics['max_value'] = max(gennocdsdictionary.values())
        gene_cds_statistics['min_value'] = min(gennocdsdictionary.values())
        gene_cds_statistics['mean_value'] = sum(gennocdsdictionary.values()) / len(gennocdsdictionary)
        listr = []
        # appending all the values in the list
        for value in gennocdsdictionary.values():
           listr.append(value)
        # # # calculating standard deviation using np.std
           std = np.std(listr)
        gene_cds_statistics['deviation'] =std
        return gene_cds_statistics
    def histogram(self,genesDict):
      #  genesDict    = {}
       # genesDict    = self.geneLength(list,'gene')
        plt.hist(genesDict.values(), bins=100, edgecolor='black')
        plt.xlabel('length')
        plt.ylabel('Genes Count')
        plt.show()
    def compareGenesCds(self,table):
        for i in range (len(self.df)):
             alignments1 = pairwise2.align.globalxx(self.df.loc[i, "geneCds"],self.df.loc[i, "cdsSeq"])

    def alignment(self,seq1,seq2,gene):
        for a, b in zip(seq1, seq2):
            if a == b:
                match.append('|')
                matchesdict['match'] += 1
            else:
                match.append(' ')
        print(alignments1[0][0])
        print("".join(match))
        print(alignments1[0][1])
    def getSeq(self):
        return  self.sequence
    def countAT(self,seq):
        seq.upper()
        countat= seq.count('A')
        countat +=seq.count('T')
        if countat == 0:
          return -1
        else:
          return (countat / len(seq)) * 100

    def cellWall(self):  #ex4
        count = 0
        listGenes=[]
        genesDict={}
        for i in range(len(self.cdsGenes)):
            for j in range(len(self.cdsGenes[i].qualifiers.values())):
                #if list(self.cdsGenes[i].qualifiers.values())[j][0] == 'yabE':
                     result = list(self.cdsGenes[i].qualifiers.values())[j][0].find('cell wall')
                     if result > 0:
                        count += 1
                        listGenes.append(self.cdsGenes[i]) #.qualifiers['gene'][0]
        print(self.countGenesAT(listGenes))
        return self.geneLength(listGenes,'CDS'),len(listGenes)
        #fasta_path = 'sequences_dna.fasta'
        #for record in SeqIO.parse(fasta_path, "fasta"):
    def getMeanValue(self,geneDict):
       return sum(geneDict.values()) / len(geneDict)
    def countGenesAT(self,listGenes):
      from collections import defaultdict
      countATGen = defaultdict(float)
      for i in range(len(listGenes)):
          start=listGenes[i].location.start.position
          end=listGenes[i].location.end.position
          seq=self.sequence[start:end]
          if 'gene' in self.list_of_gene_cds[i].qualifiers.keys() and self.list_of_gene_cds[i].type == 'gene':
              gene = self.list_of_gene_cds[i].qualifiers['gene'][0]
              countATGen[gene]=self.countAT(seq)
          else:
            countATGen[list(self.list_of_gene_cds[i].qualifiers.keys())[0]] = self.countAT(seq)

      return countATGen
    def getDictFiveTopFiveLess(self):
        listdict=[]
        listSortdict= []
        for i in range(len(self.features_list)):
               gensDict = {}
               f = self.features_list[i]
               if 'gene' in self.features_list[i].qualifiers.keys() and f.type == 'gene':
                    start_s = self.features_list[i].location.start.position
                    end_s = self.features_list[i].location.end.position
                    genName = self.features_list[i].qualifiers['gene'][0]
                    seq     = self.sequence[start_s:end_s]
                    gensDict['name']  = genName
                    gensDict['start'] = start_s
                    gensDict['end']   = end_s
                    gensDict['%AT']   = self.countAT(seq)
                    gensDict['strand']= self.features_list[i].location.strand
                    listdict.append(gensDict)

        listSortdict.append((sorted(listdict, key=lambda k: k['%AT'],reverse=True)[:5]))
        listSortdict.append((sorted(listdict, key=lambda k: k['%AT'])[:5]))
        return listSortdict





a=Genes()
#print(a.compareGenesCds(11))
#print('mean of %AT in Genes that coding to CDS:',a.getMeanValue(a.countGenesAT(a.list_of_gene_cds)))#mean of %AT in Genes that coding to CDS
#print('mean of %AT in all Genom:',a.countAT(a.sequence))#mean %AT in all Genom
#a.histogram(a.countGenesAT(a.list_of_gene_cds))
#ex4
#c=a.cellWall()
#print(c)
#print('num of cellWall=',len(c))

#print(a.gene_cds_statistics(a.cellWall()))

"""
dict={}
print((a.gene_cds_statistics(False)))
print((a.gene_cds_statistics(True)))
print(a.countAT(a.getSeq()))
a.histogram(a.features_list)
a.histogram(a.list_of_gene_no_cds)
a.histogram(a.list_of_gene_cds)
"""
"""
dict = a.geneLength(a.list_of_gene_cds)

for name,i in dict.items():
    if i == 16467:
     print(dict[name],name)
     for j in range(len(a.features_list)):
         f = a.features_list[j]
         if 'gene' in f.qualifiers.keys() and f.type == 'gene':
             if f.qualifiers[f.type][0]==name:
                 print(f,f.location.start.position,f.location.end.position)

"""
#a.histogram(a.list_of_gene_no_cds)





"""

# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.


    print('Record number: {}\n============='.format(i))

features_list = record_gb.features


gensLength={}

def geneLength(featureslist):
 count=0
 for i in range(len(featureslist)):
   f = featureslist[i]
   if 'gene' in featureslist[i].qualifiers.keys() and f.type=='gene' :
    start_s = featureslist[i].location.start.position
    end_s = featureslist[i].location.end.position
    genName = featureslist[i].qualifiers[f.type][0]
    gensLength[genName] = end_s - start_s
 print('gensLength', gensLength)
 #return gensLength

geneLength(features_list)
def sortGene():
 for i in range(len(features_list)):
   f = features_list[i]
   if f.type in genbank.keys():
     genbank[f.type] += 1
   else:
     genbank[f.type] = 1
 print(genbank)
sortGene()

def gene_cds():
 list_of_gene_cds = []
 list_of_gene_no_cds = []
 genecdsdict = {}
 for i in range(len(features_list)):
  j = i + 1
  if (j<len(features_list)):
    f = features_list[i]
    r = features_list[j]
    if (f.type == 'gene' and r.type == 'CDS'):
      list_of_gene_cds.append(f)
    elif(f.type == 'gene'):
      list_of_gene_no_cds.append(f)
 gensLength(list_of_gene_cds)
 print(list_of_gene_cds)
 print(len(list_of_gene_cds),len(list_of_gene_no_cds))
gene_cds()

#gene_cds_statistics={}
def gene_cds_statistics():
  genecdsdict={}
  print(gensLength(list_of_gene_cds))
  #genecdsdict=gensLength(list_of_gene_cds)
  # print(genecdsdict)
  print(len(list_of_gene_cds), len(list_of_gene_no_cds))

#gene_cds_statistics()
  # gene_cds_statistics['max_value'] = max(genecdsdict.values())
  # gene_cds_statistics['min_value'] = min(genecdsdict.values())
  # gene_cds_statistics['mean_value'] = sum(genecdsdict.values()) / len(genecdsdict)
#   listr = []
#   # appending all the values in the list
#   for value in genecdsdict.values():
#     listr.append(value)
#   std = np.std(listr)
#   gene_cds_statistics['deviation'] = std
#   print(gene_cds_statistics)
# gene_cds()
# gene_cds_statistics()
# # for i in range(len(list_of_gene_cds)):
# #   f = list_of_gene_cds[i]
# #   if 'gene' in list_of_gene_cds[i].qualifiers.keys() and f.type == 'gene':
# #     start_s = list_of_gene_cds[i].location.start.position
# #     end_s = list_of_gene_cds[i].location.end.position
# #     genName = list_of_gene_cds[i].qualifiers[f.type][0]
# #     genecdsdict[genName] = end_s - start_s
#
#
#
# #
# #
# # gennocdsdictionary={}
# # gene_nocds_statistics={}
# # for i in range(len(list_of_gene_no_cds)):
# #   f = list_of_gene_cds[i]
# #   if 'gene' in list_of_gene_no_cds[i].qualifiers.keys() and f.type == 'gene':
# #     start_s = list_of_gene_no_cds[i].location.start.position
# #     end_s = list_of_gene_no_cds[i].location.end.position
# #     genName = list_of_gene_no_cds[i].qualifiers[f.type][0]
# #     gennocdsdictionary[genName] = end_s - start_s
# # print(gennocdsdictionary)
# #
# # gene_nocds_statistics['max_value'] =max(gennocdsdictionary.values())
# # gene_nocds_statistics['min_value'] =min(gennocdsdictionary.values())
# # gene_nocds_statistics['mean_value'] =sum(gennocdsdictionary.values()) / len(gennocdsdictionary)
# # listr = []
# # # appending all the values in the list
# # for value in gennocdsdictionary.values():
# #   listr.append(value)
# # # calculating standard deviation using np.std
# # std = np.std(listr)
# # gene_nocds_statistics['deviation'] =std
# # print(gene_nocds_statistics)
"""
