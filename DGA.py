# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 00:35:25 2019

@author: Divyanshu
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 17:12:33 2019

@author: Chirag Garg and Divyanshu Chaturvedi
"""
import pandas as pd
import numpy as np
import heapq
import pprint
import math

pp=pprint.PrettyPrinter(indent=4) 
# =============================================================================
# ppi ---->  protein-protein interaction
# =============================================================================

ppi=pd.read_csv("9606.protein.links.v11.0.csv")
gene_disease_da=pd.read_csv('gene-disease-ass.csv')
protein_data=pd.read_csv("9606.protein.info.v11.0.csv")
protein_set=set()

ppi=ppi.iloc[:,:].values
protein_info=protein_data.iloc[:,:].values
gene_disease=gene_disease_da.iloc[:,:].values
protein_gene_data=protein_info[:,0:2]
protein_gene_mapping={}
for protein in protein_gene_data:
    protein_gene_mapping[protein[0]]=protein[1]
# =============================================================================
# pp.pprint(protein_gene_mapping)
# =============================================================================
def Sort(sub_li): 
    sub_li.sort(key = lambda x: x[2],reverse=True) 
    return sub_li 
network_of_proteinprotein={}
for i in range(len(ppi)):
    protein1=ppi[i,0]
    protein2=ppi[i,1]
    score=ppi[i,2]
   
    if score<900:
        continue
    if protein1 in network_of_proteinprotein:
        protein_set.add(protein1)
        network_of_proteinprotein[protein1].append([1/score,protein2])                                          #calculating inversion score and appending if node is present
    else:
        protein_set.add(protein1)
        network_of_proteinprotein[protein1]=[[1/score,protein2]]                                                #calcualting inversion score and creating new node with edge weight as inversion score
#pp.pprint(network_of_proteinprotein)
def genesymbol_disease():
    genesymbol_disease_mapping={}
    list_of_gene_symbol_from_totaldata = gene_disease[:,1:5]
  #  print(list_of_gene_symbol_from_totaldata)
    for gene_symbol in list_of_gene_symbol_from_totaldata:
        if gene_symbol[0] not in genesymbol_disease_mapping:
            genesymbol_disease_mapping[gene_symbol[0]] = [gene_symbol[3]]
        else:
            genesymbol_disease_mapping[gene_symbol[0]].append(gene_symbol[3])
# =============================================================================
#     pp.pprint(genesymbol_disease_mapping)
# =============================================================================
    return genesymbol_disease_mapping

def query_protein_disease(query_protein='9606.ENSP00000000412'):
    gene_symbol_of_query_protein = protein_gene_mapping[query_protein]
    genesymbol_disease_data = genesymbol_disease()
    idx=gene_disease_da.loc[gene_disease_da['geneSymbol'] == gene_symbol_of_query_protein ]
    if(len(idx)==0):
        return -1,-1
    disease_of_query_protein = genesymbol_disease_data[gene_symbol_of_query_protein]
    return gene_symbol_of_query_protein, disease_of_query_protein

def KNN(query_protein,k):  # k ----> No. of neighbouring proteins
    neighbouring_proteins=[]
    KNN_ppi=network_of_proteinprotein[query_protein]
    if(len(KNN_ppi)<k):
        return -1
    heapq.heapify(KNN_ppi)
    for i in range(k):
        neighbouring_proteins.append(heapq.heappop(KNN_ppi))
    return neighbouring_proteins

def neighbouringprotein_gene_mapping(neighbouring_proteins):
    gene_symbol=[]
    protein_to_removed=[]
    for protein in neighbouring_proteins:
        gene_symbol_protein=protein_gene_mapping[protein[1]]
        idx=gene_disease_da.loc[gene_disease_da['geneSymbol'] == gene_symbol_protein ]
        if(len(idx)==0):
           protein_to_removed.append(protein)
        else:
            gene_symbol.append(gene_symbol_protein)
    if(len(protein_to_removed)!=0):
        for i in protein_to_removed:
            neighbouring_proteins.remove(i)
# =============================================================================
#     print(gene_symbol)
# =============================================================================
    return gene_symbol
    
def gene_disease_mapping(neighbouringprotein_genes):
    genesymbol_disease_data = genesymbol_disease()
    d={}
    D=set()
    for i in neighbouringprotein_genes:
        if i in genesymbol_disease_data.keys():
            d[i]=genesymbol_disease_data[i]
            for j in d[i]:
                D.add(j)
    return d,D

def jaccard_index(gene1,gene2):
    gene1=set(gene1)
    gene2=set(gene2)
    length_of_intersection = len(gene1.intersection(gene2))
    length_of_gene1,length_of_gene2=len(gene1),len(gene2)
    return length_of_intersection/(length_of_gene1+length_of_gene2)
def simpson_index(gene1,gene2):
    gene1=set(gene1)
    gene2=set(gene2)
    length_of_intersection = len(gene1.intersection(gene2))
    length_of_gene=min(len(gene1),len(gene2))
    return length_of_intersection/length_of_gene

def geometric_index(gene1,gene2):
    gene1=set(gene1)
    gene2=set(gene2)
    length_of_intersection = pow(len(gene1.intersection(gene2)),2)
    length_of_gene1,length_of_gene2=len(gene1),len(gene2)
    return length_of_intersection/(length_of_gene1*length_of_gene2)

def cosine_index(gene1,gene2):
    gene1=set(gene1)
    gene2=set(gene2)
    length_of_intersection = len(gene1.intersection(gene2))
    length_of_gene1,length_of_gene2=len(gene1),len(gene2)
    return length_of_intersection/math.sqrt(length_of_gene1*length_of_gene2)

def score(query_protein_gene,a,neighbouringprotein_genes,neighbouringprotein_diseases_mapping,D,k):
    list_of_score_associations_jaccard = []
    list_of_score_associations_simpson = []
    list_of_score_associations_geometric = []
    list_of_score_associations_cosine = []
    for disease in D:
        score_disease_jaccard=0
        score_disease_simpson=0
        score_disease_geometric=0
        score_disease_cosine=0
        for j in neighbouringprotein_diseases_mapping.keys():
            if disease in neighbouringprotein_diseases_mapping[j]:
                score_disease_jaccard+=jaccard_index(a,neighbouringprotein_diseases_mapping[j])
                score_disease_simpson+simpson_index(a,neighbouringprotein_diseases_mapping[j])
                score_disease_geometric+=geometric_index(a,neighbouringprotein_diseases_mapping[j])
                score_disease_cosine+=cosine_index(a,neighbouringprotein_diseases_mapping[j])
        score_disease_jaccard=(100*score_disease_jaccard)/k
        score_disease_simpson=(100*score_disease_simpson)/k
        score_disease_geometric=(100*score_disease_geometric)/k
        score_disease_cosine=(100*score_disease_cosine)/k
        list_of_score_associations_jaccard.append([query_protein_gene,disease,score_disease_jaccard])
        list_of_score_associations_simpson.append([query_protein_gene,disease,score_disease_simpson])
        list_of_score_associations_geometric.append([query_protein_gene,disease,score_disease_geometric])
        list_of_score_associations_cosine.append([query_protein_gene,disease,score_disease_cosine])
    return list_of_score_associations_jaccard,list_of_score_associations_simpson,list_of_score_associations_geometric,list_of_score_associations_cosine
k=6
protein_set_list=list(protein_set)
def DGA():
    for i in range(len(protein_set_list)):
        protein1=protein_set_list[i]
        query_protein = protein1
        print(query_protein)
        neighbouring_proteins = KNN(query_protein,k)
        if(neighbouring_proteins==-1):
            continue
        neighbouringprotein_genes = neighbouringprotein_gene_mapping(neighbouring_proteins)
        neighbouringprotein_diseases_mapping, D = gene_disease_mapping(neighbouringprotein_genes)
        query_protein_gene,a=query_protein_disease(query_protein)
        if(query_protein_gene==-1 and a==-1):
             print("No Disease exist for given Protein")
        else:
            list_of_score_associations_jaccard,list_of_score_associations_simpson,list_of_score_associations_geometric,list_of_score_associations_cosine=score(query_protein_gene,a,neighbouringprotein_genes,neighbouringprotein_diseases_mapping,D,k)
            print('Jaccquard Index')
            pp.pprint(Sort(list_of_score_associations_jaccard)[0:k])
            print('Simpson Index')
            pp.pprint(Sort(list_of_score_associations_simpson)[0:k])
            print('Geometric Index')
            pp.pprint(Sort(list_of_score_associations_geometric)[0:k])
            print('Cosine Index')
            pp.pprint(Sort(list_of_score_associations_cosine)[0:k])
DGA()

