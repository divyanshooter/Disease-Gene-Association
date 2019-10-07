# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 17:12:33 2019

@author: Chirag Garg and Divyanshu Chaturvedi
"""
import pandas as pd
import numpy as np
import heapq
import pprint

pp=pprint.PrettyPrinter(indent=4) 
# =============================================================================
# ppi ---->  protein-protein interaction
# =============================================================================
ppi=pd.read_csv("9606.protein.links.v11.0.csv")
gene_disease_da=pd.read_csv('gene-disease-ass.csv')
protein_data=pd.read_csv("9606.protein.info.v11.0.csv")

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

network_of_proteinprotein={}
for i in range(len(ppi)):
    protein1=ppi[i,0]
    protein2=ppi[i,1]
    score=ppi[i,2]
    if score<900:
        continue
    if protein1 in network_of_proteinprotein:
        network_of_proteinprotein[protein1].append([1/score,protein2])                                          #calculating inversion score and appending if node is present
    else:
        network_of_proteinprotein[protein1]=[[1/score,protein2]]                                                #calcualting inversion score and creating new node with edge weight as inversion score
pp.pprint(network_of_proteinprotein)
def genesymbol_disease():
    genesymbol_disease_mapping={}
    list_of_gene_symbol_from_totaldata = gene_disease[:,1:5]
    print(list_of_gene_symbol_from_totaldata)
    for gene_symbol in list_of_gene_symbol_from_totaldata:
        if gene_symbol[0] not in genesymbol_disease_mapping:
            genesymbol_disease_mapping[gene_symbol[0]] = [gene_symbol[3]]
        else:
            genesymbol_disease_mapping[gene_symbol[0]].append(gene_symbol[3])
# =============================================================================
#     pp.pprint(genesymbol_disease_mapping)
# =============================================================================
    return genesymbol_disease_mapping

def query_protein_disease(query_protein='9606.ENSP00000003084'):
    gene_symbol_of_query_protein = protein_gene_mapping[query_protein]
    genesymbol_disease_data = genesymbol_disease()
    disease_of_query_protein = genesymbol_disease_data[gene_symbol_of_query_protein]
    return gene_symbol_of_query_protein, disease_of_query_protein

def KNN(query_protein,k):  # k ----> No. of neighbouring proteins
    neighbouring_proteins=[]
    ppi=network_of_proteinprotein[query_protein]
    heapq.heapify(ppi)
    for i in range(k):
        neighbouring_proteins.append(heapq.heappop(ppi))
    return neighbouring_proteins

def neighbouringprotein_gene_mapping(neighbouring_proteins):
    gene_symbol=[]
    for protein in neighbouring_proteins:
        gene_symbol.append(protein_gene_mapping[protein[1]])
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

def score(query_protein_gene,a,neighbouringprotein_genes,neighbouringprotein_diseases_mapping,D,k):
    list_of_score_associations = []
    for disease in D:
        score_disease=0
        for j in neighbouringprotein_diseases_mapping.keys():
            if disease in neighbouringprotein_diseases_mapping[j]:
                score_disease+=jaccard_index(a,neighbouringprotein_diseases_mapping[j])
        
        score_disease=(100*score_disease)/k
        list_of_score_associations.append([query_protein_gene,disease,score_disease])
    return list_of_score_associations
query_protein = '9606.ENSP00000003084'
k=6
neighbouring_proteins = KNN(query_protein,k)
neighbouringprotein_genes = neighbouringprotein_gene_mapping(neighbouring_proteins)
neighbouringprotein_diseases_mapping, D = gene_disease_mapping(neighbouringprotein_genes)
query_protein_gene,a=query_protein_disease()
pp.pprint(score(query_protein_gene,a,neighbouringprotein_genes,neighbouringprotein_diseases_mapping,D,k))