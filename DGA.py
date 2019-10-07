import pandas as pd
import numpy as np
import heapq
ppi=pd.read_csv("9606.protein.links.v11.0.csv")
gene_disease_da=pd.read_csv('gene-disease-ass.csv')
protein_data=pd.read_csv("9606.protein.info.v11.0.csv")
ppi=ppi.iloc[:,:].values
protein_info=protein_data.iloc[:,:].values
gene_disease=gene_disease_da.iloc[:,:].values

print(protein_info[:,1])


network={}

for i in range(len(ppi)):
    u=ppi[i,0]
    v=ppi[i,1]
    score=ppi[i,2]
    if score<900:
        continue
    if u in network:
        network[u].append([1/score,v])                                          #calculating inversion score and appending if node is present
         
    else:
        network[u]=[[1/score,v]]                                                #calcualting inversion score and creating new node with edge weight as inversion score

def KNN(q,k):
    protein=[]
    ppi=network[q]
    heapq.heapify(ppi)

    
    for i in range(k):
        protein.append(heapq.heappop(ppi))
    return protein

k_nearest=KNN('9606.ENSP00000000233',6)

    
