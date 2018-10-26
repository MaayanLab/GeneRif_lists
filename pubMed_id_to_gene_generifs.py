#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 11:04:14 2018

@author: megan wojciechowicz, maayan lab 

"""

from collections import Counter

gene_rif_file_names = ['Saccharomyces_cerevisiae_generifs.tsv']
pub_med_lists = ['Saccharomyces_cerevisiae_GO_biological_process_pubMed_Ids.txt']
gene_mapping_file = ['Saccharomyces_cerevisiae.gene_info']


for i, gene_rif_file in enumerate(gene_rif_file_names):
    outputfile = open(pub_med_lists[i].split('_pubMed_Ids')[0]+'_generif_gene_lists.tsv','w')
    
    # parse gene rif file and create dictionary with key:pubmed id --> values:list of gene ids
    gene_rif_dict = {}
    with open(gene_rif_file) as infile:
        for line in infile:
           line = line.strip('\n')
           line = line.split('\t')  
           gene_id = line[1]
           pubmed_id = line[2]
           if pubmed_id not in gene_rif_dict:
               gene_rif_dict[pubmed_id] = [gene_id] 
           else:
               gene_rif_dict[pubmed_id].append(gene_id)
    
    # parse mapping file to create dictionary that converts gene ids to gene symbols 
    gene_mapping = {}
    with open(gene_mapping_file[i]) as infile1:
        infile1.readline() # skip header
        for line in infile1:
            line = line.strip('\n')
            line = line.split('\t')
            gene_id = line[1]
            symbol = line[2]
            gene_mapping[gene_id]=symbol

            
    # parse pubMed query file and convert pubMed ids associated with each term --> gene ids 
    with open(pub_med_lists[i]) as infile2:
        for line in infile2:
            line = line.strip('\n')
            line = line.split('\t')
            term = line[0]
            p_ids = line[1:]
            genes = []
            for p_id in p_ids:
                try:
                    for x in gene_rif_dict[p_id]:
                        genes.append(x)
                except: 
                    pass
    
            # only keep lists with over 5 genes, order genes by count
            # if over 100 genes in a list, keep the top 100 only 
            if len(genes) >= 5:
                gene_count_dict = Counter(genes)
                sorted_genes = []
                for key, value in sorted(gene_count_dict.items(), key=lambda item: (item[1],item[0])):
                    sorted_genes.append(key)
                # find top 100 genes if list of genes is >100
                if len(sorted_genes)>100:
                    print('over 100 genes')
                    top_100_genes = sorted_genes[0:99]
                    count_100th_gene = gene_count_dict[sorted_genes[99]]
                    for key, value in gene_count_dict.items():
                        if key not in top_100_genes:
                            if value == count_100th_gene:
                                top_100_genes.append(key)
                    sorted_genes = top_100_genes
                
                # convert gene ids to gene symbols
                sorted_gene_symbols = []
                for x in sorted_genes:
                    try:
                        gene_symbol = gene_mapping[x]
                        if gene_symbol != 'NEWENTRY':
                            sorted_gene_symbols.append(gene_symbol)
                        else:
                            sorted_gene_symbols.append('ENTREZID:'+x)           
                    except:
                        sorted_gene_symbols.append('ENTREZID:'+x)
                        
                print('gene set length:'+str(len(sorted_gene_symbols ))+'\n')   
                outputfile.write(term +'\t'+ '\t'.join(sorted_gene_symbols)+'\n')
    outputfile.close()
                                


