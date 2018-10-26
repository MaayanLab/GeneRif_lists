#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 15 11:38:36 2018

@author: megan wojciechowicz, maayan lab 

"""
import pronto
from pronto.relationship import Relationship
import eutils_query as EQ
import copy
import numpy as np
import pandas as pd

# This function takes a list of parent terms and finds corresponding children
# Terms with no children are kept in the returned list 
# Input: ontology object, list of parent terms 
# Returns: list of terms 
def find_children_terms(ont, parents):
    new_terms = []
    for parent in parents:
        children = []
        for term in ont:
            ID = term.id            
            try:
                isa = ont[ID].relations[Relationship('is_a')]
                for x in isa:
                    if x.name == parent:
                        children.append(term.name)
                        new_terms.append(term.name)
            except:
                pass
        if len(children)==0:
            new_terms.append(parent)       
    return(new_terms)

# This function loops through each term in a list of terms and checks if there are also parent terms in the list
# If a parent term exists, all pubMed ids from child term are added to the parent term 
# Input: ontology object, list of terms, 2D list of pubMed ids associated with each term 
# Returns: 2D list of pubMed ids 
def back_propogation(ont,terms, ids):
    new_ids = copy.deepcopy(ids)
    for idx_term,term in enumerate(terms):
        for o in ont:
            if term == o.name:
                parents = [x.name for x in o.rparents()]
                for parent in parents:
                    if parent in terms:
                        print('parent:'+parent+'\t'+'child:'+term)
                        idx_parent = terms.index(parent)
                        for x in ids[idx_term]:
                            new_ids[idx_parent].append(x)
    return(new_ids)                                     



 

############################# FOR GO ontology only #############################

# This function subsets data in association files from http://geneontology.org/gene-associations/
# Code copied from GITHUB: https://github.com/MaayanLab/EnrichrPythonScripts/blob/master/WormEnrichr/GO/process_GO_wrombase.ipynb by Zachary Flamholz
# Input: dataframe , GO class type ('P'->biological process,'C'->cellular component,'F'->molecular function)
# Returns: subsetted dataframe 
def clean_data (df, go_type): 
    matrix = np.chararray((df.shape[0], 17), itemsize=150, unicode=True)
    for i, row in enumerate(df.itertuples()):
        lst = row[1].split('\t')
        matrix[i, :] = lst
    df_clean = pd.DataFrame(data=matrix)
    df_clean.columns = ["DB", "DB gene ID", "Gene symbol", "Qualifier", "GO ID", "Reference", 
                        "Evidence code", "Evidence from", "GO class", "attribute", "Locus tag",
                        "gene/protein", "tax id", "date", "Assigned by", "additional information", "empty"]
    df_clean= df_clean[df_clean["GO class"] == go_type]
    df_clean= df_clean[df_clean["Evidence code"] != 'IEA']# remove any annotation assigned by electronic matching and with the NOT qualifier which is used to specify a gene is not associated with a term
    df_clean= df_clean[df_clean["Qualifier"] != 'NOT']
    return(df_clean)

# This function takes subsetted gene association dataframe and extracts all GO ids
# Maps GO ids to GO terms 
# Input: dataframe , dictionary with GO ids --> GO terms 
# Returns list of unique GO terms with GO ids    e.g. ['GO term1 (GO:id1)', 'GO term2 (GO:id2)']
def get_GO_terms(df, lookup_table):  
    process = []
    for x in df['GO ID']:
        try:
            name = lookup_table[x]
            process.append(name + ' ' +'('+x+')')
        except: 
            pass
    unique_terms = np.unique(process)
    return(unique_terms)

############################# ^ FOR GO ontology only ^ #############################





############################## Upload ontology files ###############################

# fly anatomy ontology file: http://www.obofoundry.org/ontology/fbbt.html 
fly_anatomy_ont = pronto.Ontology('/Users/maayanlab/Downloads/fbbt.obo.txt')
# fly phenotype ontology: http://www.obofoundry.org/ontology/dpo.html
fly_phenotype_ont = pronto.Ontology('/Users/maayanlab/Downloads/fbcv.obo.txt')

# worm anatomy ontology file: http://www.obofoundry.org/ontology/wbbt.html
worm_anatomy_ont = pronto.Ontology('/Users/maayanlab/Downloads/wbbt.owl')
# worm phenptype file: http://www.obofoundry.org/ontology/wbphenotype.html
worm_phenotype_ont = pronto.Ontology('/Users/maayanlab/Downloads/wbphenotype.obo.txt')

# zebrafish anatomy ontology file: http://www.obofoundry.org/ontology/zfa.html
zebrafish_anatomy_ont = pronto.Ontology('/Users/maayanlab/Downloads/zfa.obo.txt')
# zebrafish phenotype ontology file: https://zfin.org/downloads/gene_expression_phenotype.txt
zebrafish_phenotype_ont = pd.DataFrame.from_csv('/Users/maayanlab/Downloads/phenotype_fish_2018.10.19.txt',sep= '\t', header=1).reset_index()

# yeast cellular component ontology : http://www.geneontology.org/ontology/subsets/goslim_yeast.obo
yeast_anatomy_ont = pronto.Ontology('/Users/maayanlab/Downloads/goslim_yeast.obo.txt')
# yeast phenotype ontology file: https://www.yeastgenome.org/search?category=download&page=0&topic=Genotype+and+phenotype&year=2017
yeast_phenotype_ont = pd.DataFrame.from_csv('phenotype_data.20170114.tab',sep= '\t', header=None).reset_index()

# GO ontology: http://snapshot.geneontology.org/ontology/go.obo
go_ont = pronto.Ontology('/Users/maayanlab/Downloads/go-basic.obo')
go_lookup = {}
for term in go_ont:
    go_lookup[term.id]=term.name
    
############################# ^ Upload ontology files ^ ##############################
    





############################## Parent terms of interest ##############################
# parent terms of interest
fly_anatomy_parents=[
                    'tagma',
                     'egg',
                     'embryo',
                     'organism',
                     'appendage',
                     'organ system',
                     'adipose system',
                     'circulatory system',
                     'digestive system',
                     'excretory system',
                     'integumentary system',
                     'muscle system',
                     'nervous system',
                     'reproductive system',
                     'tracheal system',
                     'portion of tissue',
                     'multi-tissue structure',
                     'sense organ',
                     'cell',             
                     'acellular anatomical structure',
                     'developing material anatomical entity'
                     ]         
worm_anatomy_parents=[
                    'axis',
                    'body region',
                    'cell component',
                    'dorsal nerve cord',
                    'extracellular component',
                    'ganglion',
                    'gubernaculum',
                    'lateral nerve cord',
                    'nerve ring',
                    'organ',
                    'organism',
                    'pharyngeal nerve',
                    'pharyngeal nerve process',
                    'pharyngeal segment',
                    'pharyngeal-intestinal valve',
                    'somatic gonad ',
                    'stoma',
                    'vas deferens',
                    'vas deferens valve region',
                    'ventral nerve cord'
                    ]
fly_phenotype_parents=[
                        'fertile',
                        'phenotype',
                        'viable',
                        'wild-type']
worm_phenotype_parents=[
                        'body region pigmentation variant',
                         'cell development variant',
                         'cell morphology variant',
                         'cell physiology variant',
                         'cell pigmentation variant',
                         'electrophysiology variant',
                         'organ system development variant',
                         'organ system morphology variant',
                         'organ system physiology variant',
                         'organ system pigmentation variant',
                         'organism behavior variant',
                         'organism development variant',
                         'organism morphology variant',
                         'organism physiology variant',
                         'organism pigmentation variant',
                         'pericellular component development variant',
                         'pericellular component morphology variant',
                         'pericellular component physiology variant',
                         'population fitness variant'
                         ]

############################# ^ Parent terms of interest ^ ############################







######################## Get terms to search PubMed ####################################
      
# fly terms
fly_anatomy_terms = [x for x in np.unique(find_children_terms(fly_anatomy_ont, fly_anatomy_parents))]
fly_phenotype_terms = [x for x in np.unique(find_children_terms(fly_phenotype_ont, fly_phenotype_parents))]
# fly gene and GO:terms --> GO Term processes names associated with fly
# file: http://geneontology.org/gene-associations/gene_association.fb.gz
fly_gene_associtation_biological_process = clean_data(pd.read_csv('/Users/maayanlab/Downloads/fb.gaf', sep = '`', skiprows=16, header=None),'P')[['Gene symbol', 'GO ID']]
fly_biological_process_terms = get_GO_terms(fly_gene_associtation_biological_process, go_lookup)


# worm terms
worm_anatomy_terms = find_children_terms(worm_anatomy_ont, worm_anatomy_parents)
worm_phenotype_terms = [x for x in np.unique(find_children_terms(worm_phenotype_ont, worm_phenotype_parents))]
# worm gene and GO:terms --> GO Term processes names associated with worm
# file: http://geneontology.org/gene-associations/gene_association.wb.gz
worm_gene_associtation_biological_process = clean_data(pd.read_csv('/Users/maayanlab/Downloads/wb.gaf', sep = '~', skiprows=14, header=None),'P')[['Gene symbol', 'GO ID']]
worm_biological_process_terms = get_GO_terms(worm_gene_associtation_biological_process, go_lookup)


# zebrafish terms
zebrafish_anatomy_terms=[]
for z in zebrafish_anatomy_ont:
    if z.name == 'zebrafish anatomical entity':
       zebrafish_anatomy_terms = [x.name for x in z.rchildren(level=3)]
zebrafish_phenotype_terms = []
for i,x in enumerate(zebrafish_phenotype_ont['Affected Structure or Process 1 superterm Name']):
    zebrafish_phenotype_terms.append(x)
zebrafish_phenotype_terms = np.unique(zebrafish_phenotype_terms)
# zebrafish gene and GO:terms --> GO Term processes names associated with zebrafish
# file: http://geneontology.org/gene-associations/gene_association.zfin.gz
zebrafish_gene_associtation_biological_process = clean_data(pd.read_csv('/Users/maayanlab/Downloads/zfin.gaf', sep = '~', skiprows=16, header=None),'P')[['Gene symbol', 'GO ID']]
zebrafish_biological_process_terms = get_GO_terms(zebrafish_gene_associtation_biological_process, go_lookup)


# yeast terms 
yeast_anatomy_terms=[]
for y in yeast_anatomy_ont:
    if y.name == 'cellular_component':
        yeast_anatomy_terms = [x.name for x in y.rchildren(level=-1)]
yeast_phenotype_terms = np.unique([x.split(':')[0] for x in yeast_phenotype_ont.loc[:,9]])
# yeast gene and GO:terms --> GO Term processes names associated with yeast
# file: http://geneontology.org/gene-associations/gene_association.sgd.gz
yeast_gene_associtation_biological_process = clean_data(pd.read_csv('/Users/maayanlab/Downloads/sgd.gaf', sep = '~', skiprows=18, header=None),'P')[['Gene symbol', 'GO ID']]
yeast_biological_process_terms = get_GO_terms(yeast_gene_associtation_biological_process, go_lookup)

#######################  ^ Get terms to search PubMed ^ ###################################





##############################  Get PubMed Results  #########################################


# get pubmed results, perform back propogation, and write to file 
def pubmed_results_to_file(organism,terms, type_terms,ontology_file, backpropogate):
    ids = EQ.get_pubmed_ids(organism,[x.split(' (')[0] for x in terms])# GO Terms in parentheses must be removed for term search  
    if backpropogate =='yes':
        ids = back_propogation(ontology_file,terms, ids)
    writefile = open(organism.replace(' ', '_')+'_'+type_terms+'_pubMed_Ids.txt',"w")
    for i,x in enumerate (terms):
        if len(ids[i])>0:
            ids2 = np.unique(ids[i])
            writefile.write(str(x).replace(' ','_')+'\t'+'\t'.join([str(k) for k in ids2])+'\n')
        

#pubmed_results_to_file('Drosophila melanogaster',fly_anatomy_terms,'anatomy',fly_anatomy_ont,'yes')
#pubmed_results_to_file('Drosophila melanogaster',fly_phenotype_terms,'phenotype',fly_phenotype_ont,'yes')
######pubmed_results_to_file('Drosophila melanogaster',fly_pathway_terms,'GO_biological_process',go_ont,'yes')
#pubmed_results_to_file('Drosophila melanogaster',fly_biological_process_terms,'GO_biological_process',go_ont,'yes')
            
#pubmed_results_to_file('Caenorhabditis elegans',worm_anatomy_terms,'anatomy',worm_anatomy_ont,'yes')
#pubmed_results_to_file('Caenorhabditis elegans',worm_phenotype_terms,'phenotype',worm_phenotype_ont,'yes')   
#pubmed_results_to_file('Caenorhabditis elegans',worm_biological_process_terms,'GO_biological_process',go_ont,'yes')  
            
#pubmed_results_to_file('Danio_rerio',zebrafish_anatomy_terms,'anatomy',zebrafish_anatomy_ont,'yes')
#pubmed_results_to_file('Danio_rerio',zebrafish_phenotype_terms,'phenotype',zebrafish_phenotype_ont,'no') # no back propogation with a simple list of terms (no ontology tree)
#pubmed_results_to_file('Danio_rerio',zebrafish_biological_process_terms,'GO_biological_process',go_ont,'yes')  
            
#pubmed_results_to_file('Saccharomyces_cerevisiae',yeast_anatomy_terms,'cellular_component',yeast_anatomy_ont,'yes')
#pubmed_results_to_file('Saccharomyces_cerevisiae',yeast_phenotype_terms,'phenotype',yeast_phenotype_ont,'no')  # no back propogation with a simple list of terms (no ontology tree)
#pubmed_results_to_file('Saccharomyces_cerevisiae',yeast_biological_process_terms,'GO_biological_process',go_ont,'yes')

############################## ^ Get PubMed Results ^ ######################################### 