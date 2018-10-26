#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 16 15:28:09 2018

@author: megan wojciechowicz, maayan lab 
"""
import urllib
import xml.etree.ElementTree as ET
import time 

# Searches pubMed and extracts pubMed ids associated with query 
# query = 'organims name' AND 'term'
# fucntion takes an organism and list of terms as input  
def get_pubmed_ids(organism_name, terms_list):
    # PubMed API URL
    all_pubmed_ids = []
    count = 0
    for term in terms_list:
        time.sleep(0.5)
        count = count + 1
        print('term:' + str(count))
        term_ids = []
        pubmed_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term='+organism_name.replace(' ','+')+'+AND+'+term.replace(' ','+')
        pubmed_data = ET.fromstring(urllib.request.urlopen(pubmed_url).read())
        for pubMed_id in pubmed_data.iter(tag='Id'):
            try:
                term_ids.append(str(pubMed_id.text))
            except:
                pass
        all_pubmed_ids.append(term_ids)
    return(all_pubmed_ids)     
