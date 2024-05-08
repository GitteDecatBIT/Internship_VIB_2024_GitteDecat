#!/usr/bin/env python

import collections
from biocypher._logger import logger
from template_package.adapters.CROssBAR_IRefIndex_pypath_url import url as irefindex_url
import pypath.share.curl as curl
import re
import pandas as pd


logger.info("Running input script for IRefIndex")

def irefindex_interactions():

    IRefIndexPhysicalInteraction = collections.namedtuple(
        'IRefIndexPhysicalInteraction',
        (
            'partner_a',
            'partner_b',
            'pmid', 
            'method',
            'organism', 
        ),
    )

    logger.info("Downloading interations from IRefIndex")
    interactions = []
    url = irefindex_url.get("irefindex").get("url")
    #'https://storage.googleapis.com/irefindex-data/archive/release_20.0/psi_mitab/MITAB2.6/7227.mitab.08-28-2023.txt.zip'
    c = curl.Curl(url, silent = False, large = True, slow = True)
    f = next(iter(c.result.values()))
    nul = f.readline()

    for l in f:
        l = l.split('\t')
        
        # PARTNER_A : finalReference A 
        input_partner_a = l[38]
        parts = input_partner_a.split(":") # Splitting the string at ":"
        if len(parts) > 1:
            partner_a = parts[1] # Selecting the second part after ":"
        else:
            partner_a = ""
        #print("partner_a: {}".format(partner_a))

        # PARTNER_B: FinalReference B 
        input_partner_b = l[39] 
        parts = input_partner_b.split(":")  # Splitting the string at ":"        
        if len(parts) > 1:  
            partner_b = parts[1]    # Selecting the second part after ":
        else:   
            partner_b = ""  
        #print("partner_b: {}".format(partner_b))    

        # PUBMED_ID
        input_pmid = l[8]
        pattern_pmid= r'\d+'
        pmid = re.findall(pattern_pmid, input_pmid)
        #print("pmid: {}".format(pmid))

        # METHOD
        input_method= l[6]
        pattern_method= r'\((.*?)\)'
        match_method= re.search(pattern_method, input_method)
        if match_method:
            method = match_method.group(1) # Extracting the text between brackets
        else:
            method = ""
        #print("method: {}".format(method))

        # ORGANISM
        input_organism= l[10]
        pattern_organism= r'taxid:(\d+)'
        match_organism= re.search(pattern_organism, input_organism)
        if match_organism:
            organism = match_organism.group(1)
        else:
            organism = ""
        #print("organism:{}".format(organism))
        
    
        interactions.append(
            IRefIndexPhysicalInteraction(
                partner_a = partner_a,
                partner_b = partner_b,
                pmid= pmid, 
                method= method, 
                organism= organism,  
            )
        
        )

    logger.info("Getting information for the IRefIndex database: partner_a, partner_b, pumed_ids, method and organism")



def irefindex_species() -> dict[int,str]:
    organism = {}

    url = irefindex_url.get("irefindex").get("url")
    #'https://storage.googleapis.com/irefindex-data/archive/release_20.0/psi_mitab/MITAB2.6/7227.mitab.08-28-2023.txt.zip'
    c = curl.Curl(url, silent = False, large = True, slow = True)
    f = next(iter(c.result.values()))
    nul = f.readline()

    for l in f:
            l = l.split('\t')
            # ORGANISM
            input_organism= l[10]
            pattern_organism= r'taxid:(\d+)'
            match_organism= re.search(pattern_organism, input_organism)
            if match_organism:
                organism = match_organism.group(1) 
            else:
                organism = ""

            return organism 


def irefindex_partner_a():
    partner_a= {}
    url = irefindex_url.get("irefindex").get("url")
    c = curl.Curl(url, silent = False, large = True, slow = True)
    f = next(iter(c.result.values()))
    nul = f.readline()


    for l in f:
        l = l.split('\t')
        # PARTNER_A : finalReference A 
        input_partner_a = l[38]
        parts = input_partner_a.split(":") # Splitting the string at ":"
        if len(parts) > 1:
            partner_a = parts[1] # Selecting the second part after ":"
        else:
            partner_a = ""
        #print("partner_a: {}".format(partner_a))

        return partner_a


def irefindex_partner_b():
    partner_b= {}
    url = irefindex_url.get("irefindex").get("url")
    c = curl.Curl(url, silent = False, large = True, slow = True)
    f = next(iter(c.result.values()))
    nul = f.readline()

    for l in f:
        l = l.split('\t')
        # PARTNER_B: FinalReference B 
        input_partner_b = l[39] 
        parts = input_partner_b.split(":")  # Splitting the string at ":"        
        if len(parts) > 1:  
            partner_b = parts[1]    # Selecting the second part after ":
        else:   
            partner_b = ""  
        #print("partner_b: {}".format(partner_b))    
            return partner_b
  
def irefindex_pmids():
    pmid= {}
    url = irefindex_url.get("irefindex").get("url")
    c = curl.Curl(url, silent = False, large = True, slow = True)
    f = next(iter(c.result.values()))
    nul = f.readline()

    for l in f:
        l = l.split('\t')

        # PUBMED_ID
        input_pmid = l[8]
        pattern_pmid= r'\d+'
        pmid = re.findall(pattern_pmid, input_pmid)
        #print("pmid: {}".format(pmid))

        return pmid


def irefindex_method():
    method = {}
    url = irefindex_url.get("irefindex").get("url")
    c = curl.Curl(url, silent = False, large = True, slow = True)
    f = next(iter(c.result.values()))
    nul = f.readline()

    for l in f:
        l = l.split('\t')
        # METHOD
        input_method= l[6]
        pattern_method= r'\((.*?)\)'
        match_method= re.search(pattern_method, input_method)
        if match_method:
            method = match_method.group(1) # Extracting the text between brackets
        else:
            method = ""
        #print("method: {}".format(method))
        return method 

