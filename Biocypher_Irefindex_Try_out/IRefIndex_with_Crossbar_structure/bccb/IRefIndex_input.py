#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2022
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: Dénes Türei (turei.denes@gmail.com)
#           Nicolàs Palacio
#           Sebastian Lobentanzer
#           Erva Ulusoy
#           Olga Ivanova
#           Ahmet Rifaioglu
#           Erva Ulusoy
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

from typing import List, Optional
from numbers import Number

import collections

from IRefIndex_pypath_url import url as irefindex_url
#from pypath_url import url 
import pypath.share.curl as curl

import re
import pandas as pd

print("running")


def irefindex_interactions():

    """
    Downloads and processes Physical multi-validated BioGRID interactions.
    Keeps only the "low throughput" interactions.
    Returns list of interactions.

    Args:
        organism: NCBI Taxonomy ID of organism.
        htp_limit: Exclude interactions only from references cited at more
        than this number of interactions.
    """

    IRefIndexPhysicalInteraction = collections.namedtuple(
        'IRefIndexPhysicalInteraction',
        (
            'partner_a',
            'partner_b',
            'pmid',
            'method',
            'organism'
        ),
    )

    #organism = str(organism)
    interactions = []
    refc = []
    url = irefindex_url.get("irefindex").get("url")
    #'https://storage.googleapis.com/irefindex-data/archive/release_20.0/psi_mitab/MITAB2.6/7227.mitab.08-28-2023.txt.zip'
    c = curl.Curl(url, silent = False, large = True, slow = True)
    f = next(iter(c.result.values()))
    nul = f.readline()

    for l in f:
        l = l.split('\t')

        """
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

        #print(organism)
        """

        #PARTNER_A
        input_partner_a = {
            'FinalReferenceA': [
                l[38]    
            ]
        }

        # Create DataFrame
        df_partner_a = pd.DataFrame(input_partner_a)

        # Split the 'ID' column into two columns
        df_partner_a[['Database', 'Accession']] = df_partner_a['FinalReferenceA'].str.split(':', expand=True)

        # Display the DataFrame
        print(df_partner_a)



        #PARTNER_B
        input_partner_b = {
            'FinalReferenceB': [
                l[39]    
            ]
        }

        # Create DataFrame
        df_partner_b = pd.DataFrame(input_partner_b)

        # Split the 'ID' column into two columns
        df_partner_b[['Database', 'Accession']] = df_partner_b['FinalReferenceB'].str.split(':', expand=True)

        # Display the DataFrame
        print(df_partner_b)


        #PMID_ID
        input_pmid = {
            'pmids': [
                l[8]    
            ]
        }

        # Create DataFrame
        df_pmid = pd.DataFrame(input_pmid)

        # Split the 'ID' column into two columns
        df_pmid[['Database', 'Accession']] = df_pmid['pmids'].str.split(':', expand=True)

        # Display the DataFrame
        print(df_pmid)

        
        #METHOD
        input_method = {
            'method': [
                l[6]    
            ]
        }

        # Create DataFrame
        df_method = pd.DataFrame(input_method)        

        # Extract text within brackets
        df_method['method'] = df_method['method'].str.extract(r'\((.*?)\)')         

        # Remove the extracted text from the original column
        df_method['method'] = df_method['method'].str.replace(r'\(.*?\)', '')          

        # Display the DataFrame
        print(df_method)


        # ORGANISM
        input_organism = {
            'tax_id': [
                l[10]    
            ]
        }
        # Create DataFrame
        df_organism = pd.DataFrame(input_organism)

        # _organismExtract text within parentheses
        df_organism['Species'] = df_organism['tax_id'].str.extract(r'\((.*?)\)')

        # _organismRemove the extracted text from the original column
        df_organism['tax_id'] = df_organism['tax_id'].str.replace(r'\(.*?\)', '')

        # Display the DataFrame
        print(df_organism)

        
        interactions.append(
            IRefIndexPhysicalInteraction(
                partner_a = df_partner_a.Accession,
                partner_b = df_partner_b,
                pmid= df_pmid, 
                method= df_method, 
                organism= df_organism,  
            )
        
        )
    
        refc.extend(pmid)

    refc = collections.Counter(refc)

    #print(interactions)

# Call the irefindex_interactions function
irefindex_interactions()