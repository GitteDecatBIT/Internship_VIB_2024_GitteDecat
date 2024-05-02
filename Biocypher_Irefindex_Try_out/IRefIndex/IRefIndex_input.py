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

        interactions.append(
            IRefIndexPhysicalInteraction(
                partner_a = partner_a,
                partner_b = partner_b,
                pmid= pmid, 
                method= method, 
                organism= organism,  
            )
        
        )
    
        refc.extend(pmid)

    refc = collections.Counter(refc)

    #print(interactions)

# Call the irefindex_interactions function
#irefindex_interactions()