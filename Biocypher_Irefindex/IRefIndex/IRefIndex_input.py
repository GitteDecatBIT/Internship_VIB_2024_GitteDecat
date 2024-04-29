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

from IRefIndex.pypath_url import urls
import pypath.resources.urls as urls
import pypath.share.curl as curl

import re



def irefindex_interactions(
        organism: int = 7227,
        htp_limit: Optional[Number] = 1,
        ltp: bool = True,
    ) -> List[tuple]:



    
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
        ),
    )

    organism = str(organism)
    interactions = []
    refc = []
    url = urls.urls['inrefindex']['all']
    c = curl.Curl(url, silent = False, large = True, slow = True)
    f = next(iter(c.result.values()))
    nul = f.readline()

    for l in f:

        l = l.split('\t')

        if len(l) > 17:

                interactions.append(
                    IRefIndexPhysicalInteraction(
                        partner_a = l[38],
                        partner_b = l[39],
                        pmid = l[8],
                    )
                )
                refc.append(l[8])

    refc = collections.Counter(refc)

    if htp_limit is not None:

        interactions = [i for i in interactions if refc[i[2]] <= htp_limit]

    return interactions


def irefindex_all_interactions(
        organism: int = 7227,
        htp_limit: Optional[Number] = 1,
        ltp: bool = True,
    ) -> List[tuple]:
    """
    Downloads and processes all BioGRID interactions.
    Keeps only the "low throughput" interactions.
    Returns list of interactions.

    Args:
        organism: NCBI Taxonomy ID of organism.
        htp_limit: Exclude interactions only from references cited at
            more than this number of interactions.
    """


    IRefIndexInteraction = collections.namedtuple(
        'IRefIndexInteraction',
        (
            'partner_a',
            'partner_b',
            'pmid',
            'method',
        ),
    )

    organism = str(organism)
    interactions = []
    refc = []
    mv_dict=collections.defaultdict(list)

    for i in irefindex_all_interactions(organism, htp_limit, ltp):

        mv_dict[i.partner_a].append(i.partner_b)

    url = urls.urls['irefindex']['all']
    c = curl.Curl(url, silent = False, large = True, slow = True)
    f = next(iter(c.result.values()))
    nul = f.readline()

    for l in f:

        l = l.split('\t')


        interactions.append(
            IRefIndexInteraction(
            partner_a = l[0],
            partner_b = l[1],
            pmid = l[8],
            method = l[6],
                       
                    )
                )

        refc.append(l[14])

    refc = collections.Counter(refc)

    if htp_limit is not None:

        interactions = [i for i in interactions if refc[i[2]] <= htp_limit]

    return interactions