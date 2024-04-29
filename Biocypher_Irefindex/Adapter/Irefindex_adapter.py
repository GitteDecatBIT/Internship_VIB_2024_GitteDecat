import os
import sys
import pandas as pd
import numpy as np
import time
import collections
from typing import Union

from pathlib import Path
from time import time


from pypath.inputs import biogrid, uniprot
from pypath.share import curl, settings

from tqdm import tqdm # progress bar

from biocypher._logger import logger
from pypath.resources import urls
from contextlib import ExitStack

from bioregistry import normalize_curie

from enum import Enum

global adapter_name
adapter_name = "irefindex" # might be useful for future

class IRefIndexEdgeFields(Enum):
    SOURCE = "source"
    PUBMED_IDS = "pmid"
    METHOD = "method"
    #INTERACTION__TYPE = "interaction_type"
    #EXPERIMENTAL_SYSTEM = "experimental_system"
    # type of interaction (e.g., physical interaction, regulatory interaction), the experimental evidence supporting the interaction, or the confidence score of the interaction.


class IRefIndex:
    def __init__(self, 
                 output_dir = None, 
                 export_csvs = False, 
                 split_output = False, 
                 cache=False, 
                 debug=False, 
                 retries=6,
                 organism=7227, 
                 irefindex_fields: Union[None, list[IRefIndexEdgeFields]] = None, 
                 add_prefix = True, 
                 test_mode = False, 
                 aggregate_pubmed_ids: bool = True,
                 aggregate_methods: bool = True):
        """
        Downloads and processes BioGRID data

            Args:
                export_csvs: Flag for whether or not create csvs of outputs of databases
                split_csvs: whether or not to split output csv files to multiple parts
                cache: if True, it uses the cached version of the data, otherwise
                forces download.
                debug: if True, turns on debug mode in pypath.
                retries: number of retries in case of download error.
                organism: taxonomy id number of selected organism, if it is None, downloads all organism data.
                biogrid_fields: biogrid fields to be used in the graph.
                add_prefix: if True, add prefix to uniprot ids
                test_mode: if True, take small sample from data for testing
                aggregate_pubmed_ids: if True, collects all pubmed ids that belongs to same protein pair
                aggregate_methods = if True, collects all experiemental methods that belongs to same protein pair
                
        """
        
        self.export_csvs = export_csvs
        self.split_output = split_output
        self.swissprots = list(uniprot._all_uniprots("*", True))
        self.cache = cache
        self.debug = debug
        self.retries = retries        
        self.organism = organism
        self.irefindex_fields = irefindex_fields
        self.add_prefix = add_prefix
        self.test_mode = test_mode
        
        self.aggregate_dict = {IRefIndexEdgeFields.PUBMED_IDS.value:aggregate_pubmed_ids,
                              IRefIndexEdgeFields.METHOD.value:aggregate_methods}

        if export_csvs:
            Path(output_dir).mkdir(parents=True, exist_ok=True)
            self.output_dir = output_dir
    
    def export_dataframe(self, dataframe, data_label):
        Path(self.output_dir).mkdir(parents=True, exist_ok=True)

        # TODO: activate this block after adding n_rows_in_file setting to config file
        # if self.split_output:
        #     n_chunks = round(len(dataframe) / n_rows_in_file)
            # output_path =  os.path.join(self.output_dir, data_label)
            # for id, chunk in  enumerate(np.array_split(dataframe, n_chunks)):
            #     chunk.to_csv(os.path.join(output_path, f"crossbar_ppi_data_{data_label}_{id+1}.csv"), index=False)
        # else:
        output_path = os.path.join(self.output_dir, f"irefindex_data_{data_label}.csv")
        dataframe.to_csv(output_path, index=False)

        return output_path
