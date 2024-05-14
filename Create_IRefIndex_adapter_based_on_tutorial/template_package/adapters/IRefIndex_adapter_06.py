#!/usr/bin/env python
import random
import string
from enum import Enum, auto
from itertools import chain
from typing import Optional
from biocypher._logger import logger


# extra from CROssBAR

from time import time

from pypath.share import curl, settings
import re
from typing import Union
#from create_knowledge_graph_IrefIndex import paths
import os 
import collections
import pandas as pd
import numpy as np
from pypath.inputs import uniprot
from bioregistry import normalize_curie
from tqdm import tqdm
from pathlib import Path
#############################################################################################################


logger.debug(f"Loading module {__name__}.")

class IRefIndexNodeType(Enum):
    PROTEIN= auto()

class IRefIndexEdgeFields(Enum):
    #PROTEIN= auto()
    #SOURCE = "source"
    PUBMED_IDS = "pmid"
    TAXON= "taxon"
    METHODS = "method"
    #PROTEIN_PROTEIN_INTERACTION = "protein_protein_interaction"
    #INTERACTION_TYPE = "interaction_type" # column [11]
    #INTERACTION_SOURCE = "interaction_source" # column [12]
    

    
class IRefIndexAdapter:
    """
    Example BioCypher adapter. Generates nodes and edges for creating a
    knowledge graph.

    Args:
        node_types: List of node types to include in the result.
        node_fields: List of node fields to include in the result.
        edge_types: List of edge types to include in the result.
        edge_fields: List of edge fields to include in the result.
    """
    def __init__(self, 
                 output_dir = None,
                 export_csvs = False, 
                 split_output = False, 
                 cache=False, 
                 debug=False, 
                 retries=6,
                 organism=9606, 
                 irefindex_fields: Union[None, list[IRefIndexEdgeFields]] = None, 
                 add_prefix = True, 
                 test_mode = False, 
                 aggregate_pubmed_ids: bool = True,
                 aggregate_methods: bool = True,
                 node_types: Union[None, list[IRefIndexNodeType]] = None,
                 node_fields:Union[None, list[IRefIndexEdgeFields]] = None,
                 edge_types: Union[None, list[IRefIndexEdgeFields]] = None,
                 edge_fields:Union[None, list[IRefIndexEdgeFields]] = None,
            
                 ):
        self.output_dir = output_dir
        self.export_csvs = export_csvs
        self.irefindex_fields = irefindex_fields
        self.swissprots = list(uniprot._all_uniprots("*", True))
       
        self.add_prefix = add_prefix

        self.aggregate_dict = {IRefIndexEdgeFields.PUBMED_IDS.value:aggregate_pubmed_ids,
                              IRefIndexEdgeFields.METHODS.value:aggregate_methods}
        

        self._set_types_and_fields(node_types, node_fields, edge_types, edge_fields)
         
        if export_csvs:
            Path(output_dir).mkdir(parents=True, exist_ok=True)
            self.output_dir = output_dir

    def export_dataframe(self,irefindex_df_unique):
        Path(self.output_dir).mkdir(parents=True, exist_ok=True)

        # TODO: activate this block after adding n_rows_in_file setting to config file
        # if self.split_output:
        #     n_chunks = round(len(dataframe) / n_rows_in_file)
            # output_path =  os.path.join(self.output_dir, data_label)
            # for id, chunk in  enumerate(np.array_split(dataframe, n_chunks)):
            #     chunk.to_csv(os.path.join(output_path, f"crossbar_ppi_data_{data_label}_{id+1}.csv"), index=False)
        # else:
        output_path = os.path.join(self.output_dir, f"irefindex_data.csv")
        irefindex_df_unique.to_csv(output_path, index=False)

        return output_path

# column in process 
# get nodes --> zie wat er in de nodes dan moet       
    def irefindex_process(self):

        # Extracting information from IRefIndex
        logger.info("EXECUTING THE ADAPTER FILE:")
        logger.info("Extracting information from IRefIndex data")
        
        IRefIndexInteraction = collections.namedtuple(
            'IRefIndexPhysicalInteraction',
            (
                'partner_a',
                'partner_b',
                'pmid', 
                'method',
                'taxon', 
            ),
        )

        logger.info("1) Created a collection 'IRefIndexInteractions' for the variables: ")

        logger.info("2) Searching for the path the data text file is located")
        # Path to the extracted data directory
        input_file_dir = '.cache/IRefIndex/7227.mitab.08-28-2023.txt.zip.unzip'
        logger.info("--> input file dir: {}" .format(input_file_dir))
        # Iterate over the files in the extracted directory and parse each file as needed
        for filename in os.listdir(input_file_dir):
            if filename.endswith(".txt"):
                inputfile_path = os.path.join(input_file_dir, filename)
        logger.info("3) Open/read the text file ")
        logger.info(" --> inputfile:{}".format(inputfile_path))

 
        logger.info("4)Exctracting information in certain columns from the IRefIndex database")

        interactions = []

        # Open the file
        with open(inputfile_path, 'r') as file:
            logger.info("--> Getting inforamtion for partner_a, partner_b, pubmed_id, method and taxon id")
            for line in file:
                # Split the line by tab character
                l = line.split('\t')
                
                # PARTNER_A: finalReference A 
                input_partner_a = l[38]
                parts = input_partner_a.split(":")
                partner_a = parts[1] if len(parts) > 1 else ""

                input_partner_a = l[38]
                parts = input_partner_a.split(":") # Splitting the string at ":"
                if len(parts) > 1:
                    partner_a = parts[1] # Selecting the second part after ":"
                else:
                    partner_a = ""

                # PARTNER_B: FinalReference B 
                input_partner_b = l[39]
                parts = input_partner_b.split(":")
                if len(parts) > 1:
                    partner_b = parts[1]  
                else:
                    partner_b= ""

                # PUBMED_ID
                input_pmid = l[8]
                # Split the string by '|'
                parts = input_pmid.split('|')
                # Take the last part of the split
                last_part = parts[-1]

                # Split again by '.'
                numbers = last_part.split(':')

                # Take the last part of this split, which is the number
                pmid = numbers[-1]

                #pattern_pmid = r'\d+'
                #pmid = re.findall(pattern_pmid, input_pmid)



                # METHOD
                input_method= l[6]
                pattern_method= r'\((.*?)\)'
                match_method= re.search(pattern_method, input_method)
                if match_method:
                    method = match_method.group(1) # Extracting the text between brackets
                else:
                    method = ""

                # ORGANISM/taxon 
                input_taxon= l[10]
                pattern_taxon= r'taxid:(\d+)'
                match_taxon= re.search(pattern_taxon, input_taxon)
                if match_taxon:
                    taxon = match_taxon.group(1)
                else:
                    taxon = ""

                interactions.append(
                    IRefIndexInteraction(
                        partner_a=partner_a,
                        partner_b=partner_b,
                        pmid=pmid,
                        method=method,
                        taxon=taxon,
                    )
                )
                
            logger.info("--> Succesfully extracted columns from the IRefIndex database!")
            self.irefindex_ints = interactions
            #logger.info(interactions)
            return interactions
        
    

    def get_nodes(self, rename_selected_fields: Union[None, list[str]] = None):
        """
        Processor function for irefindex data. It drops duplicate and reciprocal duplicate protein pairs and collects pubmed ids of duplicated pairs. In addition, it
        maps entries to uniprot ids using gene name and tax id information in the BioGRID data. Also, it filters protein pairs found in swissprot.
        
         Args:
            rename_selected_fields : List of new field names for selected fields. If not defined, default field names will be used.
        """
        logger.info("5) Generating nodes.")

        selected_fields = self.set_edge_fields()
            
        default_field_names = {"pmid":"pubmed_ids", "taxon":"taxon" ,"method":"method"}
        
        self.irefindex_field_new_names = {}
        
        if rename_selected_fields:
            if len(selected_fields) != len(rename_selected_fields):
                raise Exception("Length of selected_fields variable should be equal to length of rename_selected_fields variable")
            
            for field_old_name, field_new_name in list(zip(selected_fields, rename_selected_fields)):

                self.irefindex_field_new_names[field_old_name] = field_new_name
                #logger.info("field_old_name:{}".format(field_old_name))
                #logger.info("field_new_name:{}".format(field_new_name))
                # geeft deze niet in terminal

            self.irefindex_field_new_names["partner_a"] = "uniprot_a"
            self.irefindex_field_new_names["partner_b"] = "uniprot_b"
        else:
            for field_old_name in selected_fields:
                #logger.info("field_old_name:{}".format(field_old_name))
                self.irefindex_field_new_names[field_old_name] = default_field_names[field_old_name]
            
            self.irefindex_field_new_names["partner_a"] = "uniprot_a"
            self.irefindex_field_new_names["partner_b"] = "uniprot_b"
        
        logger.info(" --> Changing 'partner_a' and 'partner_b' to 'uniprot_a' and 'uniprot_b'")
        logger.info("6) Started processing IRefIndex data")
        t1 = time()
                         
        # create dataframe          
        irefindex_df = pd.DataFrame.from_records(self.irefindex_ints, columns=self.irefindex_ints[0]._fields)

        #logger.info("create irefindex_dataframe:{}".format(irefindex_df)) # --> partner_a |partner_b |pmid |method |taxon
        logger.info("--> Created an irefindex_dataframe")

        
        # add source database info
        irefindex_df["source"] = "IRefIndex"
        # filter selected fields
        irefindex_df = irefindex_df[list(self.irefindex_field_new_names.keys())]
        #logger.info("filter irefindex_dataframe:{}".format(irefindex_df)) # --> pmid| taxon| method| partner_a |partner_b
        logger.info(" --> Filtered the irefindex_dataframe")
        # rename columns
        irefindex_df.rename(columns=self.irefindex_field_new_names, inplace=True)
        #logger.info("rename irefindex_dataframe:{}".format(irefindex_df))
        logger.info(" --> Renamed the headers of the irefindex_dataframe to uniprot_a and uniprot_b")
        
        # drop rows if uniprot_a or uniprot_b is not a swiss-prot protein
        irefindex_df = irefindex_df[(irefindex_df["uniprot_a"].isin(self.swissprots)) & (irefindex_df["uniprot_b"].isin(self.swissprots))]
        irefindex_df.reset_index(drop=True, inplace=True)
        #logger.info("drop in irefindex_dataframe if not in swissprot :{}".format(irefindex_df)) # --> pubmed_ids | taxon | method |uniprot_a |uniprot_b
        logger.info(" --> Dropped the rows in irefindex_dataframe if uniprot_a or uniprot_b is not in swissprot")
        
        # drop duplicates if same a x b pair exists multiple times 
        # keep the first pair and collect pubmed ids of duplicated a x b pairs in that pair's pubmed id column
        # if a x b pair has same experimental system type with b x a pair, drop b x a pair
        irefindex_df_unique = irefindex_df.dropna(subset=["uniprot_a", "uniprot_b"]).reset_index(drop=True)        
        logger.info(" --> Created a unique irefindex_dataframe, removed duplicates ")
        
        
        def aggregate_fields(element):
            #logger.info("element:{}".format(element)) # element is pubmed id 
            element = "|".join([str(e) for e in set(element.dropna())])
            if not element:
                return np.nan
            else:
                return element
        logger.info("7) Takes an element representing the Pubmed_id + drops any missing values ")
        logger.info("8) Converts remaining values to a set to remove duplicates + joins them into a single string separated by '|'. This operation effectively aggregates multiple PubMed IDs into a single string.")

        if any(list(self.aggregate_dict.values())):
            agg_field_list = [k for k, v in self.aggregate_dict.items() if v]
            
            agg_dict = {}            
            for k, v in self.irefindex_field_new_names.items():
                if k in agg_field_list:
                    agg_dict[v] = aggregate_fields
                else:                
                    agg_dict[v] = "first"
        
        irefindex_df_unique = irefindex_df_unique.groupby(["uniprot_a", "uniprot_b"], sort=False, as_index=False).aggregate(agg_dict)
        #biogrid_df_unique["pubmed_id"].replace("", np.nan, inplace=True)
        logger.info("9) Group the irefindex_dataframe based in uniprot_a and uniprot_b")
        
        if "method" in self.irefindex_field_new_names.keys():            
            irefindex_df_unique = irefindex_df_unique[~irefindex_df_unique[["uniprot_a", "uniprot_b", self.irefindex_field_new_names["method"]]].apply(frozenset, axis=1).duplicated()].reset_index(drop=True)
        else:
            irefindex_df_unique = irefindex_df_unique[~irefindex_df_unique[["uniprot_a", "uniprot_b"]].apply(frozenset, axis=1).duplicated()].reset_index(drop=True)

        if self.export_csvs:
            irefindex_output_path = self.export_dataframe(irefindex_df_unique, "irefindex")
            logger.info(f'Final IRefIndex data is written: {irefindex_output_path}') 
             # UnboundLocalError: cannot access local variable 'irefindex_output_path' where it is not associated with a value


        self.final_irefindex_ints = irefindex_df_unique
        logger.info("FINAL DATAFRAME:{}".format(self.final_irefindex_ints)) # --> processed dataframe
        
        t2 = time()
        logger.info(f'IRefIndex data is processed in {round((t2-t1) / 60, 2)} mins')

        """
        Returns a generator of node tuples for node types specified in the
        adapter constructor.
        """


        
            
    def set_edge_fields(self) -> list:
        """
        Sets biogrid edge fields
        Returns:
            selected field list
        """
        if self.irefindex_fields is None:
            return [field.value for field in IRefIndexEdgeFields]
        else:
            return [field.value for field in self.irefindex_fields]

 
    
    def add_prefix_to_id(self, prefix="uniprot", identifier=None, sep=":") -> str:
        """
        Adds prefix to uniprot id
        """
        if self.add_prefix and identifier:
            return normalize_curie( prefix + sep + str(identifier))
        
        return identifier
       


        #https://github.com/biocypher/dependency-map/blob/main/dmb/adapter.py#L342 --> line 312 ( get nodes)
        
        
    

    def get_edges(self) -> list:
        """
        Get PPI edges from biogrid data
        """
        
        # create edge list
        edge_list = []
        
        for index, row in tqdm(self.final_irefindex_ints.iterrows(), total=self.final_irefindex_ints.shape[0]):
            _dict = row.to_dict()
            
            _source = self.add_prefix_to_id(identifier = _dict["uniprot_a"])
            _target = self.add_prefix_to_id(identifier = _dict["uniprot_b"])
            
            del _dict["uniprot_a"], _dict["uniprot_b"]
            
            _props = dict()
            for k, v in _dict.items():
                if str(v) != "nan":
                    if isinstance(v, str) and "|" in v:
                        # if column has multiple entries create list
                        _props[str(k).replace(" ","_").lower()] = v.split("|")
                    else:
                        _props[str(k).replace(" ","_").lower()] = v
           

            edge_list.append((None, _source, _target, "Interacts_With", _props))
            
        #logger.info(edge_list)
        return edge_list

    def _set_types_and_fields(self, node_types, node_fields, edge_types, edge_fields):
        if node_types:
            self.node_types = node_types
        else:
            self.node_types = [type for type in IRefIndexEdgeFields]

        if node_fields:
            self.node_fields = node_fields
        else:
            self.node_fields = [
                field
                for field in chain(
                    IRefIndexEdgeFields,
        
                )
            ]

        if edge_types:
            self.edge_types = edge_types
        else:
            self.edge_types = [type for type in IRefIndexEdgeFields]

        if edge_fields:
            self.edge_fields = edge_fields
        else:
            self.edge_fields = [
                field 
                for field in chain(
                IRefIndexEdgeFields,
                )
            ]



