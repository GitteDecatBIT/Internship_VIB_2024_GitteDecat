#!/usr/bin/env python
from biocypher._logger import logger
from enum import Enum, auto
from typing import Union
from pathlib import Path
from itertools import chain
import os 
import collections
import re
import pandas as pd
import numpy as np
from time import time
from typing import Optional
from bioregistry import normalize_curie
from tqdm import tqdm # progress bar

logger.debug(f"Loading module {__name__}.")

class IRefIndexNodeType(Enum):
    PROTEIN= auto()

class IRefIndexEdgeType(Enum):
    """
    Enum for the types of the protein adapter.
    """
    PROTEIN_PROTEIN_INTERACTION = "protein_protein_interaction"

class IRefIndexEdgeFields(Enum):
    PUBMED_IDS = "pmid"
    TAXON= "taxon"
    METHODS = "method"
    RELATIONSHIP_ID = "relationship_id"

class IRefIndexAdapter:
    """
    Adapter that creaes nodes and edges for creating a knowledge graph.

    Args:
        node_types: List of node types to include in the result.
        node_fields: List of node fields to include in the result.
        edge_types: List of edge types to include in the result.
        edge_fields: List of edge fields to include in the result.
    """
    def __init__(self,output_dir = None,irefindex_fields: Union[None, list[IRefIndexEdgeFields]] = None, 
                 add_prefix = True, aggregate_pubmed_ids: bool = True,aggregate_methods: bool = True,nodes_ids= None,
                 node_types: Union[None, list[IRefIndexNodeType]] = None,node_fields:Union[None, list[IRefIndexEdgeFields]] = None,
                 edge_types: Union[None, list[IRefIndexEdgeType]] = None, edge_fields:Union[None, list[IRefIndexEdgeFields]] = None,
                 ):
        
        self.output_dir = output_dir
        self.irefindex_fields = irefindex_fields
        self.add_prefix= add_prefix 
        self.nodes_ids = nodes_ids
        self.add_prefix = add_prefix

        self.aggregate_dict = {IRefIndexEdgeFields.PUBMED_IDS.value:aggregate_pubmed_ids,
                              IRefIndexEdgeFields.METHODS.value:aggregate_methods}
        
        self._set_types_and_fields(node_types, node_fields, edge_types, edge_fields)
        
    def _set_types_and_fields(self, node_types, node_fields, edge_types, edge_fields):

        if node_types:
            self.node_types = node_types
        else:
            self.node_types = [type for type in IRefIndexNodeType]

        if node_fields:
            self.node_fields = node_fields
        else:
            self.node_fields = [field for field in chain()]

        if edge_types:
            self.edge_types = edge_types
        else:
            self.edge_types = [type for type in IRefIndexEdgeType]

        if edge_fields:
            self.edge_fields = edge_fields
        else:
            self.edge_fields = [field for field in chain(IRefIndexEdgeFields)]
      
    def irefindex_process(self):

        logger.info("EXECUTING THE ADAPTER FILE:")
        logger.info("Extracting information from IRefIndex data")
        
        IRefIndexInteraction = collections.namedtuple(
            'IRefIndexInteraction',
            (
                'partner_a',
                'partner_b',
                'pmid', 
                'method',
                'taxon', 
                'relationship_id',
            ),
        )
        
        logger.info("1) Created a collection 'IRefIndexInteractions' for the variables: ")
        logger.info("2) Data is downloaded by the resource as a zip file and will be unzipped in a certain folder --> location of the text file that is needed")
        
        input_file_dir = '.cache/IRefIndex/7227.mitab.08-28-2023.txt.zip.unzip'
        logger.info("--> Location of the text file: {}" .format(input_file_dir))

        # Iterate over the files in the extracted directory and parse each file as needed
        for filename in os.listdir(input_file_dir):
            if filename.endswith(".txt"):
                inputfile_path = os.path.join(input_file_dir, filename)
        logger.info("3) Open/read the text file ")
        logger.info(" --> inputfile:{}".format(inputfile_path))

        # parsing text file 
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
                
                # PARTNER_B: FinalReference B 
                input_partner_b = l[39]

                # skip lines that start with complex 
                if input_partner_a.startswith('complex:') or input_partner_b.startswith('complex:'):
                    continue
                if input_partner_a.startswith('refseq:') or input_partner_b.startswith('refseq:'):
                    continue
                if input_partner_a.startswith('pdb:') or input_partner_b.startswith('pdb:'):
                    continue
                if input_partner_a.startswith('dbj/embl/genbank:') or input_partner_b.startswith('dbj/embl/genbank:'):
                    continue
                if input_partner_a.startswith('flybase:') or input_partner_b.startswith('flybase:'):
                    continue
                
                #isolate id for partner_a
                parts = input_partner_a.split(":")
                partner_a = parts[1] if len(parts) > 1 else ""

                input_partner_a = l[38]
                parts = input_partner_a.split(":") # Splitting the string at ":"
                if len(parts) > 1:
                    partner_a = parts[1] # Selecting the second part after ":"
                else:
                    partner_a = ""

                #isolate id for partner_b
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

                # relationship id 
                input_relationhsip_id= l[13]
                pattern_relationship_id = r"rigid:([^|]+)"

                # Search for the pattern in the column value
                match_relationship_id = re.search(pattern_relationship_id, input_relationhsip_id)

                # Extract the matched value if found
                if match_relationship_id:
                    relationship_id = match_relationship_id.group(1)
                    
                else:
                    relationship_id =""

                interactions.append(
                    IRefIndexInteraction(
                        partner_a=partner_a,
                        partner_b=partner_b,
                        pmid=pmid,
                        method=method,
                        taxon=taxon,
                        relationship_id=relationship_id,
                    )
                )
            #logger.info(partner_a)  
            #logger.info(partner_b)
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

        logger.info("6) Started processing IRefIndex data")
            
        default_field_names = {"pmid":"pubmed_ids", "taxon":"taxon" ,"method":"method", "relationship_id":"relationship_id"}
        
        self.irefindex_field_new_names = {}
        
        if rename_selected_fields:
            if len(selected_fields) != len(rename_selected_fields):
                raise Exception("Length of selected_fields variable should be equal to length of rename_selected_fields variable")
            
            for field_old_name, field_new_name in list(zip(selected_fields, rename_selected_fields)):
                self.irefindex_field_new_names[field_old_name] = field_new_name
                

            self.irefindex_field_new_names["partner_a"] = "uniprot_a"
            self.irefindex_field_new_names["partner_b"] = "uniprot_b"

        else:
            for field_old_name in selected_fields:
                self.irefindex_field_new_names[field_old_name] = default_field_names[field_old_name]
            
            self.irefindex_field_new_names["partner_a"] = "uniprot_a"
            self.irefindex_field_new_names["partner_b"] = "uniprot_b"
        
        logger.info("--> Changing 'partner_a' and 'partner_b' to 'uniprot_a' and 'uniprot_b'")
        
        t1 = time()
                        
        # create dataframe     
        irefindex_data_without_headers = self.irefindex_ints[1:]     
        irefindex_df = pd.DataFrame.from_records(irefindex_data_without_headers, columns=self.irefindex_ints[0]._fields)
        #logger.info("initial dataframe:{}".format(irefindex_df))

        #logger.info("create irefindex_dataframe:{}".format(irefindex_df)) # --> partner_a |partner_b |pmid |method |taxon
        logger.info("--> Created an irefindex_dataframe")

        
        # add source database info
        irefindex_df["source"] = "IRefIndex"

        # filter selected fields
        irefindex_df = irefindex_df[list(self.irefindex_field_new_names.keys())]
        #logger.info("filter irefindex_dataframe:{}".format(irefindex_df)) # --> pmid| taxon| method| partner_a |partner_b
        logger.info("--> Filtered the irefindex_dataframe")
        
        # rename columns
        irefindex_df.rename(columns=self.irefindex_field_new_names, inplace=True)
        #logger.info("rename irefindex_dataframe:{}".format(irefindex_df))
        logger.info("--> Renamed the headers of the irefindex_dataframe to uniprot_a and uniprot_b")
        
        #drop complexes and refseq 
        

        # drop rows if uniprot_a or uniprot_b is not a swiss-prot protein
        # this is to strict --> removes id like A0A021WW32 & X2JKM9(trembl) --> not reviewed by humans
        #with swissprot= none you do look at reviewed and unreviewed 
        #irefindex_df = irefindex_df[(irefindex_df["uniprot_a"].isin(self.swissprots)) & (irefindex_df["uniprot_b"].isin(self.swissprots))]
        #irefindex_df.reset_index(drop=True, inplace=True)
        #logger.info("drop in irefindex_dataframe if not in swissprot :{}".format(irefindex_df)) # --> pubmed_ids | taxon | method |uniprot_a |uniprot_b
        #logger.info(" --> Dropped the rows in irefindex_dataframe if uniprot_a or uniprot_b is not in swissprot")
        
        # drop duplicates if same a x b pair exists multiple times 
        # keep the first pair and collect pubmed ids of duplicated a x b pairs in that pair's pubmed id column
        # if a x b pair has same experimental system type with b x a pair, drop b x a pair
        irefindex_df_unique = irefindex_df.dropna(subset=["uniprot_a", "uniprot_b"]).reset_index(drop=True)        
        #logger.info(" --> Created a unique irefindex_dataframe, removed duplicates:{} ".format(irefindex_df_unique))
        
        
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
        
        # group by unique combinations of uniprot a and b
        irefindex_df_unique = irefindex_df_unique.groupby(["uniprot_a", "uniprot_b"], sort=False, as_index=False).aggregate(agg_dict)
        #biogrid_df_unique["pubmed_id"].replace("", np.nan, inplace=True)
        logger.info("9) Group the irefindex_dataframe based in uniprot_a and uniprot_b")
        
        if "method" in self.irefindex_field_new_names.keys():            
            irefindex_df_unique = irefindex_df_unique[~irefindex_df_unique[["uniprot_a", "uniprot_b", self.irefindex_field_new_names["method"]]].apply(frozenset, axis=1).duplicated()].reset_index(drop=True)
        else:
            irefindex_df_unique = irefindex_df_unique[~irefindex_df_unique[["uniprot_a", "uniprot_b"]].apply(frozenset, axis=1).duplicated()].reset_index(drop=True)

        self.final_irefindex_ints = irefindex_df_unique
        logger.info("FINAL DATAFRAME:\n{}".format(self.final_irefindex_ints)) 
          
        # Extract unique UniProt IDs from both columns
        nodes_a = irefindex_df_unique["uniprot_a"].unique()
        nodes_b = irefindex_df_unique["uniprot_b"].unique()
        self.nodes_ids = set(nodes_a) | set(nodes_b)
        
        node_id_to_taxon = dict(zip(irefindex_df_unique["uniprot_a"], irefindex_df_unique["taxon"]))
        node_id_to_taxon.update(zip(irefindex_df_unique["uniprot_b"], irefindex_df_unique["taxon"]))

        node_id_to_pubmed_id = dict(zip(irefindex_df_unique["uniprot_a"], irefindex_df_unique["pubmed_ids"]))
        node_id_to_pubmed_id.update(zip(irefindex_df_unique["uniprot_b"], irefindex_df_unique["pubmed_ids"]))

        node_id_to_method = dict(zip(irefindex_df_unique["uniprot_a"], irefindex_df_unique["method"]))
        node_id_to_method.update(zip(irefindex_df_unique["uniprot_b"], irefindex_df_unique["method"]))
        self.nodes = []
         
        t3= time()
        if IRefIndexNodeType.PROTEIN in self.node_types:
            logger.info("Getting pubmed ids")
            for node_id in self.nodes_ids:
        # Create Protein object for the current node_id with corresponding taxon
                taxon = node_id_to_taxon.get(node_id, None)        
                pubmed_id= node_id_to_pubmed_id.get(node_id, None)
                method= node_id_to_method.get(node_id, None)
                self.nodes.append(Protein(node_id=node_id, taxon=taxon,pubmed_id= pubmed_id,method =method,  fields=self.node_fields))
                #self.nodes.append(Protein(node_id=node_id, taxon=taxon, pubmed_id= pubmed_id))

        for node in self.nodes:
            yield (node.get_id(), node.get_label(), node.get_properties())

        t4= time()
        logger.info(f'Getting information like pubmed_id and taxon from dataframe in {round((t4-t3) / 60, 2)} mins')
        t2 = time()
        logger.info(f'IRefIndex data is processed in {round((t2-t1) / 60, 2)} mins')
         
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
    
        
        return identifier # --> none 
       


        #https://github.com/biocypher/dependency-map/blob/main/dmb/adapter.py#L342 --> line 312 ( get nodes)
        
        

    def get_edges(self, label: str = "protein_protein_interaction") -> list[tuple]:
        """
        Get  edges from merged data
        Args:
            label: label of protein-protein interaction edges
        """

        # create edge list
        edge_list = []
        for index, row in tqdm(self.final_irefindex_ints.iterrows(), total=self.final_irefindex_ints.shape[0]):
            _dict = row.to_dict()

            _source = self.add_prefix_to_id(identifier=str(row["uniprot_a"]))
            _target = self.add_prefix_to_id(identifier=str(row["uniprot_b"]))

            del _dict["uniprot_a"], _dict["uniprot_b"]

            _props = {}
            for k, v in _dict.items():
                if str(v) != "nan":
                    if isinstance(v, str) and "|" in v:
                        _props[str(k).replace(" ", "_").lower()] = v.replace(
                            "'", "^"
                        ).split("|")
                    else:
                        _props[str(k).replace(" ", "_").lower()] = str(
                            v
                        ).replace("'", "^")

            edge_list.append((None, _source, _target, label, _props))

        #logger.info(edge_list)
        return edge_list
    
    # relationship id in irefindex database --> interaction id
        



class Node:
    """
    Base class for nodes.
    """

    def __init__(self):
        self.id = None
        self.label = None
        self.properties = {}

    def get_id(self):
        """
        Returns the node id.
        """
        return self.id

    def get_label(self):
        """
        Returns the node label.
        """
        return self.label

    def get_properties(self):
        """
        Returns the node properties.
        """
        return self.properties


class Protein(Node):
    """
    Generates instances of proteins.
    """

    def __init__(self, node_id:str, taxon:str, pubmed_id: str, method: str, fields: Optional[list] = None):
        self.fields = fields
        self.id = node_id
        self.label = "uniprot_protein"
        self.properties = self._generate_properties(taxon, pubmed_id,method)
        self.taxon = taxon
        self.pubmed_id =pubmed_id
        self.method = method 
        
    def _generate_properties(self,taxon,pubmed_id,method):
        properties = {}

        ## pmid
        if self.fields is not None and IRefIndexEdgeFields.PUBMED_IDS in self.fields:
           properties["pmid"] = pubmed_id
       
        ## taxon
        if self.fields is not None and IRefIndexEdgeFields.TAXON in self.fields:
           properties["taxon"] = taxon

        ## method
        if self.fields is not None and IRefIndexEdgeFields.TAXON in self.fields:
           properties["method"] = method

        return properties