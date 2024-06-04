#!/usr/bin/env python
from biocypher._logger import logger
from enum import Enum, auto
from typing import Union
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
    PROTEIN_IREFINDEX= auto()

class IRefIndexNodeFields(Enum):
    PUBMED_IDS = "pmid"
    TAXON= "taxon"
    METHODS = "method"

class IRefIndexEdgeType(Enum):
    """
    Enum for the types of the protein adapter.
    """
    PROTEIN_PROTEIN_INTERACTION = "protein_protein_interaction"

class IRefIndexEdgeFields(Enum):
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
    def __init__(self,
                 output_dir = None,
                 irefindex_fields: Union[None, list[IRefIndexNodeFields, IRefIndexEdgeFields]] = None, 
                 add_prefix = True, 
                 aggregate_pubmed_ids: bool = True,
                 aggregate_methods: bool = True,
                 nodes_ids= None,
                 node_types: Union[None, list[IRefIndexNodeType]] = None,
                 node_fields:Union[None, list[IRefIndexNodeFields]] = None,
                 edge_types: Union[None, list[IRefIndexEdgeType]] = None, 
                 edge_fields:Union[None, list[IRefIndexEdgeFields]] = None,
                 ):
        
        self.output_dir = output_dir
        self.irefindex_fields = irefindex_fields
        self.add_prefix= add_prefix 
        self.nodes_ids = nodes_ids
        self.add_prefix = add_prefix
        self.aggregate_dict = {IRefIndexNodeFields.PUBMED_IDS.value:aggregate_pubmed_ids,
                              IRefIndexNodeFields.METHODS.value:aggregate_methods}
        
        self._set_types_and_fields(node_types, node_fields, edge_types, edge_fields)
        
    def _set_types_and_fields(self, node_types,node_fields, edge_types, edge_fields):

        if node_types:
            self.node_types = node_types
        else:
            self.node_types = [type for type in IRefIndexNodeType]
        if node_fields:
            self.node_fields = node_fields
        else:
            self.node_fields = [field for field in chain(IRefIndexNodeFields)]
        if edge_types:
            self.edge_types = edge_types
        else:
            self.edge_types = [type for type in IRefIndexEdgeType]

        if edge_fields:
            self.edge_fields = edge_fields
        else:
            self.edge_fields = [field for field in chain(IRefIndexEdgeFields)]
      
    def irefindex_process(self):

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
                
        # Get the current directory
        current_directory = os.getcwd()

        # Iterate through the files in the current directory
        for root, dirs, files in os.walk(current_directory):
            for file in files:
                if file.endswith(".txt"):
                    inputfile= os.path.join(root, file)
        logger.info("This is the input file that is used: {}".format(inputfile))
            

        # parsing text file 

        interactions = []
        
        # Open the file
        with open(inputfile, 'r') as file:
            logger.info("Getting inforamtion for partner_a, partner_b, pubmed_id, method, taxon id and relationship id")
            for line in file:
                # Split the line by tab character
                line = line.split('\t')
                
                # PARTNER_A: finalReference A 
                input_partner_a = line[38]
                
                # PARTNER_B: FinalReference B 
                input_partner_b = line[39]

                # skip lines that start with complex 
                if input_partner_a.startswith('complex:') or input_partner_b.startswith('complex:'):
                    continue
               
                if input_partner_a.startswith('pdb:') or input_partner_b.startswith('pdb:'):
                    continue

                if input_partner_a.startswith('flybase:') or input_partner_b.startswith('flybase:'):
                    continue
                
                #isolate id for partner_a
                parts = input_partner_a.split(":")
                partner_a = parts[1] if len(parts) > 1 else ""

                input_partner_a = line[38]
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
                input_pmid = line[8]
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
                input_method= line[6]
                pattern_method= r'\((.*?)\)'
                match_method= re.search(pattern_method, input_method)
                if match_method:
                    method = match_method.group(1) # Extracting the text between brackets
                else:
                    method = ""

                # ORGANISM/taxon 
                input_taxon= line[10]
                pattern_taxon= r'taxid:(\d+)'
                match_taxon= re.search(pattern_taxon, input_taxon)
                if match_taxon:
                    taxon = match_taxon.group(1)
                else:
                    taxon = ""

                # relationship id 
                input_relationhsip_id= line[13]
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
                # yield

            logger.info("--> Succesfully extracted information from the IRefIndex database!")
            self.irefindex_ints = interactions
            return interactions
        


    def get_nodes(self, rename_selected_fields: Union[None, list[str]] = None):
        """
        Processor function for irefindex data. It drops duplicate and reciprocal duplicate protein pairs and collects pubmed ids of duplicated pairs. In addition, it
        maps entries to uniprot ids using gene name and tax id information in the BioGRID data. Also, it filters protein pairs found in swissprot.
        
         Args:
            rename_selected_fields : List of new field names for selected fields. If not defined, default field names will be used.
        """
        logger.info("Generating nodes.")

        selected_node_fields  = self.set_node_fields()

        selected_edge_fields  = self.set_edge_fields()
        selected_fields= selected_node_fields + selected_edge_fields

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
        
        t1 = time()
                        
        # create dataframe  
        logger.info(" Creating a dataframe from the input")   
        irefindex_data_without_headers = self.irefindex_ints[1:]     
        irefindex_df = pd.DataFrame.from_records(irefindex_data_without_headers, columns=self.irefindex_ints[0]._fields)

        # add source database info
        irefindex_df["source"] = "IRefIndex"

        # filter selected fields
        irefindex_df = irefindex_df[list(self.irefindex_field_new_names.keys())]
        
        # rename columns
        irefindex_df.rename(columns=self.irefindex_field_new_names, inplace=True)        

        # drop duplicates if same a x b pair exists multiple times 
        # keep the first pair and collect pubmed ids of duplicated a x b pairs in that pair's pubmed id column
        # if a x b pair has same experimental system type with b x a pair, drop b x a pair
        irefindex_df_unique = irefindex_df.dropna(subset=["uniprot_a", "uniprot_b"]).reset_index(drop=True)  
                
        def aggregate_fields(element):
            element = "|".join([str(e) for e in set(element.dropna())])
            if not element:
                return np.nan
            else:
                return element

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

        # Drop any remaining duplicates just in case
        irefindex_df_unique = irefindex_df_unique.drop_duplicates(subset=["uniprot_a", "uniprot_b"])

        from create_knowledge_graph import get_taxon_id
        taxon_id = get_taxon_id()

        # Check if the taxon_id is present in the DataFrame
        if taxon_id not in irefindex_df_unique['taxon'].values:
            error_message = f"Taxon ID {taxon_id} not found in the DataFrame."
            logger.error(error_message)
            raise ValueError(error_message)
        
        irefindex_df_unique = irefindex_df_unique[irefindex_df_unique['taxon'] == taxon_id]

        
        if "method" in self.irefindex_field_new_names.keys():            
            irefindex_df_unique = irefindex_df_unique[~irefindex_df_unique[["uniprot_a", "uniprot_b", self.irefindex_field_new_names["method"]]].apply(frozenset, axis=1).duplicated()].reset_index(drop=True)
        else:
            irefindex_df_unique = irefindex_df_unique[~irefindex_df_unique[["uniprot_a", "uniprot_b"]].apply(frozenset, axis=1).duplicated()].reset_index(drop=True)

        self.final_irefindex_ints = irefindex_df_unique
        logger.info("FINAL DATAFRAME:\n{}".format(self.final_irefindex_ints)) 
                

        # Function to determine if an ID is a UniProt ID
        def is_uniprot_id(node_id):
            return bool(re.match(r'^([A-N,R-Z][0-9]([A-Z][A-Z, 0-9][A-Z, 0-9][0-9]){1,2})|([O,P,Q][0-9][A-Z, 0-9][A-Z, 0-9][A-Z, 0-9][0-9])(\.\d+)?$', node_id))
        
        # Function to determine if an ID is a refseq ID
        def is_refseq_id(node_id):
            return bool(re.match(r'^(((AC|AP|NC|NG|NM|NP|NR|NT|NW|WP|XM|XP|XR|YP|ZP)_\d+)|(NZ_[A-Z]{2,4}\d+))(\.\d+)?', node_id))
        
        def is_entrez_id(node_id):
            return bool(re.match(r'^[A-Z]+[0-9]+(\.\d+)?$', node_id))


        # Extract unique IDs from both columns
        nodes_a = irefindex_df_unique["uniprot_a"].unique()
        nodes_b = irefindex_df_unique["uniprot_b"].unique()
        nodes_ids = set(nodes_a) | set(nodes_b)

        # Map node IDs to their properties
        node_id_to_taxon = dict(zip(irefindex_df_unique["uniprot_a"], irefindex_df_unique["taxon"]))
        node_id_to_taxon.update(zip(irefindex_df_unique["uniprot_b"], irefindex_df_unique["taxon"]))

        node_id_to_pubmed_id = dict(zip(irefindex_df_unique["uniprot_a"], irefindex_df_unique["pubmed_ids"]))
        node_id_to_pubmed_id.update(zip(irefindex_df_unique["uniprot_b"], irefindex_df_unique["pubmed_ids"]))

        node_id_to_method = dict(zip(irefindex_df_unique["uniprot_a"], irefindex_df_unique["method"]))
        node_id_to_method.update(zip(irefindex_df_unique["uniprot_b"], irefindex_df_unique["method"]))

        self.nodes = []

        # Example field list, replace with actual field list from your context
        self.node_fields = ["pubmed_ids", "taxon", "method"]

        #from create_knowledge_graph_IREFINDEX_06 import get_taxon_id
        #taxon_id = get_taxon_id()
#
        #filtered_nodes = [node for node in nodes_ids if node_id_to_taxon.get(node) == taxon_id]

        if IRefIndexNodeType.PROTEIN_IREFINDEX in self.node_types:
            for node_id in nodes_ids:
                taxon = node_id_to_taxon.get(node_id, None)
                pubmed_id = node_id_to_pubmed_id.get(node_id, None)
                method = node_id_to_method.get(node_id, None)

                if is_uniprot_id(node_id):
                    protein_type = "uniprot_protein"

                elif is_refseq_id(node_id):
                    protein_type = "refseq_protein"

                elif is_entrez_id(node_id):
                    protein_type = "entrez_protein"

                self.nodes.append(Protein(node_id=node_id, protein_type=protein_type, taxon=taxon, pubmed_id=pubmed_id, method=method, fields=self.node_fields))

        for node in self.nodes:
            yield (node.get_id(), node.get_label(), node.get_properties())

        t2 = time()
        logger.info(f'IRefIndex data is processed in {round((t2-t1) / 60, 2)} mins')
         
    def set_node_fields(self) -> list:
        """
        Returns:
            selected field list
        """
        if self.irefindex_fields is None:
            return [field.value for field in IRefIndexNodeFields]
        else:
            return [field.value for field in self.irefindex_fields]
 
    def set_edge_fields(self) -> list:
        """
        Returns:
            selected field list
        """
        if self.irefindex_fields is None:
            return [field.value for field in IRefIndexEdgeFields]
        else:
            return [field.value for field in self.irefindex_fields]
        
    def add_prefix_to_id(self, identifier=None, sep=":") -> str:
        """
        Adds prefix to the protein id based on the identifier type
        """ 
        # Function to determine if an ID is a UniProt ID
        def is_uniprot_id(identifier):
            return bool(re.match(r'^([A-N,R-Z][0-9]([A-Z][A-Z, 0-9][A-Z, 0-9][0-9]){1,2})|([O,P,Q][0-9][A-Z, 0-9][A-Z, 0-9][A-Z, 0-9][0-9])(\.\d+)?(\-\d+)?$', identifier))
        
        # Function to determine if an ID is a refseq ID
        def is_refseq_id(identifier):
            return bool(re.match(r'^(((AC|AP|NC|NG|NM|NP|NR|NT|NW|WP|XM|XP|XR|YP|ZP)_\d+)|(NZ_[A-Z]{2,4}\d+))(\.\d+)?$', identifier))
        
        def is_entrez_id(identifier):
            return bool(re.match(r'^[A-Z]+[0-9]+(\.\d+)?$', identifier))
        
        if self.add_prefix and identifier:
            if is_uniprot_id(identifier):
                prefix = "uniprot"
            elif is_refseq_id(identifier):
                prefix = "refseq"
            elif is_entrez_id(identifier):
                prefix = "entrez"
            else:
                prefix= ""

            return normalize_curie(prefix + sep + str(identifier))
        
        return identifier
        
        

    def get_edges(self, label: str = "protein_protein_interaction") -> list[tuple]:
        """
        Get  edges from merged data
        Args:
            label: label of protein-protein interaction edges --> Must be same as in the schema_config.yaml file 
        """
        logger.info("Generating edges.")

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

        
        return edge_list
        
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

    def __init__(self, node_id:str, protein_type: str, taxon:str, pubmed_id: str, method: str, fields: Optional[list] = None):
        self.fields = fields
        self.id = node_id
        self.label = protein_type
        self.properties = self._generate_properties(taxon,pubmed_id,method)
        self.taxon = taxon
        self.pubmed_id =pubmed_id
        self.method = method 
        
    def _generate_properties(self, taxon, pubmed_id, method):
        properties = {}

        if self.fields is None or "pubmed_ids" in self.fields:
            properties["pubmed_ids"] = pubmed_id

        if self.fields is None or "taxon" in self.fields:
            properties["taxon"] = taxon

        if self.fields is None or "method" in self.fields:
            properties["method"] = method

        return properties
    
