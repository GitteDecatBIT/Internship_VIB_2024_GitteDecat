#!/usr/bin/env python
import random
import string
from enum import Enum, auto
from itertools import chain
from typing import Optional
from biocypher._logger import logger


# extra from CROssBAR
from contextlib import ExitStack
from time import time
from template_package.adapters.CROssBAR_IRefIndex_pypath_url import url as irefindex_url
from pypath.share import curl, settings
import re
from typing import Union
#from create_knowledge_graph_IrefIndex import paths
import os 
import collections
import pandas as pd
import numpy as np
from pypath.inputs import uniprot
#############################################################################################################


logger.debug(f"Loading module {__name__}.")

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
                 export_csv = False, 
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
                 node_types: Union[None, list[IRefIndexEdgeFields]] = None,
                 node_fields:Union[None, list[IRefIndexEdgeFields]] = None,
                 edge_types: Union[None, list[IRefIndexEdgeFields]] = None,
                 edge_fields:Union[None, list[IRefIndexEdgeFields]] = None,
                 rename_selected_fields: Union[None, list[str]] = None,
            
                 ):
        self.export_csv = export_csv
        self.swissprots = list(uniprot._all_uniprots("*", True))
        self.irefindex_fields = irefindex_fields

        self.aggregate_dict = {IRefIndexEdgeFields.PUBMED_IDS.value:aggregate_pubmed_ids,
                              IRefIndexEdgeFields.METHODS.value:aggregate_methods}
        

        self._set_types_and_fields(node_types, node_fields, edge_types, edge_fields)


        
    def get_irefindex_data(self):

        # Extracting information from IRefIndex
        logger.info("Extracting information from IRefIndex")
        
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

        logger.info("Created a collection of the variables: IRefIndexInteractions")

        # Path to the extracted data directory
        input_file_dir = '.cache/IRefIndex/7227.mitab.08-28-2023.txt.zip.unzip'
        logger.info("input file dir: {}" .format(input_file_dir))
        # Iterate over the files in the extracted directory and parse each file as needed
        for filename in os.listdir(input_file_dir):
            if filename.endswith(".txt"):
                inputfile_path = os.path.join(input_file_dir, filename)
        
        logger.info("inputfile:{}".format(inputfile_path))

 
        logger.info("Exctracting columns from the IRefIndex database")

        interactions = []

        # Open the file
        with open(inputfile_path, 'r') as file:
            logger.info("Getting inforamtion: partner_a, partner_b, pubmed_id, method and taxon id")
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
                pattern_pmid = r'\d+'
                pmid = re.findall(pattern_pmid, input_pmid)

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
                
            logger.info("Succesfully extracted columns from the IRefIndex database!")
            self.irefindex_ints = interactions
            return interactions
        
    

    def irefindex_process(self, rename_selected_fields: Union[None, list[str]] = None) -> None:
        """
        Processor function for irefindex data. It drops duplicate and reciprocal duplicate protein pairs and collects pubmed ids of duplicated pairs. In addition, it
        maps entries to uniprot ids using gene name and tax id information in the BioGRID data. Also, it filters protein pairs found in swissprot.
        
         Args:
            rename_selected_fields : List of new field names for selected fields. If not defined, default field names will be used.
        """
        
        selected_fields = self.set_edge_fields()
            
        default_field_names = { "pmid":"pubmed_ids", "taxon":"taxon" ,"method":"method"}
        
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
                logger.info("field_old_name:{}".format(field_old_name))
                self.irefindex_field_new_names[field_old_name] = default_field_names[field_old_name]
            
            self.irefindex_field_new_names["partner_a"] = "uniprot_a"
            self.irefindex_field_new_names["partner_b"] = "uniprot_b"
        
        
        logger.info("Started processing BioGRID data")
        t1 = time()
                         
        # create dataframe          
        irefindex_df = pd.DataFrame.from_records(self.irefindex_ints, columns=self.irefindex_ints[0]._fields)

        logger.info("create irefindex_dataframe:{}".format(irefindex_df)) # --> partner_a     partner_b        pmid method taxon


        
        # add source database info
        irefindex_df["source"] = "IRefIndex"
        # filter selected fields
        irefindex_df = irefindex_df[list(self.irefindex_field_new_names.keys())]
        logger.info("filter irefindex_dataframe:{}".format(irefindex_df)) # --> pmid taxon                          method                    partner_a     partner_b
        # rename columns
        irefindex_df.rename(columns=self.irefindex_field_new_names, inplace=True)
        logger.info("rename irefindex_dataframe:{}".format(irefindex_df))
        
        
        # drop rows if uniprot_a or uniprot_b is not a swiss-prot protein
        irefindex_df = irefindex_df[(irefindex_df["uniprot_a"].isin(self.swissprots)) & (irefindex_df["uniprot_b"].isin(self.swissprots))]
        irefindex_df.reset_index(drop=True, inplace=True)
        logger.info("drop in irefindex_dataframe if not in swissprot :{}".format(irefindex_df)) # --> pubmed_ids taxon                          method                    uniprot_a     uniprot_b

        
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
        
        irefindex_df_unique = irefindex_df_unique.groupby(["uniprot_a", "uniprot_b"], sort=False, as_index=False).aggregate(agg_dict)
        #biogrid_df_unique["pubmed_id"].replace("", np.nan, inplace=True)
        
        if "method" in self.irefindex_field_new_names.keys():            
            irefindex_df_unique = irefindex_df_unique[~irefindex_df_unique[["uniprot_a", "uniprot_b", self.irefindex_field_new_names["method"]]].apply(frozenset, axis=1).duplicated()].reset_index(drop=True)
        else:
            irefindex_df_unique = irefindex_df_unique[~irefindex_df_unique[["uniprot_a", "uniprot_b"]].apply(frozenset, axis=1).duplicated()].reset_index(drop=True)

        if self.export_csvs:
            irefindex_output_path = self.export_dataframe(irefindex_df_unique, "biogrid")
            logger.info(f'Final BioGRID data is written: {irefindex_output_path}')

        self.final_irefindex_ints = irefindex_df_unique
        
        t2 = time()
        logger.info(f'BioGRID data is processed in {round((t2-t1) / 60, 2)} mins')
            
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

    # maybe put the parsing in a differen def????
    def get_nodes(self):
        """
        Returns a generator of node tuples for node types specified in the
        adapter constructor.
        """
        logger.info("Generating nodes.")
        
        self.nodes = []

        if IRefIndexEdgeFields.PROTEIN in self.node_types:
            [self.nodes.append(Protein(fields=self.node_fields))for _ in range(100)]


        for node in self.nodes:
            yield (node.get_id(), node.get_label(), node.get_properties())


        #https://github.com/biocypher/dependency-map/blob/main/dmb/adapter.py#L342 --> line 312 ( get nodes)
        
        
    



    # duplicates are found!!!! see CRosBAR file for this 
            
    def get_edges(self, probability: float = 0.5):
        """
        Returns a generator of edge tuples for edge types specified in the
        adapter constructor.

        Args:
            probability: Probability of generating an edge between two nodes.
        """

        logger.info("Generating edges.")

        if not self.nodes:
            raise ValueError("No nodes found. Please run get_nodes() first.")

        for node in self.nodes:
            if random.random() < probability:
                other_node = random.choice(self.nodes)

                # generate random relationship id by choosing upper or lower
                # letters and integers, length 10, and joining them
                relationship_id = "".join(
                    random.choice(string.ascii_letters + string.digits)
                    for _ in range(10)
                )

                # determine type of edge from other_node type
                if (
                    isinstance(other_node, Protein)
                    and IRefIndexEdgeFields.PROTEIN_PROTEIN_INTERACTION
                    in self.edge_types
                ):
                    edge_type = IRefIndexEdgeFields.PROTEIN_PROTEIN_INTERACTION.value
                else:
                    continue

                yield (
                    relationship_id,
                    node.get_id(),
                    other_node.get_id(),
                    edge_type,
                    {"example_property": "example_value"}, #### ????
                )


        """
        Returns the number of nodes generated by the adapter.
        """
        return len(list(self.get_nodes()))

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
## specify what is uploaded to the exel 
    #:ID, sequence, description, taxon, id, preferred id (comes from schema_config.yaml), :LABEL


class Protein(Node):
    """
    Generates instances of proteins.
    """

    def __init__(self, fields: Union[None, list[IRefIndexEdgeFields]] = None):
        self.fields = fields
        self.id = self._get_id()
        self.label = "uniprot_protein"
        self.properties = self._get_properties()


    # ?????
    # id is a function and join function does not work in biocypher _batch_writer.py on line 643
        # lines.append(self.delim.join(line) + "\n")
        #line = [n.get_id()]
        #_id = node.get_id()

    def _get_id(self):
        """
        Get uniprot id.
        """
        id_partner_a= interactions.partner_a
        

        return id_partner_a


    def _get_properties(self):
        properties = {}

        ## pmid
        if self.fields is not None and IRefIndexAdapterProteinNodeField.TAXON in self.fields:
            properties["pmid"] = irefindex_pmids()
        
        ## method
        #if self.fields is not None and IRefIndexAdapterProteinField.TAXON in self.fields:
        #    properties["method"] = irefindex_method
        
        ## organism(taxon)
        if self.fields is not None and IRefIndexAdapterProteinNodeField.TAXON in self.fields:
            properties["taxon"] = irefindex_species()

        ## source ??
            
        return properties


# output in excel: 
    # :ID /pmid/taxon/id/preffered_id/ :LABEL
    # preffered id is specified in the chema_config.yaml
    # label comes from???? 

