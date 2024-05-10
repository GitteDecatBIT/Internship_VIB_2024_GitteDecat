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
#############################################################################################################


logger.debug(f"Loading module {__name__}.")



class IRefIndexAdapterNodeType(Enum):
    """
    Define types of nodes the adapter can provide.
    """

    PROTEIN = auto()
    


class IRefIndexAdapterProteinNodeField(Enum):
    """
    Define possible fields the adapter can provide for proteins.
    """

    PUBMED_ID = "pmid"
    TAXON = "taxon"
    #SOURCE = "source"
    

class IRefIndexAdapterEdgeType(Enum):
    """
    Enum for the types of the protein adapter.
    """

    PROTEIN_PROTEIN_INTERACTION = "protein_protein_interaction"
    


class IRefIndexAdapterProteinProteinEdgeField(Enum):
    """
    Define possible fields the adapter can provide for protein-protein edges.
    """

    INTERACTION_TYPE = "interaction_type" # column [11]
    INTERACTION_SOURCE = "interaction_source" # column [12]
    METHOD = "method"

    

    

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
    
    def __init__(
        self,
        node_types: Union[None, list[IRefIndexAdapterNodeType]] = None,
        node_fields:Union[None, list[IRefIndexAdapterProteinNodeField]] = None,
        edge_types: Union[None, list[IRefIndexAdapterEdgeType]] = None,
        edge_fields:Union[None, list[IRefIndexAdapterProteinProteinEdgeField]] = None,
        input_dir= None, 
        #organism= organism, 
        
    ):
        self._set_types_and_fields(node_types, node_fields, edge_types, edge_fields)

    # maybe put the parsing in a differen def????
    def get_nodes(self):
        """
        Returns a generator of node tuples for node types specified in the
        adapter constructor.
        """
        logger.info("Generating nodes.")
        
        self.nodes = []

        if IRefIndexAdapterNodeType.PROTEIN in self.node_types:
            [self.nodes.append(Protein(fields=self.node_fields))for _ in range(100)]

        for node in self.nodes:
            yield (node.get_id(), node.get_label(), node.get_properties())
        
        #https://github.com/biocypher/dependency-map/blob/main/dmb/adapter.py#L342 --> line 312 ( get nodes)
            
        # Extracting information from IRefIndex
        logger.info("Extracting information from IRefIndex")
        
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

        logger.info("Created a collection of the variables: IRefIndexPhysicalInteractions")

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
            logger.info("Getting inforamtion: partner_a, partner_b, pubmed_id, method and organism")
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

                # ORGANISM
                input_organism= l[10]
                pattern_organism= r'taxid:(\d+)'
                match_organism= re.search(pattern_organism, input_organism)
                if match_organism:
                    organism = match_organism.group(1)
                else:
                    organism = ""

                interactions.append(
                    IRefIndexPhysicalInteraction(
                        partner_a=partner_a,
                        partner_b=partner_b,
                        pmid=pmid,
                        method=method,
                        organism=organism,
                    )
                )
            #logger.info(interactions)
        
            logger.info("Succesfully extracted columns from the IRefIndex database!")
            return interactions
           
        




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
                    and IRefIndexAdapterEdgeType.PROTEIN_PROTEIN_INTERACTION
                    in self.edge_types
                ):
                    edge_type = IRefIndexAdapterEdgeType.PROTEIN_PROTEIN_INTERACTION.value
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
            self.node_types = [type for type in IRefIndexAdapterNodeType]

        if node_fields:
            self.node_fields = node_fields
        else:
            self.node_fields = [
                field
                for field in chain(
                    IRefIndexAdapterProteinNodeField,
        
                )
            ]

        if edge_types:
            self.edge_types = edge_types
        else:
            self.edge_types = [type for type in IRefIndexAdapterEdgeType]

        if edge_fields:
            self.edge_fields = edge_fields
        else:
            self.edge_fields = [
                field 
                for field in chain(
                IRefIndexAdapterProteinProteinEdgeField,
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

    def __init__(self, fields: Union[None, list[IRefIndexAdapterProteinNodeField]] = None):
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
        id_partner_a= partner_a
        

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
