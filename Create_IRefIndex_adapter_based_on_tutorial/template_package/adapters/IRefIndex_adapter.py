import random
import string
from enum import Enum, auto
from itertools import chain
from typing import Optional
from biocypher._logger import logger


# extra from CROssBAR
from contextlib import ExitStack
from time import time
from template_package.adapters.CROssBAR_IRefIndex_input import irefindex_interactions,irefindex_species
from template_package.adapters.CROssBAR_IRefIndex_pypath_url import url as irefindex_url
from pypath.share import curl, settings
import re
from typing import Union

logger.debug(f"Loading module {__name__}.")


class ExampleAdapterNodeType(Enum):
    """
    Define types of nodes the adapter can provide.
    """

    PROTEIN = auto()
    DISEASE = auto()


class ExampleAdapterProteinField(Enum):
    """
    Define possible fields the adapter can provide for proteins.
    """

    ID = "id"
    SEQUENCE = "sequence"
    DESCRIPTION = "description"
    TAXON = "taxon"


class ExampleAdapterDiseaseField(Enum):
    """
    Define possible fields the adapter can provide for diseases.
    """

    ID = "id"
    NAME = "name"
    DESCRIPTION = "description"


class ExampleAdapterEdgeType(Enum):
    """
    Enum for the types of the protein adapter.
    """

    PROTEIN_PROTEIN_INTERACTION = "protein_protein_interaction"
    PROTEIN_DISEASE_ASSOCIATION = "protein_disease_association"


class ExampleAdapterProteinProteinEdgeField(Enum):
    """
    Define possible fields the adapter can provide for protein-protein edges.
    """

    INTERACTION_TYPE = "interaction_type"
    INTERACTION_SOURCE = "interaction_source"


class ExampleAdapterProteinDiseaseEdgeField(Enum):
    """
    Define possible fields the adapter can provide for protein-disease edges.
    """

    ASSOCIATION_TYPE = "association_type"
    ASSOCIATION_SOURCE = "association_source"

class IRefIndexEdgeFields(Enum):
    SOURCE = "source"
    PUBMED_IDS = "pmid"
    METHOD = "method"

    

class ExampleAdapter:
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
        node_types: ExampleAdapterNodeType.PROTEIN,
        node_fields: ExampleAdapterProteinField,
        edge_types: ExampleAdapterEdgeType.PROTEIN_PROTEIN_INTERACTION,
        edge_fields: Union[None, list[IRefIndexEdgeFields]] = None,
        test_mode = False, # change to false if you want the intire dataset
        organism= irefindex_species,
    ):
        self._set_types_and_fields(node_types, node_fields, edge_types, edge_fields)
        self.test_mode = test_mode
        self.organism = organism


 

    def download_irefindex_data(self):
        
        """
        Wrapper function to download IRefIndex data using pypath; used to access
        settings.
            
        To do: Make arguments of irefindex_all_interactions selectable for user. 
        
        """

        t0 = time()

        with ExitStack() as stack:   
            stack.enter_context(settings.context(retries=self.retries))
            
            if self.debug:                
                stack.enter_context(curl.debug_on())
            if not self.cache:
                stack.enter_context(curl.cache_off())

            if self.organism is None:
                species =irefindex_species()
                self.tax_ids = list(species.keys())
            else:
                self.tax_ids = [self.organism]                      
            # download irefindex data
            self.irefindex_ints = irefindex_interactions()


            logger.info(f"This is the link of IRefIndex data we downloaded:{irefindex_url}. Please check if it is up to date")    
            logger.debug("Started downloading IRefIndex data")

        # devides into a smaller subset of data 
        if self.test_mode:
            self.irefindex_ints = self.irefindex_ints[:100]            
                    

        t1 = time()
        logger.info(f'IRefIndex data is downloaded in {round((t1-t0) / 60, 2)} mins')

    download_irefindex_data(self)
            
    def get_nodes(self):
        """
        Returns a generator of node tuples for node types specified in the
        adapter constructor.
        """

        logger.info("Generating nodes.")

        self.nodes = []

        if ExampleAdapterNodeType.PROTEIN in self.node_types:
            [self.nodes.append(Protein(fields=self.node_fields))for _ in range(100)]

        #if ExampleAdapterNodeType.DISEASE in self.node_types:
         #   [self.nodes.append(Disease(fields=self.node_fields)) for _ in range(100)]

        for node in self.nodes:
            yield (node.get_id(), node.get_label(), node.get_properties())

    def get_edges(self, probability: float = 0.3):
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
                    and ExampleAdapterEdgeType.PROTEIN_PROTEIN_INTERACTION
                    in self.edge_types
                ):
                    edge_type = ExampleAdapterEdgeType.PROTEIN_PROTEIN_INTERACTION.value
                elif (
                    isinstance(other_node, Disease)
                    and ExampleAdapterEdgeType.PROTEIN_DISEASE_ASSOCIATION
                    in self.edge_types
                ):
                    edge_type = ExampleAdapterEdgeType.PROTEIN_DISEASE_ASSOCIATION.value
                else:
                    continue

                yield (
                    relationship_id,
                    node.get_id(),
                    other_node.get_id(),
                    edge_type,
                    {"example_property": "example_value"}, #### ????
                )

    def get_node_count(self):
        """
        Returns the number of nodes generated by the adapter.
        """
        return len(list(self.get_nodes()))

    def _set_types_and_fields(self, node_types, node_fields, edge_types, edge_fields):
        if node_types:
            self.node_types = node_types
        else:
            self.node_types = [type for type in ExampleAdapterNodeType]

        if node_fields:
            self.node_fields = node_fields
        else:
            self.node_fields = [
                field
                for field in chain(
                    ExampleAdapterProteinField,
                    ExampleAdapterDiseaseField,
                )
            ]

        if edge_types:
            self.edge_types = edge_types
        else:
            self.edge_types = [type for type in ExampleAdapterEdgeType]

        if edge_fields:
            self.edge_fields = edge_fields
        else:
            self.edge_fields = [field for field in chain()]


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

    def __init__(self, fields: Optional[list] = None):
        self.fields = fields
        self.id = self._generate_id()
        self.label = "uniprot_protein"
        self.properties = self._generate_properties()

    def _generate_id(self):
        """
        Generate a random UniProt-style id.
        """
        lets = [random.choice(string.ascii_uppercase) for _ in range(3)]
        nums = [random.choice(string.digits) for _ in range(3)]

        # join alternating between lets and nums
        return "".join([x for y in zip(lets, nums) for x in y])

    def _generate_properties(self):
        properties = {}

        ## random amino acid sequence
        if (
            self.fields is not None
            and ExampleAdapterProteinField.SEQUENCE in self.fields
        ):
            # random int between 50 and 250
            l = random.randint(50, 250)

            properties["sequence"] = "".join(
                [random.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(l)],
            )

        ## random description
        if (
            self.fields is not None
            and ExampleAdapterProteinField.DESCRIPTION in self.fields
        ):
            properties["description"] = " ".join(
                [random.choice(string.ascii_lowercase) for _ in range(10)],
            )

        ## taxon
        if self.fields is not None and ExampleAdapterProteinField.TAXON in self.fields:
            properties["taxon"] = irefindex_species

        return properties


class Disease(Node):
    """
    Generates instances of diseases.
    """

    def __init__(self, fields: Optional[list] = None):
        self.fields = fields
        self.id = self._generate_id()
        self.label = "do_disease"
        self.properties = self._generate_properties()

    def _generate_id(self):
        """
        Generate a random disease id.
        """
        nums = [random.choice(string.digits) for _ in range(8)]

        return f"DOID:{''.join(nums)}"

    def _generate_properties(self):
        properties = {}

        ## random name
        if self.fields is not None and ExampleAdapterDiseaseField.NAME in self.fields:
            properties["name"] = " ".join(
                [random.choice(string.ascii_lowercase) for _ in range(10)],
            )

        ## random description
        if (
            self.fields is not None
            and ExampleAdapterDiseaseField.DESCRIPTION in self.fields
        ):
            properties["description"] = " ".join(
                [random.choice(string.ascii_lowercase) for _ in range(10)],
            )

        return properties


    

