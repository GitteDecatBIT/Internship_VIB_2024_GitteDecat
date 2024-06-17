#!/usr/bin/env python
from biocypher._logger import logger
from enum import Enum, auto
from typing import Union
from itertools import chain
import re
from time import time
from typing import Optional
from bioregistry import normalize_curie
from tqdm import tqdm  # progress bar
from dataclasses import dataclass


logger.debug(f"Loading module {__name__}.")

@dataclass
class Interaction:
    partner_a: str
    partner_b: str
    pmid: set[str]
    method: set[str]
    taxon_a: set[str]
    taxon_b: set[str]
    taxon_id: set[str]
    relationship_id: set[str]

class IRefIndexNodeType(Enum):
    PROTEIN = auto()

class IRefIndexNodeFields(Enum):
    PUBMED_IDS = "pmid"
    TAXON = "taxon"
    METHODS = "method"

class IRefIndexEdgeType(Enum):
    PROTEIN_PROTEIN_INTERACTION = "protein_protein_interaction"

class IRefIndexEdgeFields(Enum):
    RELATIONSHIP_ID = "relationship_id"


# Function to determine if an ID is a UniProt/Refseq/Entrez ID
    
def find_protein_type(node_id: str) -> Optional[str]:
    if bool(re.match(
            r"^([A-N,R-Z][0-9]([A-Z][A-Z, 0-9][A-Z, 0-9][0-9]){1,2})|([O,P,Q][0-9][A-Z, 0-9][A-Z, 0-9][A-Z, 0-9][0-9])(\.\d+)?(\-\d+)?$",
            node_id,
            )
            ):
        return "uniprot"

    elif bool(
                re.match(
                    r"^(((AC|AP|NC|NG|NM|NP|NR|NT|NW|WP|XM|XP|XR|YP|ZP)_\d+)|(NZ_[A-Z]{2,4}\d+))(\.\d+)?$",
                    node_id,
                )
            ):
        return "refseq"

    elif bool(re.match(r"^[A-Z]+[0-9]+(\.\d+)?$", node_id)):
        return "entrez"

    return None

class IRefIndexAdapter:
    """
    Adapter that creaes nodes and edges for creating a knowledge graph.

    Args:
        node_types: List of node types to include in the result.
        node_fields: List of node fields to include in the result.
        edge_types: List of edge types to include in the result.
        edge_fields: List of edge fields to include in the result.
    """

    def __init__(
        self,
        output_dir=None,
        irefindex_fields: Union[
            None, list[IRefIndexNodeFields, IRefIndexEdgeFields]
        ] = None,
        add_prefix=True,
        aggregate_pubmed_ids: bool = True,
        aggregate_methods: bool = True,
        nodes_ids=None,
        node_types: Union[None, list[IRefIndexNodeType]] = None,
        node_fields: Union[None, list[IRefIndexNodeFields]] = None,
        edge_types: Union[None, list[IRefIndexEdgeType]] = None,
        edge_fields: Union[None, list[IRefIndexEdgeFields]] = None,
    ):
        self.output_dir = output_dir
        self.irefindex_fields = irefindex_fields
        self.add_prefix = add_prefix
        self.nodes_ids = nodes_ids
        self.add_prefix = add_prefix
        self.aggregate_dict = {
            IRefIndexNodeFields.PUBMED_IDS.value: aggregate_pubmed_ids,
            IRefIndexNodeFields.METHODS.value: aggregate_methods,
        }
        self.interactions = {}

        self._set_types_and_fields(node_types, node_fields, edge_types, edge_fields)

    def _set_types_and_fields(self, node_types, node_fields, edge_types, edge_fields):
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

    def irefindex_process(self,taxon_id, paths):
        logger.info("Extracting information from IRefIndex data")
        # Get input file
        inputfile= paths[0]
        logger.info("This is the input file that is used: {}".format(inputfile))

        # Open the file
        with open(inputfile, "r") as file:
            logger.info(
                "Getting information for partner_a, partner_b, pubmed_id, method, taxon id and relationship id"
            )

            for line in file:
                # Split the line by tab character
                line = line.split("\t")
            
                # PARTNER_A: finalReference A
                input_partner_a = line[38]
                # PARTNER_B: FinalReference B
                input_partner_b = line[39]

                # skip lines that start with complex
                if input_partner_a.startswith("complex:") or input_partner_b.startswith(
                    "complex:"
                ):
                    continue

                if input_partner_a.startswith("pdb:") or input_partner_b.startswith(
                    "pdb:"
                ):
                    continue

                if input_partner_a.startswith("flybase:") or input_partner_b.startswith(
                    "flybase:"
                ):
                    continue

                # isolate id for partner_a
                parts = input_partner_a.split(":")  
                if len(parts) > 1:
                    partner_a = parts[1]  
                else:
                    partner_a = ""

                # isolate id for partner_b
                parts = input_partner_b.split(":")
                if len(parts) > 1:
                    partner_b = parts[1]
                else:
                    partner_b = ""

                # PUBMED_ID
                input_pmid = line[8]
                parts = input_pmid.split("|")
                last_part = parts[-1]
                numbers = last_part.split(":")
                pmid = numbers[-1]

                # METHOD
                input_method = line[6]
                pattern_method = r"\((.*?)\)"
                match_method = re.search(pattern_method, input_method)
                if match_method:
                    method = match_method.group(
                        1
                    )  # Extracting the text between brackets
                else:
                    method = ""

                # taxon_a
                input_taxon = line[9]
                pattern_taxon = r"taxid:(\d+)"
                match_taxon = re.search(pattern_taxon, input_taxon)
                if match_taxon:
                    taxon_a = match_taxon.group(1)
                else:
                    taxon_a = ""
                
                # taxon_b
                input_taxon = line[10]
                pattern_taxon = r"taxid:(\d+)"
                match_taxon = re.search(pattern_taxon, input_taxon)
                if match_taxon:
                    taxon_b = match_taxon.group(1)
                else:
                    taxon_b = ""

                # relationship id
                input_relationhsip_id = line[13]
                pattern_relationship_id = r"rigid:([^|]+)"
                match_relationship_id = re.search(
                    pattern_relationship_id, input_relationhsip_id
                )

                # Extract the matched value if found
                if match_relationship_id:
                    relationship_id = match_relationship_id.group(1)

                else:
                    relationship_id = ""

                if taxon_id != "*":
                    if (taxon_a, taxon_b) != (taxon_id, taxon_id):
                        continue
                

                if (partner_a, partner_b) not in self.interactions:
                    self.interactions[(partner_a, partner_b)] = Interaction(
                        partner_a=partner_a,
                        partner_b=partner_b,
                        pmid={pmid} if pmid.strip() != "" else set(),
                        method={method} if method != "" else set(),
                        taxon_a={taxon_a} if taxon_a != "" else set(),
                        taxon_b={taxon_b} if taxon_b != "" else set(),
                        taxon_id = {taxon_id},
                        relationship_id={relationship_id} if relationship_id != "" else set(),
                    )
                else:
                    self.interactions[(partner_a, partner_b)].pmid.add(pmid)
                    if method != "":
                        self.interactions[(partner_a, partner_b)].method.add(method)
                    if taxon_a != "":
                        self.interactions[(partner_a, partner_b)].taxon_a.add(taxon_a)
                    if taxon_b != "":
                        self.interactions[(partner_a, partner_b)].taxon_b.add(taxon_b)
                    if taxon_id != "":
                        self.interactions[(partner_a, partner_b)].taxon_id.add(taxon_id)
                    if relationship_id != "":
                        self.interactions[(partner_a, partner_b)].relationship_id.add(relationship_id)

            logger.info(
                "--> Succesfully extracted information from the IRefIndex database!"
            )

            return self.interactions

    def get_nodes(self):
        """
        Processor function for irefindex data. It drops duplicate and reciprocal duplicate protein pairs and collects pubmed ids of duplicated pairs. In addition, it
        maps entries to uniprot ids using gene name and tax id information in the BioGRID data. Also, it filters protein pairs found in swissprot.

         Args:
            rename_selected_fields : List of new field names for selected fields. If not defined, default field names will be used.
        """
        logger.info("Generating nodes.")

        t1 = time()

        node_ids = set()
        self.nodes= []
        for (partners, interaction) in self.interactions.items():
            for node_id in partners:
                if node_id not in node_ids:
                    node_ids.add(node_id)

                    if node_id == interaction.partner_a:
                        taxon_id = interaction.taxon_a
                    if node_id == interaction.partner_b:
                        taxon_id = interaction.taxon_b

                    # for node_id in node_ids:
                    #taxon_id = interaction.taxon_id
                    pubmed_id = interaction.pmid
                    method = interaction.method
                    protein_type = f"{find_protein_type(node_id)}_protein"
                    
                    node_id = self.add_prefix_to_identifier(node_id)

                    yield (
                        node_id, protein_type, {
                            "pubmed_ids": "|".join(pubmed_id),
                            "taxon_id": "|".join(taxon_id),
                            "method": "|".join(method),
                        }
                    )
        t2 = time()
        logger.info(f"Nodes were generated in {round((t2-t1) / 60, 2)} mins")
            
    
    def add_prefix_to_identifier(self, identifier=None, sep=":") -> str:
        """
        Adds prefix to the protein id based on the identifier type
        """

        if self.add_prefix and identifier:
            prefix = find_protein_type(identifier) or ""

            return normalize_curie(prefix + sep + str(identifier))

        return identifier

    def get_edges(self, label: str = "protein_protein_interaction"):
        """
        Get  edges from merged data
        Args:
            label: label of protein-protein interaction edges --> Must be same as in the schema_config.yaml file
        """
        logger.info("Generating edges.")
        t3 = time()


        # create edge list
        for row in tqdm(
            self.interactions.values()
        ):
            _source = self.add_prefix_to_identifier(identifier= row.partner_a)
            _target = self.add_prefix_to_identifier(identifier=row.partner_b)

            props = {
                "pubmed_ids": [val.replace("'", "^") for val in row.pmid],
                "method": [val.replace("'", "^") for val in row.method],
                "taxon_a": [val.replace("'", "^") for val in row.taxon_a],
                "taxon_b": [val.replace("'", "^") for val in row.taxon_b],
                "relationship_id": [val.replace("'", "^") for val in row.relationship_id],
            }

            for m in props.get("method", []):
                assert "|" not in m

            yield (None, _source, _target, label, props)
        t4 = time()
        logger.info(f"Edges were generated in {round((t4-t3) / 60, 2)} mins")



