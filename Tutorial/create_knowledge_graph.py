#!/usr/bin/env python
from biocypher import BioCypher, Resource
from adapter.IRefIndex_adapter import IRefIndexAdapter,IRefIndexNodeType, IRefIndexEdgeType, IRefIndexEdgeFields
from adapter.Uniprot_adapter import Uniprot,UniprotNodeType,UniprotNodeField,UniprotEdgeType,UniprotIDField
from biocypher._logger import logger
import yaml
import os



# Instantiate the BioCypher interface
# You can use `config/biocypher_config.yaml` to configure the framework or
# supply settings via parameters below
bc = BioCypher(biocypher_config_path= r"config/biocypher_config.yaml",
               schema_config_path= r"config/schema_config.yaml")

# Whether to cache data by pypath for future usage
CACHE = True

# Flag for exporting node and edge files as csv format
export_as_csv = True

# dirs
output_dir_path = "YOUR_PATH"


# uniprot configuration
uniprot_node_types = [
    UniprotNodeType.PROTEIN,
    #UniprotNodeType.GENE,
    UniprotNodeType.ORGANISM,
]

uniprot_node_fields = [
    UniprotNodeField.PRIMARY_GENE_NAME,
    UniprotNodeField.LENGTH,
    UniprotNodeField.MASS,
    UniprotNodeField.ORGANISM,
    UniprotNodeField.ORGANISM_ID,
    UniprotNodeField.PROTEIN_NAMES,
    UniprotNodeField.PROTEOME,
    UniprotNodeField.PROTEIN_GENE_NAMES,
    UniprotNodeField.ENSEMBL_TRANSCRIPT_IDS,
    UniprotNodeField.ENSEMBL_GENE_IDS,
    UniprotNodeField.ENTREZ_GENE_IDS,
    #UniprotNodeField.VIRUS_HOSTS,
    #UniprotNodeField.KEGG_IDS,
    UniprotNodeField.SEQUENCE,
    UniprotNodeField.PROTT5_EMBEDDING,
    UniprotNodeField.ESM2_EMBEDDING,
]

uniprot_edge_types = [
    UniprotEdgeType.PROTEIN_TO_ORGANISM,
    UniprotEdgeType.GENE_TO_PROTEIN,
]

uniprot_id_type = [
    UniprotIDField.GENE_ENTREZ_ID,
]

uniprot_adapter = Uniprot(
        organism="7227",
        node_types=uniprot_node_types,
        node_fields=uniprot_node_fields,
        edge_types=uniprot_edge_types,
        id_fields=uniprot_id_type,
        test_mode=False,
    )

uniprot_adapter.download_uniprot_data(cache=CACHE, retries=6)
uniprot_adapter.download_prott5_embeddings()
uniprot_adapter.retrieve_esm2_embeddings()
uniprot_nodes = uniprot_adapter.get_nodes()
uniprot_edges = uniprot_adapter.get_edges()

bc.write_nodes(uniprot_nodes)
bc.write_edges(uniprot_edges)
bc.write_import_call()

if export_as_csv:
    uniprot_adapter.export_data_to_csv(path=output_dir_path,
                                    node_data=uniprot_nodes,
                                    edge_data=uniprot_edges)
"""
############# IREFINDEX #############
############# Download and cache resources #############

# path to the config file
config_file_path = "/home/guest/Github/Internship_VIB_2024_GitteDecat/Create_IRefIndex_adapter_based_on_tutorial/config/biocypher_config.yaml"
# does not work for others with absolute path 

# Check if the config file exists
if not os.path.exists(config_file_path):
    print("Config file not found!")
    exit()

# Read the config file
with open(config_file_path, "r") as file:
    config = yaml.safe_load(file)

# Extract the URL from the config
url = config["biocypher"]["url"]

# Check if URL is present in the config
if not url:
    print("URL not found in the config!")
    exit()

logger.info("This is the link of IRefIndex data that is downloaded:{}, pleae check if it is up to date".format(url))


############# Specify resource #############
resource = Resource(
    name="IRefIndex",  # Name of the resource
    url_s= url,
    lifetime=7,  # seven days cache lifetime
)

logger.info("The url is correctly specified in the resource: {}" .format(resource))


paths = bc.download(resource)  # Downloads to '.cache' by default
logger.info("Data is downloaded in this path: {}" .format(paths))
# You can use the list of paths returned to read the resource into your adapter

# Choose node types to include in the knowledge graph.
# These are defined in the adapter (`adapter.py`).

node_types = IRefIndexNodeType.PROTEIN,
edge_type = IRefIndexEdgeType.PROTEIN_PROTEIN_INTERACTION
edge_fields= [IRefIndexEdgeFields.PUBMED_IDS,
              IRefIndexEdgeFields.TAXON, 
              IRefIndexEdgeFields.METHODS]
############# Create a protein adapter instance #############


irefindex_adapter = IRefIndexAdapter(
    node_types=node_types,
    edge_types=edge_type,
    edge_fields=edge_fields, 
    )

irefindex_adapter.irefindex_process()
irefindex_adapter.set_edge_fields()
irefindex_adapter.add_prefix_to_id()

############# Create a knowledge graph from the adapter #############
bc.write_nodes(irefindex_adapter.get_nodes())
bc.write_edges(irefindex_adapter.get_edges())



"""
############# UNIPROT #############


############# Write admin import statement #############
bc.write_import_call()

############# Print summary #############
bc.summary()

