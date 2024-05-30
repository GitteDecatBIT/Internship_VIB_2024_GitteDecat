#!/usr/bin/env python
from biocypher import BioCypher, Resource
from template_package.adapters.IRefIndex_adapter_13 import IRefIndexAdapter,IRefIndexNodeType, IRefIndexEdgeType, IRefIndexEdgeFields
from biocypher._logger import logger
import yaml
import os

#ruff



# Instantiate the BioCypher interface
# You can use `config/biocypher_config.yaml` to configure the framework or
# supply settings via parameters below
bc = BioCypher(biocypher_config_path= r"config/biocypher_config.yaml",
               schema_config_path= r"config/schema_config.yaml")


############# Download and cache resources #############
# specify the taxon_id
"""
all : *
homo sapiens: 9606
Mus musculus: 10090
Saccharomyces cerevisiae S288C: 559292
Escherichia: 562
Rattus norvegicus: 10116
Saccharomyces cerevisiae:4932
Drosophila melanogaster: 7227
Caenorhabditis elegans: 6239)

"""
taxon_id= "9606"

# path to the config file
config_file_path = os.path.join(os.getcwd(), 'config', 'biocypher_config.yaml')

# Check if the config file exists
if not os.path.exists(config_file_path):
    print("Config file not found!")
    exit()

# Read the config file
with open(config_file_path, "r") as file:
    config = yaml.safe_load(file)

if taxon_id == "*":
    url= config["biocypher"]["url"]["all"]
elif taxon_id == "9606":
    url= config["biocypher"]["url"]["homo_sapiens"]
elif taxon_id == "10090":
    url= config["biocypher"]["url"]["mus_musculus"]
elif taxon_id == "559292":
    url= config["biocypher"]["url"]["saccharomyces_cerevisiae_S288C"]
elif taxon_id == "562":
    url= config["biocypher"]["url"]["escherichia"]
elif taxon_id == "10116":
    url= config["biocypher"]["url"]["rattus_norvegicus"]
elif taxon_id == "4932":
    url= config["biocypher"]["url"]["saccharomyces_cerevisiae"]
elif taxon_id == "7227":
    url= config["biocypher"]["url"]["drosophila_melanogaster"]
elif taxon_id == "6239":
    url= config["biocypher"]["url"]["caenorhabditis_elegans"]
else: 
    print("Taxon ID not recognized!")

# Check if URL is present in the config
if not url:
    print("URL not found for the specified organism in the config!")
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
adapter = IRefIndexAdapter(
    node_types=node_types,
    edge_types=edge_type,
    edge_fields=edge_fields, 
    )

adapter.irefindex_process()
adapter.set_edge_fields()
adapter.add_prefix_to_id()

############# Create a knowledge graph from the adapter #############
bc.write_nodes(adapter.get_nodes())
bc.write_edges(adapter.get_edges())

############# Write admin import statement #############
bc.write_import_call()

############# Print summary #############
bc.summary()




#urls = "https://storage.googleapis.com/irefindex-data/archive/release_20.0/psi_mitab/MITAB2.6/7227.mitab.08-28-2023.txt.zip"