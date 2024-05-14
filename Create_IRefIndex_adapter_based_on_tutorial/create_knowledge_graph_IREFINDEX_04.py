#!/usr/bin/env python
from biocypher import BioCypher, Resource
from template_package.adapters.IRefIndex_adapter_06 import IRefIndexAdapter
from biocypher._logger import logger
import yaml
import os



# Instantiate the BioCypher interface
# You can use `config/biocypher_config.yaml` to configure the framework or
# supply settings via parameters below
bc = BioCypher()


############# Download and cache resources #############

# path to the config file
config_file_path = "/home/guest/Github/Internship_VIB_2024_GitteDecat/Create_IRefIndex_adapter_based_on_tutorial/config/biocypher_config.yaml"

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

logger.info("This is the link of IRefIndex data that is downloaded:{}".format(url))


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


############# Create a protein adapter instance #############
adapter = IRefIndexAdapter()



adapter.irefindex_process()
adapter.set_edge_fields()
adapter.add_prefix_to_id()
#adapter.get_nodes()


############# Create a knowledge graph from the adapter #############
bc.write_nodes(adapter.get_nodes())
bc.write_edges(adapter.get_edges())


############# Write admin import statement #############
bc.write_import_call()

############# Print summary #############
bc.summary()




#urls = "https://storage.googleapis.com/irefindex-data/archive/release_20.0/psi_mitab/MITAB2.6/7227.mitab.08-28-2023.txt.zip"