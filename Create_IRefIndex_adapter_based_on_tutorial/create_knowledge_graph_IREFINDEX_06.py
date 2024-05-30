#!/usr/bin/env python
from biocypher import BioCypher, Resource
from template_package.adapters.IRefIndex_adapter_13 import IRefIndexAdapter,IRefIndexNodeType, IRefIndexEdgeType, IRefIndexEdgeFields
from biocypher._logger import logger


#ruff

# Instantiate the BioCypher interface
# You can use `config/biocypher_config.yaml` to configure the framework or
# supply settings via parameters below
bc = BioCypher(biocypher_config_path= r"config/biocypher_config.yaml",
               schema_config_path= r"config/schema_config.yaml")


############# Download and cache resources #############

"""
Create a url that will donwload data for a given organism
all organisms: All
homo sapiens: 9606
Mus musculus: 10090
Saccharomyces cerevisiae S288C: 559292
Escherichia: 562
Rattus norvegicus: 10116
Saccharomyces cerevisiae:4932
Drosophila melanogaster: 7227
Caenorhabditis elegans: 6239)

"""
taxon_id= "All"
release_version= "release_20.0" 

url= "https://storage.googleapis.com/irefindex-data/archive/{}/psi_mitab/MITAB2.6/{}.mitab.08-28-2023.txt.zip".format(release_version, taxon_id)

logger.info("This is the link of IRefIndex data that is downloaded:{}".format(url))


############# Specify resource #############
resource = Resource(
    name="IRefIndex",  # Name of the resource
    url_s= url,
    lifetime=0,  #days cache lifetime
)
# PROBLEM: when setting lifetime to 0 it always downloads is again 

paths = bc.download(resource)  # Downloads to '.cache' by default
logger.info("Data is downloaded in this path: {}" .format(paths)) 

# Choose node types to include in the knowledge graph.
# These are defined in the adapter file.

node_types = IRefIndexNodeType.PROTEIN,
edge_type = IRefIndexEdgeType.PROTEIN_PROTEIN_INTERACTION
edge_fields= [IRefIndexEdgeFields.PUBMED_IDS,
              IRefIndexEdgeFields.TAXON, 
              IRefIndexEdgeFields.METHODS]

############# Create a adapter instance #############
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
