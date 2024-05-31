#!/usr/bin/env python
from biocypher import BioCypher, Resource
from template_package.adapters.IRefIndex_adapter_14 import IRefIndexAdapter,IRefIndexNodeType, IRefIndexNodeFields, IRefIndexEdgeType, IRefIndexEdgeFields
from biocypher._logger import logger


# Define the taxon_id and the release version 
#This can be another taxon_id then in the list, if it is not in the list, 
# a file will be downloaded that contains alle the organisms and will be filtered on the specified taxon_id

"""
Create a url that will download data for a given organism
all organisms: All
homo sapiens: 9606
Mus musculus: 10090
Saccharomyces cerevisiae S288C: 559292
Escherichia: 562
Rattus norvegicus: 10116
Saccharomyces cerevisiae:4932
Drosophila melanogaster: 7227
Caenorhabditis elegans: 6239

"""
release_version = "release_20.0"
taxon_id = "9771"

taxon_ids = {
    'homo sapiens': '9606',
    'Mus musculus': '10090',    
    'Saccharomyces cerevisiae S288C': '559292',
    'Escherichia': '562',
    'Rattus norvegicus': '10116',
    'Saccharomyces cerevisiae': '4932',
    'Drosophila melanogaster': '7227',
    'Caenorhabditis elegans': '6239'
}

############# Download data #############

if taxon_id not in taxon_ids.values():
    logger.info(f"Taxon ID {taxon_id} not recognized. Downloading file for 'all organisms'.")
    url = "https://storage.googleapis.com/irefindex-data/archive/{}/psi_mitab/MITAB2.6/All.mitab.08-28-2023.txt.zip".format(release_version)
else:
    logger.info(f"Downloading file for taxon ID: {taxon_id}")
    url = "https://storage.googleapis.com/irefindex-data/archive/{}/psi_mitab/MITAB2.6/{}.mitab.08-28-2023.txt.zip".format(release_version, taxon_id)

logger.info("This is the link of IRefIndex data that is downloaded:{}".format(url))


# Main execution part (only runs when script1.py is executed directly)
if __name__ == "__main__":
    # Instantiate the BioCypher interface
    # You can use `config/biocypher_config.yaml` to configure the framework or
    # supply settings via parameters below
    bc = BioCypher(biocypher_config_path= r"config/biocypher_config.yaml",
                   schema_config_path= r"config/schema_config.yaml")


    ############# Specify resource #############
    resource = Resource(
        name="IRefIndex",  # Name of the resource
        url_s=url,
        lifetime=0,  # days cache lifetime
    )
    # PROBLEM: when setting lifetime to 0 it always downloads it again

    paths = bc.download(resource)  # Downloads to '.cache' by default
    logger.info("Data is downloaded in this path: {}".format(paths))

    # Choose node types to include in the knowledge graph.
    # These are defined in the adapter file.

    node_types = IRefIndexNodeType.PROTEIN,
    node_fields = [IRefIndexNodeFields.PUBMED_IDS,
                   IRefIndexNodeFields.TAXON,
                   IRefIndexNodeFields.METHODS]
    edge_type = IRefIndexEdgeType.PROTEIN_PROTEIN_INTERACTION
    edge_fields = IRefIndexEdgeFields.RELATIONSHIP_ID

    ############# Create an adapter instance #############
    adapter = IRefIndexAdapter(
        node_types=node_types,
        node_fields=node_fields,
        edge_types=edge_type,
        edge_fields=edge_fields,
    )

    adapter.irefindex_process()
    adapter.set_node_fields()
    adapter.add_prefix_to_id()

    ############# Create a knowledge graph from the adapter #############
    
    bc.write_nodes(adapter.get_nodes())
    bc.write_edges(adapter.get_edges())

    ############# Write admin import statement #############
    bc.write_import_call()

    ############# Print summary #############
    bc.summary()

def get_taxon_id():
    return taxon_id

