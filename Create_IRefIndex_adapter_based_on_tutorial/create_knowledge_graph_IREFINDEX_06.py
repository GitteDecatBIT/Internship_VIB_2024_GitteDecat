#!/usr/bin/env python
from biocypher import BioCypher, Resource
from template_package.adapters.IRefIndex_adapter_14 import IRefIndexAdapter,IRefIndexNodeType, IRefIndexNodeFields, IRefIndexEdgeType, IRefIndexEdgeFields
from biocypher._logger import logger

############# Download data #############
# Define the taxon_id and the release version 
# if you want to specify a taxon_id that is not in the list below, the file "All" will be downloaded and filtered on taxon_id

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
taxon_id = "7227"

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

if taxon_id not in taxon_ids.values():
    logger.info(f"Taxon ID {taxon_id} not recognized. Downloading file for 'all organisms'.")
    url = "https://storage.googleapis.com/irefindex-data/archive/{}/psi_mitab/MITAB2.6/All.mitab.08-28-2023.txt.zip".format(release_version)
else:
    logger.info(f"Downloading file for taxon ID: {taxon_id}")
    url = "https://storage.googleapis.com/irefindex-data/archive/{}/psi_mitab/MITAB2.6/{}.mitab.08-28-2023.txt.zip".format(release_version, taxon_id)

logger.info("This is the link of IRefIndex data that is downloaded:{}".format(url))


# Main execution part 
if __name__ == "__main__":
    # Instantiate the BioCypher interface
    bc = BioCypher(biocypher_config_path= r"config/biocypher_config.yaml",
                   schema_config_path= r"config/schema_config.yaml")


    ############# Specify resource #############
    resource = Resource(
        name="IRefIndex",  # Name of the resource
        url_s=url,
        lifetime=0,  
    )

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

