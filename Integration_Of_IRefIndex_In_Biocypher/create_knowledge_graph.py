#!/usr/bin/env python
import sys
from biocypher import BioCypher, Resource
from template_package.adapters.IRefIndex_adapter import (
    IRefIndexAdapter,
    IRefIndexNodeType,
    IRefIndexNodeFields,
    IRefIndexEdgeType,
    IRefIndexEdgeFields,
)
from biocypher._logger import logger

############# Download data #############
# Define the taxon_id and the release version
# if you want to specify a taxon_id that is not in the list below, the file "All" will be downloaded and filtered on taxon_id

def get_taxon_id():
    """Get taxon_id from command line arguments."""
    if len(sys.argv) < 2:
        raise ValueError("Please provide a taxon_id as a command line argument.")
    return sys.argv[1]

# Main execution part
if __name__ == "__main__":
    
    release_version = "release_20.0"

    # Set of known taxon IDs
    known_taxon_ids = {
        "9606",  # homo sapiens
        "10090",  # Mus musculus
        "559292",  # Saccharomyces cerevisiae S288C
        "562",  # Escherichia
        "10116",  # Rattus norvegicus
        "4932",  # Saccharomyces cerevisiae
        "7227",  # Drosophila melanogaster
        "6239",  # Caenorhabditis elegans
    }
    taxon_id = get_taxon_id()

    if taxon_id not in known_taxon_ids:
        logger.info(
            f"Taxon ID {taxon_id} not recognized. Downloading file for 'all organisms'."
        )
        url = "https://storage.googleapis.com/irefindex-data/archive/{}/psi_mitab/MITAB2.6/All.mitab.08-28-2023.txt.zip".format(
            release_version
        )
    else:
        logger.info(f"Downloading file for taxon ID: {taxon_id}")
        url = "https://storage.googleapis.com/irefindex-data/archive/{}/psi_mitab/MITAB2.6/{}.mitab.08-28-2023.txt.zip".format(
            release_version, taxon_id
        )

    logger.info("This is the link of IRefIndex data that is downloaded:{}".format(url))

    # Instantiate the BioCypher interface
    bc = BioCypher(
        biocypher_config_path=r"config/biocypher_config.yaml",
        schema_config_path=r"config/schema_config.yaml",
    )

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
    node_types = (IRefIndexNodeType.PROTEIN,)
    node_fields = [
        IRefIndexNodeFields.PUBMED_IDS,
        IRefIndexNodeFields.TAXON,
        IRefIndexNodeFields.METHODS,
    ]
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

    ############# Create a knowledge graph from the adapter #############
    bc.write_nodes(adapter.get_nodes())
    bc.write_edges(adapter.get_edges())

    ############# Write admin import statement #############
    bc.write_import_call()

    ############# Print summary #############
    bc.summary()


def get_taxon_id():
    return taxon_id
