from biocypher import BioCypher, Resource
from template_package.adapters.IRefIndex_adapter import (
    ExampleAdapter,
    ExampleAdapterNodeType,
    ExampleAdapterEdgeType,
    ExampleAdapterProteinField,
    ExampleAdapterDiseaseField,
    IRefIndexEdgeFields,
)

# Instantiate the BioCypher interface
# You can use `config/biocypher_config.yaml` to configure the framework or
# supply settings via parameters below
bc = BioCypher()

# Download and cache resources (change the directory in the options if needed)
urls = "https://storage.googleapis.com/irefindex-data/archive/release_20.0/psi_mitab/MITAB2.6/7227.mitab.08-28-2023.txt.zip"
resource = Resource(
    name="IRefIndex",  # Name of the resource
    url_s=urls,  # URL to the resource(s)
    lifetime=7,  # seven days cache lifetime
)
paths = bc.download(resource)  # Downloads to '.cache' by default
print(paths)
# You can use the list of paths returned to read the resource into your adapter

# Choose node types to include in the knowledge graph.
# These are defined in the adapter (`adapter.py`).
node_types = [
    ExampleAdapterNodeType.PROTEIN,
]

# Choose protein adapter fields to include in the knowledge graph.
# These are defined in the adapter (`adapter.py`).
node_fields = [
    # Proteins
    ExampleAdapterProteinField.ID,
    ExampleAdapterProteinField.SEQUENCE,
    ExampleAdapterProteinField.DESCRIPTION,
    ExampleAdapterProteinField.TAXON,

]

edge_types = [
    ExampleAdapterEdgeType.PROTEIN_PROTEIN_INTERACTION,

]

edge_fields = [
    IRefIndexEdgeFields.SOURCE,
    IRefIndexEdgeFields.PUBMED_IDS,
    IRefIndexEdgeFields.METHOD,
]

# Create a protein adapter instance
adapter = ExampleAdapter(
    node_types=node_types,
    node_fields=node_fields,
    edge_types=edge_types,
    edge_fields= edge_fields,

    # we can leave edge fields empty, defaulting to all fields in the adapter
)


# Create a knowledge graph from the adapter
bc.write_nodes(adapter.get_nodes())
bc.write_edges(adapter.get_edges())

# Write admin import statement
bc.write_import_call()

# Print summary
bc.summary()
