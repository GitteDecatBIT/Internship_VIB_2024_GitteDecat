from BIOGRID.ppi_adapter import (
    PPI
)


from biocypher import BioCypher

bc = BioCypher(biocypher_config_path= r"config/biocypher_config.yaml",
               schema_config_path= r"config/schema_config.yaml",
)

# Whether to cache data by pypath for future usage
CACHE = True

# Flag for exporting node and edge files as csv format
export_as_csv = True

# dirs
output_dir_path = "YOUR_PATH"

# PPI
ppi_adapter = PPI(organism=None, 
                  cache=CACHE, 
                  output_dir=output_dir_path,
                  export_csv=export_as_csv)

ppi_adapter.download_intact_data()
ppi_adapter.download_biogrid_data()
ppi_adapter.download_string_data()

ppi_adapter.intact_process()
ppi_adapter.biogrid_process()
ppi_adapter.string_process()

bc.write_edges(ppi_adapter.get_ppi_edges())



# Write import call and other post-processing
bc.write_import_call()
bc.summary()
