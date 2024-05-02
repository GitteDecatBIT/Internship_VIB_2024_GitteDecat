from IRefIndex.Irefindex_adapter import (
    IRefIndex
)
# PPI
irefindex_adapter = IRefIndex(organism=None, 
                                cache=CACHE, 
                                output_dir=output_dir_path,
                                export_csv=export_as_csv)

irefindex_adapter.download_irefindex_data()


irefindex_adapter.irefindex_process()


bc.write_edges(irefindex_adapter.get_ppi_edges())