from ....data.multimodal import tenx_3k_pbmc_multiome
from ....tools.decorators import dataset


@dataset("10X PBMC Multiome")
def tenx_multiome(test=False):
    adata = tenx_3k_pbmc_multiome.load_tenx_multiome(test=test)

    adata.uns["species"] = "homo_sapiens"
    adata.uns["version"] = "GRCh38"
    adata.uns["release"] = "103"
    return adata