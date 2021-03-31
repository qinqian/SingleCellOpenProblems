from .. import utils
from ...utils import loader

import os
import scanpy as sc
import scprep
import tempfile

URL = "https://ndownloader.figshare.com/files/27379250"

@loader
def load_tenx_multiome(test=False):
    """Download 10X multiome data from Figshare."""
    with tempfile.TemporaryDirectory() as tempdir:
        filepath = os.path.join(tempdir, "pbmc3kmultiome.h5ad")
        scprep.io.download.download_url(URL, filepath)
        adata = sc.read_h5ad(filepath)

    if test:
        adata = utils.subset_joint_data(adata)

    return adata
