from ....tools.decorators import method


def _rp_simple(tss_to_peaks, adata):
    import numpy as np
    import scipy

    # the coordinates of the current peaks
    peaks = adata.uns["mode2_var"]

    decay = 10000

    def Sg(x): return 2 ** (-x)
    gene_distance = 15 * decay

    weights = []
    for ri, r in tss_to_peaks.iterrows():
        wi = 0
        summit_chr, tss_summit_start, tss_summit_end = r[:3]
        tss_extend_chr, tss_extend_start, tss_extend_end = r[4:7]

        # print(summit_chr, summit_start, summit_end, extend_chr, extend_start, extend_end)
        sel_chr = [pi for pi in peaks if pi[0] == tss_extend_chr]
        sel_peaks = [
            pi
            for pi in sel_chr
            if int(pi[1]) >= tss_extend_start and int(pi[2]) <= tss_extend_end
        ]

        # print('# peaks in chromosome', len(sel_chr), '# of peaks around tss', len(sel_peaks))
        # if len(sel_peaks) > 0:
        # print(sel_peaks)

        # if peaks then this is take them into account, one by one
        for pi in sel_peaks:
            summit_peak = int((int(pi[2]) + int(pi[1])) / 2)
            distance = np.abs(tss_summit_start - summit_peak)
            # print(pi, distance, Sg(distance / decay))
            wi += Sg(distance / decay)

        weights.append(wi)

    tss_to_peaks["weight"] = weights

    gene_peak_weight = scipy.sparse.csr_matrix(
        (
            tss_to_peaks.weight.values,
            (tss_to_peaks.thickEnd.astype("int32").values, tss_to_peaks.name.values),
        ),
        shape=(adata.shape[1], adata.uns["mode2_var"].shape[0]),
    )

    adata.obsm["gene_score"] = adata.obsm["mode2"] @ gene_peak_weight.T


@method(
    method_name="RP_simple",
    paper_name="""Integrative analyses of single-cell transcriptome and regulome using MAESTRO.""",
    paper_url="https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02116-x",
    paper_year=2020,
    code_version="1.0",
    code_url="https://github.com/liulab-dfci/MAESTRO",
    image="openproblems-python-extras",
)
def rp_simple(adata, n_top_genes=500):
    from .beta import _atac_genes_score

    adata = _atac_genes_score(adata, top_genes=n_top_genes, method="rp_simple")
    return adata
