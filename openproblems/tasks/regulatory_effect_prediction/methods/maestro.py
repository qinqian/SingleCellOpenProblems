from ....tools.decorators import method

import pandas as pd


def _rp_simple(tss_to_peaks, adata, log_each=500):
    import numpy as np
    import scipy

    # the coordinates of the current peaks
    peaks = adata.uns["mode2_var"]

    decay = 10000

    def Sg(x):
        return 2 ** (-x)

    # gene_distance = 15 * decay

    weights = []

    tss_to_peaks = tss_to_peaks.drop_duplicates("itemRgb")

    # print(tss_to_peaks.shape)
    # print(tss_to_peaks.head())

    for ri, r in tss_to_peaks.iterrows():
        wi = 0
        summit_chr, tss_summit_start, tss_summit_end = r[:3]
        tss_extend_chr, tss_extend_start, tss_extend_end = r[4:7]

        # print(summit_chr, summit_start, summit_end,
        #       extend_chr, extend_start, extend_end)
        sel_chr = [pi for pi in peaks if pi[0] == tss_extend_chr]
        sel_peaks = [
            pi
            for pi in sel_chr
            if int(pi[1]) >= tss_extend_start and int(pi[2]) <= tss_extend_end
        ]

        # print('# peaks in chromosome', len(sel_chr),
        #       '# of peaks around tss', len(sel_peaks))
        # if len(sel_peaks) > 0:
        # print(sel_peaks)

        # if peaks then this is take them into account, one by one
        for pi in sel_peaks:
            summit_peak = int((int(pi[2]) + int(pi[1])) / 2)
            distance = np.abs(tss_summit_start - summit_peak)
            # print(pi, distance, Sg(distance / decay))
            wi += Sg(distance / decay)

        if log_each is not None and len(weights) % log_each == 0:
            if len(weights) > 0:
                print(
                    "# weights calculated so far",
                    len(weights),
                    "out of",
                    tss_to_peaks.shape[0],
                )
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


def _rp_enhanced(tss_to_peaks, adata, log_each=500):
    import numpy as np
    import scipy

    # prepare the exonic ranges
    exon_ranges = []
    for exons in adata.var["exons"]:
        exon_ranges.append([e.start, e.end] for e in exons)
    adata.var["exon.ranges"] = exon_ranges
    exon_coordinates_by_gene = adata.var["exon.ranges"].to_dict()

    tss_to_peaks["exon.ranges"] = tss_to_peaks["itemRgb"].map(exon_coordinates_by_gene)
    tss_to_peaks = tss_to_peaks.drop_duplicates("itemRgb").reset_index(drop=True)

    # the coordinates of the current peaks
    peaks = adata.uns["mode2_var"]

    decay = 10000

    def Sg(x):
        return 2 ** (-x)

    print("calculating weights per gene...")
    weights = []
    for ri, r in tss_to_peaks.iterrows():
        wi = 0
        summit_chr, tss_summit_start, tss_summit_end = r[:3]
        tss_extend_chr, tss_extend_start, tss_extend_end = r[4:7]

        # gene_name = r[-2]
        exon_ranges = r[-1]

        # print(summit_chr, tss_summit_start, tss_summit_end,
        #       tss_extend_chr, tss_extend_start, tss_extend_end,
        #      gene_name, exon_ranges)
        sel_chr = [pi for pi in peaks if pi[0] == tss_extend_chr]
        sel_peaks = [
            pi
            for pi in sel_chr
            if int(pi[1]) >= tss_extend_start and int(pi[2]) <= tss_extend_end
        ]
        # check whether the peak overlaps with a given exon
        if not pd.isnull(exon_ranges):
            sel_peak_summits = [(int(pi[1]) + int(pi[2])) / 2.0 for pi in sel_peaks]
            peak_in_exons = [
                np.sum([ps >= ex[0] and ps <= ex[1] for ex in exon_ranges]) >= 1
                for ps in sel_peak_summits
            ]
        else:
            peak_in_exons = [False for pi in sel_peaks]

        # if sum(peak_in_exons) > 0:
        #     print ('exon / peak overlap found!')
        #     print(ri, peak_in_exons)
        # if peaks then this is take them into account, one by one
        for pi, peak_in_exon in zip(sel_peaks, peak_in_exons):
            # the current peak is part of an exon
            # if peak_in_exon:
            #     print(pi)
            summit_peak = int((int(pi[2]) + int(pi[1])) / 2)
            distance = np.abs(tss_summit_start - summit_peak)
            # print(pi, distance, Sg(distance / decay))
            wi += Sg(distance / decay) if not peak_in_exon else 1.0

        if log_each is not None and len(weights) % log_each == 0:
            if len(weights) > 0:
                print(
                    "# weights calculated so far",
                    len(weights),
                    "out of",
                    tss_to_peaks.shape[0],
                )

        weights.append(wi)

    out = tss_to_peaks.copy()
    out["weight"] = weights

    gene_peak_weight = scipy.sparse.csr_matrix(
        (
            out.weight.values,
            (out.thickEnd.astype("int32").values, out.name.values),
        ),
        shape=(adata.shape[1], adata.uns["mode2_var"].shape[0]),
    )

    adata.obsm["gene_score"] = adata.obsm["mode2"] @ gene_peak_weight.T


@method(
    method_name="RP_simple",
    paper_name="""Integrative analyses of single-cell transcriptome\
and regulome using MAESTRO.""",
    paper_url="https://pubmed.ncbi.nlm.nih.gov/32767996",
    paper_year=2020,
    code_version="1.0",
    code_url="https://github.com/liulab-dfci/MAESTRO",
    image="openproblems-python-extras",
)
def rp_simple(adata, n_top_genes=2000):
    from .beta import _atac_genes_score

    adata = _atac_genes_score(
        adata,
        top_genes=n_top_genes,
        method="rp_simple",
    )
    return adata


@method(
    method_name="RP_enhanced",
    paper_name="""Integrative analyses of single-cell transcriptome\
and regulome using MAESTRO.""",
    paper_url="https://pubmed.ncbi.nlm.nih.gov/32767996",
    paper_year=2020,
    code_version="1.0",
    code_url="https://github.com/liulab-dfci/MAESTRO",
    image="openproblems-python-extras",
)
def rp_enhanced(adata, n_top_genes=2000):
    from .beta import _atac_genes_score

    adata = _atac_genes_score(adata, top_genes=n_top_genes, method="rp_enhanced")
    return adata
