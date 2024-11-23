import pandas as pd

from .KarUtils import get_metadata_file_path


def get_genes_in_region(chrom, start, end, gene_file=get_metadata_file_path('gtf_protein_coding.bed')):
    """
    report all genes with intersection with the region (inclusive)
    :param chrom: chr{1-22, X, Y}
    :param start:
    :param end:
    :param gene_file: bed file containing genes and location
    :return:
    """
    gene_df = pd.read_csv(gene_file, sep='\t', header=None, names=['chrom', 'start', 'end', 'gene'])
    chrom = chrom[:3].lower() + chrom[3:]

    intersection_gene_df = gene_df[gene_df['chrom'] == chrom]
    intersection_gene_df = intersection_gene_df[intersection_gene_df['start'] <= end]
    intersection_gene_df = intersection_gene_df[intersection_gene_df['end'] >= start]
    intersection_gene_df = intersection_gene_df.sort_values(by='start')
    intersection_gene_df = intersection_gene_df.drop_duplicates(subset=['gene'], keep='first')
    return intersection_gene_df['gene'].tolist()


def get_DDG_overlapped_genes(input_gene_list, DDG_file=get_metadata_file_path('DDG2P_14_11_2023.csv')):
    DDG_df = pd.read_csv(DDG_file, sep='\t')
    overlapped_gene_df = DDG_df[DDG_df['gene symbol'].isin(input_gene_list)]
    return overlapped_gene_df


def format_gene_disease_omim(filtered_DDG_df):
    """
    :param filtered_DDG_df:
    :return: zipped gene, disease, allelic requirement, mutation consequence, confidence
    """
    # for the same gene (under same gene OMIM), we can have multiple diseases (different disease OMIM)
    col_of_interest = ['gene symbol', 'gene mim', 'disease name', 'allelic requirement', 'mutation consequence', 'confidence category', 'organ specificity list']
    vals = [{col: row[col] for col in col_of_interest} for _, row in filtered_DDG_df.iterrows()]
    return vals


def test_if_DDG_has_duplicate():
    DDG_file = 'Metadata/DDG2P_14_11_2023.csv'
    DDG_df = pd.read_csv(DDG_file, sep='\t')
    print(DDG_df.shape)
    filtered_DDG_df = DDG_df.drop_duplicates(subset=['gene symbol'], keep=False)
    print(filtered_DDG_df.shape)
    overlapping_rows = DDG_df.groupby('gene symbol').filter(lambda x: len(x) > 1)
    print(overlapping_rows)


def get_band_location(chrom, nt_idx, cyto_file=get_metadata_file_path('cytoBand.txt')):
    cyto_df = pd.read_csv(cyto_file, sep='\t', header=None, names=['chrom', 'start', 'end', 'band_name', 'stain'])
    chrom = chrom[:3].lower() + chrom[3:]
    for index, row in cyto_df.iterrows():
        if row['chrom'] != chrom:
            continue
        if row['start'] <= nt_idx <= row['end']:
            return row['band_name']
    raise RuntimeError('no band found')


########################CODE ORIGINALLY IN INTERPRETER###########################
def gather_breakpoints(breakpoints: []):
    all_breakpoints = []
    for c_breakpoint in breakpoints:
        if c_breakpoint[1] is None:
            continue
        all_breakpoints.append((c_breakpoint[0], c_breakpoint[1]))
    return all_breakpoints


def get_genes_near_breakpoints(breakpoints: [(str, int)], proximity=50000):
    breakpoint_ranges = []
    for c_breakpoint in breakpoints:
        breakpoint_ranges.append((c_breakpoint[0],
                                  c_breakpoint[1] - proximity,
                                  c_breakpoint[1] + proximity))
    genes_in_regions = []
    for breakpoint_range in breakpoint_ranges:
        genes_in_regions += get_genes_in_region(*breakpoint_range)
    return genes_in_regions


def report_on_genes_based_on_breakpoints(breakpoints):
    breakpoints = gather_breakpoints(breakpoints)
    genes_near_bp = get_genes_near_breakpoints(breakpoints)
    DDG_df = get_DDG_overlapped_genes(genes_near_bp)
    formated_DDG_genes = format_gene_disease_omim(DDG_df)
    return genes_near_bp, formated_DDG_genes


def report_cnv_genes_on_region(chrom, start, end, proximity=50000):
    if start > end:
        temp = start
        start = end
        end = temp
    start = max(0, start - proximity)
    end = end + proximity
    genes = get_genes_in_region('chr' + chrom, start, end)
    DDG_df = get_DDG_overlapped_genes(genes)
    formated_DDG_genes = format_gene_disease_omim(DDG_df)
    return genes, formated_DDG_genes


def test_get_genes():
    genes = get_genes_in_region('Chr22', 25200725, 25560371)
    print(genes)
    df = get_DDG_overlapped_genes(genes)
    a, b = format_gene_disease_omim(df)
    for idx, a_itr in enumerate(a):
        print(a_itr, b[idx])


def test_get_band():
    x = get_band_location('Chr22', 25200725)
    print(x)


if __name__ == "__main__":
    test_get_band()
