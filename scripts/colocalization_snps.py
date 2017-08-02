import os
import datetime
from classes.ensembl_client import EnsemblRestClient
from classes.gene_database import GeneDatabase
from classes.enrichr import EnrichR
import scipy.stats as stats
from statsmodels.sandbox.stats.multicomp import multipletests


def load_lines(file_name):
    # read lines from file and returns them in a list
    f = open(file_name, 'r')
    ids = [line.strip() for line in f]
    f.close()
    return ids


def get_regions_from_ensembl_snps(snps):
    # returns regions in this format: {'<assembly1>': {'<chrNN>': [<pos1>, <pos2>, ...], ...}, ...}
    regions = {}

    for snp_id in snps:
        data = snps[snp_id]
        mappings = data.get('mappings', [])
        for mapping in mappings:
            pos = mapping['start']
            chr = str(mapping['seq_region_name'])
            assembly = mapping['assembly_name']
            reg_assembly = regions.get(assembly, {})  # create a new one if not exists
            regions[assembly] = reg_assembly  # assign to ensure creation
            reg_chr = reg_assembly.get(chr, [])  # create a new one if not exists
            reg_assembly[chr] = reg_chr  # assign to ensure creation
            reg_chr.append(pos)

    return regions


def is_in_region(region, chr, start, end, wsize):
    # checks whether a feature located at chr:start-end overlaps a region
    dist = wsize / 2
    positions = region.get(chr, [])
    for pos in positions:
        if abs(pos - start) < dist or abs(pos - end) < dist:
            return pos
    return None


def count_genes_in_region(genes, genes_db, region, wsize, add_to_list=None):
    # count genes in list that overlaps a region
    count = 0
    for gene_id in genes:
        gene_data = genes_db.get_by_id(gene_id)
        if is_in_region(region, gene_data.chr, gene_data.start, gene_data.end, wsize) is not None:
            count += 1
            if type(add_to_list) is list:
                add_to_list.append(gene_data)
    return count


def create_snp_regions(snps_ids):
    # build snps regions using ensembl REST API
    client = EnsemblRestClient()
    snps = client.get_snps(snps_ids)
    return get_regions_from_ensembl_snps(snps)


def create_gene_db(taxonomy_id, mart_file):
    # create a gene DB using an ensembl mart export
    genes_db = GeneDatabase(taxonomy_id)
    genes_db.load_mart_export(mart_file)
    return genes_db


def calc_genes_in_region_table(region, genes_db, gene_ids, wsize):
    """
    Calculates 2x2 table of gene counts: [[b1, n1], [b2, n2]]
    b1: number of selected genes that match a region, n1: total of selected genes
    b2: number of background genes that match a region, n2: total of background genes
    :param region: dict of {chr: [p1, p2, ...]} as returned by 'get_regions_from_ensembl_snps'
    :param genes_db:
    :param gene_ids: list of selected genes (for example: diff expressed in a GEO study)
    :param wsize: window size used to calculate region match (centered in a region)
    :return: a tuple with the contingency table (2x2) and matching genes array
    """
    control_gene_ids = genes_db.get_difference(gene_ids)
    matching_genes = []
    b1 = count_genes_in_region(gene_ids, genes_db, region, wsize, add_to_list=matching_genes)
    n1 = len(gene_ids)
    b2 = count_genes_in_region(control_gene_ids, genes_db, region, wsize)
    n2 = len(control_gene_ids)
    table = [[b1, n1], [b2, n2]]
    return table, matching_genes


def human_format(num):
    # returns an amount in a human readable format
    magnitude = 0
    while abs(num) >= 1000:
        magnitude += 1
        num /= 1000.0
        # add more suffixes if you need them
    return '%i%s' % (num, ['', 'K', 'M', 'G', 'T', 'P'][magnitude])


def test1():
    base_path = '/home/victor/Escritorio/Genotipado_Alternativo/colocalizacion'

    snps_ids = load_lines(os.path.join(base_path, 'MS.txt'))
    regions = create_snp_regions(snps_ids)
    genes_db = create_gene_db('9606', os.path.join(base_path, 'GRCh38/mart_export.txt.gz'))

    gene_ids = load_lines(os.path.join(base_path, 'sp140_genes.txt'))
    wsize = 500000

    table, match_genes = calc_genes_in_region_table(regions.get('GRCh38'), genes_db, gene_ids, wsize)
    oddsratio, pvalue = stats.fisher_exact(table)

    print 'b: %i, n: %i' % tuple(table[0])
    print 'B: %i, N: %i' % tuple(table[1])
    print 'oddsratio: %f, pvalue: %f' % (oddsratio, pvalue)


def enrichr_db_test(file_name, region, genes_db, wsize, pvalue_thr=0.05):
    """
    Runs, for each record of a 'enrich_db table', a fisher test using 'calc_genes_in_region_table' contingency table
    :param file_name: full path to enrich_db table
    :param region: dict of {chr: [p1, p2, ...]} as returned by 'get_regions_from_ensembl_snps'
    :param genes_db:
    :param wsize: window size used to calculate region match (centered in a region)
    :param pvalue_thr: pvalue threshold for FDR
    :return: a list of tuples (lib_name, record, b1, n1, b2, n2, oddsratio, pval, corr_pval, matching_genes)
    """
    data_lib = EnrichR.load_library(file_name)
    lib_name = os.path.basename(file_name[:-7])  # remove '.txt.gz'
    results = []

    for record in data_lib:
        # extract genes highlighted by the study
        gene_names = EnrichR.extract_gene_list(data_lib, record)
        gene_ids = []

        for name in gene_names:
            # translate gene names to gene_id (ensembl ids)
            gene = genes_db.get_by_name(name)
            if gene:
                gene_ids.append(gene.id)
            else:
                pass  # gene name not found

        # calculate contingency table for the genes in the study
        t, match_genes = calc_genes_in_region_table(region, genes_db, gene_ids, wsize)
        oddsratio, pvalue = stats.fisher_exact(t)

        results.append((lib_name, record, t[0][0], t[0][1], t[1][0], t[1][1], oddsratio, pvalue, match_genes))

    # multiple test correction using FDR
    pvals = map(lambda r: r[7], results)
    vals = multipletests(pvals, alpha=pvalue_thr, method='fdr_bh')

    # add corrected p-value to results and filter by FDR threshold
    results_corr = []
    for i, res in enumerate(results):
        if vals[0][i]:
            # if test passed -> append to results_corr list and put matching_genes at the end of the tuple
            results_corr.append(res[:-1] + ((vals[1][i]).item(), res[-1],))

    return results_corr


if __name__ == "__main__":
    print 'Started:', datetime.datetime.now().isoformat()

    # test1()

    base_path = '/home/victor/Escritorio/Genotipado_Alternativo/colocalizacion'
    snps_ids = load_lines(os.path.join(base_path, 'MS.txt'))
    regions = create_snp_regions(snps_ids)
    genes_db = create_gene_db('9606', os.path.join(base_path, 'GRCh38/mart_export.txt.gz'))

    enrichr_path = os.path.join(base_path, 'enrichr')
    lib_files = EnrichR.list_libraries(enrichr_path)
    lib_files = filter(lambda n: n.startswith('Single_Gene_Perturbations_from_GEO'), lib_files)
    wsizes = [1e6, 500000.0, 250000.0, 100000.0, 50000.0, 20000.0, 10000.0]

    for wsize in wsizes:
        wsize_str = human_format(wsize)
        lib_results = {}

        for name in lib_files:
            lib_name = name[:-7]  # remove '.txt.gz'
            res = enrichr_db_test(os.path.join(enrichr_path, name), regions.get('GRCh38'), genes_db, wsize)
            print '%i matches in %s, [%s]' % (len(res), lib_name, datetime.datetime.now().isoformat())
            lib_results[lib_name] = res

        f = open(os.path.join(base_path, 'output_enrichr_%s.txt' % wsize_str), 'w')
        for lib_name in lib_results:
            for res in lib_results[lib_name]:
                f.write('%s\t%s\t%i\t%i\t%i\t%i\t%f\t%f\t%f\t' % res[:9])
                # add matching genes at the end of the row (comma separated)
                f.write('%s\n' % ','.join(map(lambda g: g.name, res[9])))
        f.close()

    print 'Finished:', datetime.datetime.now().isoformat()
