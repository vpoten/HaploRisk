import os
import datetime
from classes.ensembl_client import EnsemblRestClient
from classes.gene_database import GeneDatabase
import scipy.stats as stats


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
    dist = wsize / 2
    positions = region.get(chr, [])
    for pos in positions:
        if abs(pos - start) < dist or abs(pos - end) < dist:
            return pos
    return None


def count_genes_in_region(genes, genes_db, region, wsize):
    count = 0
    for gene_id in genes:
        gene_data = genes_db.get_by_id(gene_id)
        if is_in_region(region, gene_data.chr, gene_data.start, gene_data.end, wsize):
            count += 1
    return count


if __name__ == "__main__":
    print 'Started:', datetime.datetime.now().isoformat()

    base_path = '/home/victor/Escritorio/Genotipado_Alternativo/colocalizacion'

    taxonomy_id = '9606'
    snps_ids = load_lines(os.path.join(base_path, 'MS.txt'))
    gene_ids = load_lines(os.path.join(base_path, 'sp140_genes.txt'))

    client = EnsemblRestClient()
    snps = client.get_snps(snps_ids)
    regions = get_regions_from_ensembl_snps(snps)

    genes_db = GeneDatabase(taxonomy_id)
    genes_db.load_mart_export(os.path.join(base_path, 'GRCh38/mart_export.txt.gz'))

    control_gene_ids = genes_db.get_difference(gene_ids)

    wsize = 500000
    b1 = count_genes_in_region(gene_ids, genes_db, regions.get('GRCh38'), wsize)
    n1 = len(gene_ids)
    b2 = count_genes_in_region(control_gene_ids, genes_db, regions.get('GRCh38'), wsize)
    n2 = len(control_gene_ids)
    print 'b: %i, n: %i' % (b1, n1)
    print 'B: %i, N: %i' % (b2, n2)
    oddsratio, pvalue = stats.fisher_exact([[b1, n1], [b2, n2]])
    print 'oddsratio: %f, pvalue: %f' % (oddsratio, pvalue)

    print 'Finished:', datetime.datetime.now().isoformat()
