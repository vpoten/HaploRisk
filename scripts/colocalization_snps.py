import os
import datetime
from classes.ensembl_client import EnsemblRestClient
from classes.gene_database import GeneDatabase


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


if __name__ == "__main__":
    print 'Started:', datetime.datetime.now().isoformat()

    # snps_ids = load_lines(file_name)
    taxonomy_id = '9606'
    snps_ids = ['rs12722489', 'rs6897932', 'rs6498169', 'rs6604026', 'rs10984447']

    client = EnsemblRestClient()
    snps = client.get_snps(snps_ids)
    regions = get_regions_from_ensembl_snps(snps)

    genes_db = GeneDatabase(taxonomy_id)
    genes_db.load_mart_export('/home/victor/Escritorio/Genotipado_Alternativo/colocalizacion/GRCh38/mart_export.txt.gz')

    print 'Finished:', datetime.datetime.now().isoformat()
