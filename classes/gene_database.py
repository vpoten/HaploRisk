import gzip
import collections

Gene = collections.namedtuple('Gene', ['id', 'name', 'chr', 'start', 'end', 'strand', 'type'])


class GeneDatabase(object):
    valid_chrs = [str(i) for i in range(1, 23)] + ['X', 'Y']

    def __init__(self, taxonomy_id):
        """
        :param taxonomy_id:
        """
        self.taxonomy_id = taxonomy_id
        self.gene_index = {}
        self.chr_index = {}

    def load_mart_export(self, mart_file):
        """
        Load a biomart export TSV (gzipped) with format:
        Gene_stable_ID, Chromosome/scaffold_name, Gene_start_(bp), Gene_end_(bp), Strand, Gene_name, Gene_type
        :param mart_file:
        :return:
        """
        f = gzip.open(mart_file, 'r')
        f.readline()  # skip header

        for line in f:
            toks = line.strip().split('\t')
            chr = toks[1]

            if chr not in self.valid_chrs:
                continue

            gene = Gene(id=toks[0], name=toks[5], chr=chr, start=int(toks[2]), end=int(toks[3]),
                        strand=int(toks[4]), type=toks[6])
            self.gene_index[gene.id] = gene
            chr_list = self.chr_index.get(gene.chr, [])  # create if new
            self.chr_index[gene.chr] = chr_list
            chr_list.append(gene)

        f.close()

        # sort chr lists by position
        for chr in self.chr_index.keys():
            self.chr_index[chr] = sorted(self.chr_index[chr], key=lambda g: g.start)

    def get_by_id(self, id):
        return self.gene_index.get(id)

    def get_chr_genes(self, chr):
        return self.chr_index.get(chr)

    def get_difference(self, gene_ids):
        diff_set = []
        for db_id in self.gene_index:
            if db_id not in gene_ids:
                diff_set.append(db_id)
        return diff_set
