import gzip
import os


class EnrichR(object):
    @staticmethod
    def load_library(file_name):
        data = {}
        f = gzip.open(file_name, 'r') if file_name.endswith('.gz') else open(file_name, 'r')

        for line in f:
            toks = line.split('\t')
            record = toks[0]
            genes = []
            for gene_res in toks[2:]:
                gene_res = gene_res.strip().split(',')
                # build for each gene a tuple (gene, weight)
                genes.append((gene_res[0], 1.0 if len(gene_res) == 1 else float(gene_res[1])))
            data[record] = genes

        f.close()
        return data

    @staticmethod
    def extract_gene_list(data_lib, record):
        genes = data_lib.get(record, [])
        return map(lambda g: g[0], genes)

    @staticmethod
    def list_libraries(path):
        return filter(lambda f: os.path.isfile(os.path.join(path, f)) and f.endswith('.txt.gz'), os.listdir(path))
