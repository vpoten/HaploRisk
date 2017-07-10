import os
import sys

"""
convert a ped/map file to arff
"""

# ped values
PED_VALS = ['A', 'T', 'G', 'C', '0', '?', 'N']
MISSING = ['0', '?', 'N']

# arff attributes
ATT_VALUES = '0,1,2'
CLASS_VALUES = '1,2'
PHENO_ATT_NAME = 'phenotype'

# ped file fields
F_PHENO = 5
F_FAM = 0
F_SUBJ = 1


def gen_arff(out_file, chr_snps, map_geno, phenotypes):
    writer = open(out_file, 'w')

    writer.write("@RELATION %s" % os.path.basename(out_file))
    writer.write('\n')

    for chr in chr_snps:
        # write attributes
        snp_map = chr_snps[chr]
        for id in snp_map:
            writer.write("@ATTRIBUTE '%s:%s' {%s}" % (chr, id, ATT_VALUES))

    writer.write("@ATTRIBUTE %s {%s}\n" % (PHENO_ATT_NAME, CLASS_VALUES))
    writer.write('\n')
    writer.write('@DATA\n')

    for i, pheno in enumerate(phenotypes):
        writer.write("%%%s\n" % pheno['id'])

        for chr in map_geno:
            genotypes = map_geno[chr]
            assert len(chr_snps[chr]) == len(genotypes[i])
            for it in genotypes[i]:
                writer.write("%s," % '?' if it is None else it)

        writer.write(pheno['value'])
        writer.write('\n')

    writer.close()


def read_map(map_file, chr_snps):
    f = open(map_file, 'r')
    snp_list = []

    for line in file:
        # chromosome rs_id gen_dist bp
        toks = line.split('\t')
        snp_map = chr_snps.get(toks[0], {})
        chr_snps[toks[0]] = snp_map
        data = (toks[1], toks[0], int(toks[3]))  # rs_id, chr, position
        snp_map[toks[1]] = data
        snp_list.append(data)

    f.close()
    return snp_list


if __name__ == "__main__":
    ped_file = sys.argv[0]
    # TODO
