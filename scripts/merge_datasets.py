import os
import subprocess

PLINK = 'plink'
CHRS = [str(i) for i in range(1, 23)]

INPUT_DIR = '/home/victor/Escritorio/Genotipado_Alternativo/data/WTCCC2/bed'
DATASETS = ['1958BC_illumina{chr}', 'MSUK_illumina{chr}', 'NBS_illumina{chr}']

OUTPUT_DIR = '/home/victor/Escritorio/Genotipado_Alternativo/data/WTCCC2/merged'
OUT_BASE = 'wtccc2_merged{chr}'

if __name__ == "__main__":
    for chr in CHRS:
        datasets = map(lambda dat: os.path.join(INPUT_DIR, dat.replace('{chr}', chr)), DATASETS)
        out_file = os.path.join(OUTPUT_DIR, OUT_BASE.replace('{chr}', chr))
        list_file = os.path.join(OUTPUT_DIR, OUT_BASE.replace('{chr}', chr) + '_list.txt')

        f = open(list_file, 'w')
        for i in range(1, len(DATASETS)):
            f.write('%s.bed\t%s.bim\t%s.fam\n' % (datasets[i], datasets[i], datasets[i]))
        f.close()

        res = subprocess.call(
            [PLINK, '--noweb', '--bfile', datasets[0], '--merge-list', list_file, '--make-bed', '--out', out_file])

        if res != 0:
            print 'Error processing chr %s' % chr
