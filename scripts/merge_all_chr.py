import os
import subprocess

"""
Merge all chr in one file using plink
"""

PLINK = 'plink'
CHRS = [str(i) for i in range(1, 23)]

INPUT_DIR = '/home/victor/Escritorio/Genotipado_Alternativo/data/WTCCC2/selection'
DATASET = 'wtccc2_select_sp140_{chr}'

OUTPUT_DIR = '/home/victor/Escritorio/Genotipado_Alternativo/data/WTCCC2/selection'
OUT_BASE = 'wtccc2_select_sp140_all'

if __name__ == "__main__":
    datasets = []

    for chr in CHRS:
        dataset = os.path.join(INPUT_DIR, DATASET.replace('{chr}', chr))
        if os.path.isfile(dataset+'.bed'):
            datasets.append(dataset)

    out_file = os.path.join(OUTPUT_DIR, OUT_BASE)
    list_file = out_file + '_list.txt'

    f = open(list_file, 'w')
    for i in range(1, len(datasets)):
        f.write('%s.bed\t%s.bim\t%s.fam\n' % (datasets[i], datasets[i], datasets[i]))
    f.close()

    res = subprocess.call(
        [PLINK, '--noweb', '--bfile', datasets[0], '--merge-list', list_file, '--recode', '--out', out_file])
