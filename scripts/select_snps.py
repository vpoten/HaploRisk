import os
import subprocess

"""
Select snps using plink
"""

PLINK = 'plink'
CHRS = [str(i) for i in range(1, 23)]

INPUT_DIR = '/home/victor/Escritorio/Genotipado_Alternativo/data/WTCCC2/merged'
DATASET = 'wtccc2_merged{chr}'

OUTPUT_DIR = '/home/victor/Escritorio/Genotipado_Alternativo/data/WTCCC2/selection'
OUT_BASE = 'wtccc2_select_sp140_{chr}'

SNPS_FILE = '/home/victor/Escritorio/Genotipado_Alternativo/data/epistasis_sp140/snps_epistasis.txt'

if __name__ == "__main__":
    for chr in CHRS:
        dataset = os.path.join(INPUT_DIR, DATASET.replace('{chr}', chr))
        out_file = os.path.join(OUTPUT_DIR, OUT_BASE.replace('{chr}', chr))

        res = subprocess.call(
            [PLINK, '--noweb', '--bfile', dataset, '--extract', SNPS_FILE, '--make-bed', '--out', out_file])

        if res != 0:
            print 'Error processing chr %s' % chr