import os
import datetime
import re
import sys
from classes.snp_database import SnpDatabase
from classes.tfam import Tfam

# constants
PED_FILE_NAME_REGEX = '_(\d+)'
MISSING_CHAR = '?'  # char code for missing
WINDOW_SIZES = [1e6, 500000.0, 250000.0, 100000.0, 50000.0, 20000.0, 10000.0]
HG38_POS = 'hg38_pos'  # hg38 position field

hg38_CHR_SIZES = {
    '1': 248956422,
    '2': 242193529,
    '3': 198295559,
    '4': 190214555,
    '5': 181538259,
    '6': 170805979,
    '7': 159345973,
    'X': 156040895,
    '8': 145138636,
    '9': 138394717,
    '11': 135086622,
    '10': 133797422,
    '12': 133275309,
    '13': 114364328,
    '14': 107043718,
    '15': 101991189,
    '16': 90338345,
    '17': 83257441,
    '18': 80373285,
    '20': 64444167,
    '19': 58617616,
    'Y': 57227415,
    '22': 50818468,
    '21': 46709983
}


def load_gwas_snps(path):
    """
    # Load GWAS snps from map files
    :param path: 
    :return: a map with key=chr and values=list  (id, pos, missing_par, missing_child)
    """
    gwas_snps = {}

    p = re.compile(PED_FILE_NAME_REGEX + '\.map')
    map_files = filter(lambda name: p.match(name), os.listdir(path))

    for map_file in map_files:
        chro = p.match(map_file).groups()[0]
        gwas_snps[chro] = []

        f = open(os.path.join(path, map_file), 'r')
        for line in f:
            toks = line.split('\t')
            gwas_snps[chro].append({'id': toks[1], 'pos': int(toks[3]), 'missing_par': 0, 'missing_child': 0})
        f.close()

    return gwas_snps


def load_missing_info_from_ped(path, gwas_snps, tfam):
    """
    Count amount of missing in ped files
    :param path: 
    :param gwas_snps: 
    :param tfam: 
    :return: 
    """
    p = re.compile(PED_FILE_NAME_REGEX + '\.ped')
    ped_files = filter(lambda name: p.match(name), os.listdir(path))

    # count missing data
    for ped_file in ped_files:
        chro = p.match(ped_file).groups()[0]
        chro_snps = gwas_snps[chro]

        f = open(os.path.join(path, ped_file), 'r')
        for line in f:
            toks = line.split('\t')
            subject = tfam.get_subject_from_fam(toks)
            missing_count_field = 'missing_par' if tfam.is_parent(subject) else 'missing_child'

            for i, snp in enumerate(chro_snps):
                idx = 6 + i * 2
                if toks[idx] == MISSING_CHAR:
                    snp[missing_count_field] += 1

        f.close()


def check_missing_with_plink_results(snp_db, gwas_snps):
    """
    Makes sure that amount of missing read from ped is equal to missing amount in plink result files
    :param snp_db: 
    :param gwas_snps: 
    :return: 
    """
    for chro in gwas_snps:
        chro_snps = gwas_snps[chro]
        count = 0
        for snp in chro_snps:
            data = snp_db.get_snp_data(chro, snp['id'])
            assert data['lmiss']['N_MISS'] == (snp['missing_par'] + snp['missing_child'])
            count += 1
            # add extra info to snp_db
            data['lmiss']['missing_par'] = snp['missing_par']
            data['lmiss']['missing_child'] = snp['missing_child']
        print 'chr', chro, count, 'checked'


def human_format(num):
    magnitude = 0
    while abs(num) >= 1000:
        magnitude += 1
        num /= 1000.0
        # add more suffixes if you need them
    return '%i%s' % (num, ['', 'K', 'M', 'G', 'T', 'P'][magnitude])


def chr_missing_windows(chro, snp_db, window_size, pos_field):
    chr_size = hg38_CHR_SIZES[chro]
    num_windows = int(chr_size / window_size) + (1 if int(chr_size % window_size) > 0 else 0)
    miss_fields = ['missing_par', 'missing_child']
    bins = [{'missing_par': 0, 'missing_child': 0, 'count': 0} for i in range(0, num_windows)]

    for snp in snp_db.get_chr_data(chro):
        bin_idx = int(snp[pos_field] / window_size)
        for f in miss_fields:
            bins[bin_idx][f] += snp['lmiss'][f]
        bins[bin_idx]['count'] += 1.0

    for bin in bins:
        # calculate mean on each bin
        denom = bin['count']
        if denom > 0:
            for f in miss_fields:
                bin[f] /= denom
    return bins


def generate_missing_tracks(snp_db, pos_field):
    chrs = [str(i) for i in range(1, 23)]  # valid chromosomes
    tracks = {human_format(wsize): {chro: chr_missing_windows(chro, snp_db, wsize, pos_field) for chro in chrs} for
              wsize in WINDOW_SIZES}
    return tracks


if __name__ == "__main__":
    print 'Started:', datetime.datetime.now().isoformat()

    imsgc_dbgap_base = '/home/victor/Escritorio/IMSGC_dbgap'
    ped_files_path = os.path.join(imsgc_dbgap_base, 'plink')

    tfam = Tfam(os.path.join(imsgc_dbgap_base, 'IMSGC_dbgap.tfam'))

    # downloaded from: ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/snp147Common.txt.gz
    db_snp_hg38 = '/home/victor/Escritorio/Genotipado_Alternativo/data/dbSNP/snp147.txt.gz'

    db_imsgc_snps = SnpDatabase()
    db_imsgc_snps.load_from_map_files(ped_files_path, PED_FILE_NAME_REGEX)
    db_imsgc_snps.load_missing_info(os.path.join(imsgc_dbgap_base, 'missing_output'))
    db_imsgc_snps.print_stats()
    db_imsgc_snps.add_ucsc_db_snp_position(db_snp_hg38, HG38_POS)

    gwas_snps = load_gwas_snps(ped_files_path)
    load_missing_info_from_ped(ped_files_path, gwas_snps, tfam)
    check_missing_with_plink_results(db_imsgc_snps, gwas_snps)

    print 'Finished:', datetime.datetime.now().isoformat()
