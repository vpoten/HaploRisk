import os
import datetime
import re
from classes.snp_database import SnpDatabase
from classes.tfam import Tfam

# constants
ped_file_name_regex = '_(\d+)'
missing_char = '?'  # char code for missing


def load_gwas_snps(path):
    """
    # Load GWAS snps from map files
    :param path: 
    :return: a map with key=chr and values=list  (id, pos, missing_par, missing_child)
    """
    gwas_snps = {}

    p = re.compile(ped_file_name_regex + '\.map')
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
    p = re.compile(ped_file_name_regex + '\.ped')
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
                if toks[idx] == missing_char:
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
        print 'chr', chro, count, 'checked'


if __name__ == "__main__":
    print 'Started:', datetime.datetime.now().isoformat()

    imsgc_dbgap_base = '/home/victor/Escritorio/IMSGC_dbgap'
    ped_files_path = os.path.join(imsgc_dbgap_base, 'plink')

    tfam = Tfam(os.path.join(imsgc_dbgap_base, 'IMSGC_dbgap.tfam'))

    # downloaded from: ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/snp147Common.txt.gz
    db_snp_hg38 = '/home/victor/Escritorio/Genotipado_Alternativo/data/dbSNP/snp147Common.txt.gz'
    # db_hg38 = SnpDatabase()
    # db_hg38.load_from_ucsc_db_snp(db_snp_hg38)

    db_imsgc_snps = SnpDatabase()
    db_imsgc_snps.load_from_map_files(ped_files_path, ped_file_name_regex)
    db_imsgc_snps.load_missing_info(os.path.join(imsgc_dbgap_base, 'missing_output'))
    db_imsgc_snps.print_stats()

    gwas_snps = load_gwas_snps(ped_files_path)
    load_missing_info_from_ped(ped_files_path, gwas_snps, tfam)
    check_missing_with_plink_results(db_imsgc_snps, gwas_snps)

    print 'Finished:', datetime.datetime.now().isoformat()
