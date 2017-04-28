import os
import datetime
import re
from classes.snp_database import SnpDatabase
from classes.tfam import Tfam

# constants
ped_file_name_regex = '_(\d+)'


def load_gwas_snps(path):
    """
    # Load GWAS snps from map files
    :param path: 
    :return: a map with key=chr and values=list snps
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
            gwas_snps[chro].append({'id': toks[1], 'pos': int(toks[3])})
        f.close()

    return gwas_snps


def load_missing_info_from_ped(path, gwas_snps, tfam):
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
            for i, snp in enumerate(chro_snps):
                # TODO <-----
                pass

        f.close()
    # TODO <-----

if __name__ == "__main__":
    print 'Started:', datetime.datetime.now().isoformat()

    imsgc_dbgap_base = '/home/victor/Escritorio/IMSGC_dbgap'

    tfam = Tfam(os.path.join(imsgc_dbgap_base, 'IMSGC_dbgap.tfam'))

    # downloaded from: ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/snp147Common.txt.gz
    db_snp_hg38 = '/home/victor/Escritorio/Genotipado_Alternativo/data/dbSNP/snp147Common.txt.gz'

    db_hg38 = SnpDatabase()
    db_hg38.load_from_ucsc_db_snp(db_snp_hg38)

    print 'Finished:', datetime.datetime.now().isoformat()
