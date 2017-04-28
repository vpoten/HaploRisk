import os
import datetime
from classes.snp_database import SnpDatabase


if __name__ == "__main__":
    print 'Started:', datetime.datetime.now().isoformat()

    db_snp_hg38 = '/home/victor/Escritorio/Genotipado_Alternativo/data/dbSNP/snp147Common.txt.gz'

    db = SnpDatabase()
    db.load_from_ucsc_db_snp(db_snp_hg38)

    print 'Finished:', datetime.datetime.now().isoformat()
