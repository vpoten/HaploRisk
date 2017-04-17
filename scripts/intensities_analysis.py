import os
import datetime
from classes.snp_database import SnpDatabase
from classes.tfam import Tfam
from classes.intensities import Intensities
import numpy as np
import matplotlib.pyplot as plt


if __name__ == "__main__":
    print 'Started:', datetime.datetime.now().isoformat()

    birdseed_base = '/home/victor/Escritorio/matesanz2015'
    birdseed2_base = '/home/victor/Escritorio/Genotipado_Alternativo/data/birdseed_out'

    tfam = Tfam(os.path.join(birdseed_base, 'matesanz2015.tfam'))
    db = SnpDatabase()
    db.load_from_birdseed(birdseed_base, '8090939')
    intensities = Intensities.load_birdseed_summary_intensities(
        os.path.join(birdseed2_base, 'birdseed-dev.summary.txt.gz'))

    parents_index = tfam.get_parents_index(intensities['subjects'])
    offspring_index = tfam.get_offspring_index(intensities['subjects'])

    print 'Finished:', datetime.datetime.now().isoformat()