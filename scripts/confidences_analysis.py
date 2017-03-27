import os
import datetime
from classes.snp_database import SnpDatabase


def load_missing_stats(path):
    stats = {}
    f = open(path, 'r')
    f.readline()  # skip header
    for line in f:
        toks = line.split('\t')
        chro = toks[1]
        stats[chro] = {'nparents': int(toks[2]), 'nchildren': int(toks[3]), 'avgMissPar': float(toks[4]),
                       'avgMissChild': float(toks[5])}
    f.close()
    return stats


def load_oligo_confidences(path):
    sep_char = ' '
    f = open(path, 'r')

    header = f.readline()
    subjects = map(lambda s: os.path.splitext(os.path.basename(s[1:-1]))[0], header.split(sep_char))

    confidences = {}
    for line in f:
        toks = line.split(sep_char)
        confidences[toks[0][1:-1]] = map(lambda i: float(toks[i]), range(1, len(toks)))

    f.close()
    return {'subjects': subjects, 'confidences': confidences}


if __name__ == "__main__":
    birdseed_base = '/home/victor/Escritorio/matesanz2015'
    crlmm_base = '/home/victor/Escritorio/Genotipado_Alternativo/data'

    print 'Started:', datetime.datetime.now().isoformat()

    db = SnpDatabase()
    db.load_from_birdseed(birdseed_base, '8090939')

    target_missing_stats = load_missing_stats(os.path.join(birdseed_base, 'missing_output/missing_stats.txt'))

    confs_file = os.path.join(crlmm_base, 'crlmm_out/confs.txt')
    confidences = load_oligo_confidences(confs_file)

    print 'Finished:', datetime.datetime.now().isoformat()
