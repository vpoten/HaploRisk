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


def calc_missing_threshold(chro, snp_db, oligo_confs):
    probe_ids = snp_db.get_probe_ids(chro)
    confs = oligo_confs['confidences']
    snp_thresholds = {}

    for probe_id in probe_ids:
        rs_id = snp_db.get_rs_id(chro, probe_id)
        values = confs.get(probe_id)
        lmiss_data = None

        if rs_id is not None:
            lmiss_data = snp_db.get_snp_data(chro, rs_id).get('lmiss')

        if values is None or rs_id is None or lmiss_data is None:
            continue

        pos = lmiss_data['N_MISS'] - 1
        threshold = sorted(values)[pos] if pos >= 0 else None
        snp_thresholds[probe_id] = threshold if threshold is not None else -1

    return snp_thresholds


if __name__ == "__main__":
    birdseed_base = '/home/victor/Escritorio/matesanz2015'
    crlmm_base = '/home/victor/Escritorio/Genotipado_Alternativo/data'

    print 'Started:', datetime.datetime.now().isoformat()

    db = SnpDatabase()
    db.load_from_birdseed(birdseed_base, '8090939')
    db.load_missing_info(os.path.join(birdseed_base, 'missing_output'))

    target_missing_stats = load_missing_stats(os.path.join(birdseed_base, 'missing_output/missing_stats.txt'))

    confs_file = os.path.join(crlmm_base, 'crlmm_out/confs.txt')
    confidences = load_oligo_confidences(confs_file)

    chro = '1'
    snp_thresholds = calc_missing_threshold(chro, db, confidences)

    f = open(os.path.join(crlmm_base, 'snp_thresholds.txt'), 'w')
    f.write('chro\tprobe_id\tsnp_id\tthreshold\n')
    for probe_id in snp_thresholds:
        rs_id = db.get_rs_id(chro, probe_id)
        f.write(chro + '\t' + probe_id + '\t' + rs_id + '\t' + str(snp_thresholds[probe_id]) + '\n')
    f.close()

    print 'Finished:', datetime.datetime.now().isoformat()
