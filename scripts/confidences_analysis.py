import os
import datetime
from classes.snp_database import SnpDatabase
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt


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
        threshold = None

        if pos >= 0:
            values = sorted(values)
            zscores = stats.zscore(np.array(values))
            threshold = (values[pos], zscores[pos])

        snp_thresholds[probe_id] = threshold if threshold is not None else (-1, -1)

    return snp_thresholds


def write_threshold_result(chro, snp_thresholds, snp_db):
    f = open(os.path.join(crlmm_base, 'snp_thresholds_chr%s.txt' % chro), 'w')
    f.write('chro\tprobe_id\tsnp_id\tthreshold\tnorm_thr\tf_miss\n')
    for probe_id in snp_thresholds:
        rs_id = db.get_rs_id(chro, probe_id)
        lmiss_data = snp_db.get_snp_data(chro, rs_id).get('lmiss')
        threshold = snp_thresholds[probe_id]
        f.write(chro + '\t' + probe_id + '\t' + rs_id + '\t')
        f.write(str(threshold[0]) + '\t' + str(threshold[1]) + '\t' + str(lmiss_data['F_MISS']) + '\n')
    f.close()


def violin_plot(data):
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 5))

    # generate some random test data
    all_data = [np.array(values) for values in data]

    # plot violin plot
    axes[0].violinplot(all_data,
                       showmeans=False,
                       showmedians=True)
    axes[0].set_title('violin plot')

    # plot box plot
    axes[1].boxplot(all_data)
    axes[1].set_title('box plot')

    # adding horizontal grid lines
    for ax in axes:
        ax.yaxis.grid(True)
        ax.set_xticks([y + 1 for y in range(len(all_data))])
        ax.set_xlabel('xlabel')
        ax.set_ylabel('ylabel')

    # add x-tick labels
    plt.setp(axes, xticks=[y + 1 for y in range(len(all_data))],
             xticklabels=['x1', 'x2', 'x3', 'x4'])
    plt.show()


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
    write_threshold_result(chro, snp_thresholds, db)

    violin_plot(confidences['confidences'].get('SNP_A-8460085'))

    print 'Finished:', datetime.datetime.now().isoformat()
