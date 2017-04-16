import os
import datetime
import gzip
from classes.snp_database import SnpDatabase
from classes.tfam import Tfam
from classes.intensities import Intensities
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
    f = gzip.open(path, 'r')

    header = f.readline()
    subjects = map(lambda s: os.path.splitext(os.path.basename(s[1:-1]))[0], header.split(sep_char))

    confidences = {}  # dict with key=probe_id and value=list of confidences
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


def write_threshold_result(out_dir, chro, snp_thresholds, snp_db):
    f = open(os.path.join(out_dir, 'snp_thresholds_chr%s.txt' % chro), 'w')
    f.write('chro\tprobe_id\tsnp_id\tthreshold\tnorm_thr\tf_miss\n')
    for probe_id in snp_thresholds:
        rs_id = snp_db.get_rs_id(chro, probe_id)
        lmiss_data = snp_db.get_snp_data(chro, rs_id).get('lmiss')
        threshold = snp_thresholds[probe_id]
        f.write(chro + '\t' + probe_id + '\t' + rs_id + '\t')
        f.write(str(threshold[0]) + '\t' + str(threshold[1]) + '\t' + str(lmiss_data['F_MISS']) + '\n')
    f.close()


def violin_plot(data_missing, data_no_missing):
    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(12, 5))

    # plot snps with missing
    axes[0].violinplot(data_missing.values(),
                       showmeans=False,
                       showmedians=True)
    axes[0].set_title('SNPs with missing')

    # plot snps with no missing
    axes[1].violinplot(data_no_missing.values(),
                       showmeans=False,
                       showmedians=True)
    axes[1].set_title('SNPs no missing')

    # adding horizontal grid lines
    for ax in axes:
        ax.yaxis.grid(True)

    # add x-tick labels
    for idx, val in enumerate([data_missing, data_no_missing]):
        plt.setp(axes[idx], xticks=[y + 1 for y in range(len(val))], xticklabels=val.keys())
    plt.show()


def get_zscores(probe_id, confidences):
    return stats.zscore(np.array(confidences['confidences'].get(probe_id)))


def plot_fmiss_vs_threshold(chro, snp_thresholds, snp_db):
    x = []
    y = []
    y_norm = []

    for probe_id in snp_thresholds:
        rs_id = snp_db.get_rs_id(chro, probe_id)
        lmiss_data = snp_db.get_snp_data(chro, rs_id).get('lmiss')
        if lmiss_data['N_MISS'] == 0:
            continue
        thrs = snp_thresholds[probe_id]
        x.append(lmiss_data['F_MISS'])
        y.append(thrs[0])
        y_norm.append(thrs[1])

    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 5))

    axes[0].scatter(x, y, s=[1 for i in range(len(x))])
    axes[0].set_title('Missing rate vs threshold')

    axes[1].scatter(x, y_norm, s=[1 for i in range(len(x))])
    axes[1].set_title('Missing rate vs norm. threshold')

    for ax in axes:
        ax.set_xlabel('missing rate birdseed')
        ax.set_ylabel('threshold crlmm')

    plt.show()


def oligo_vs_birdseed_missing_threshold_cmp():
    birdseed_base = '/home/victor/Escritorio/matesanz2015'
    crlmm_base = '/home/victor/Escritorio/Genotipado_Alternativo/data'

    db = SnpDatabase()
    db.load_from_birdseed(birdseed_base, '8090939')
    db.load_missing_info(os.path.join(birdseed_base, 'missing_output'))

    target_missing_stats = load_missing_stats(os.path.join(birdseed_base, 'missing_output/missing_stats.txt'))

    confs_file = os.path.join(crlmm_base, 'crlmm_out/confs.txt.gz')
    confidences = load_oligo_confidences(confs_file)

    chro = '1'
    snp_thresholds = calc_missing_threshold(chro, db, confidences)
    write_threshold_result(crlmm_base, chro, snp_thresholds, db)

    snps_missing = ['SNP_A-8358491', 'SNP_A-2187372', 'SNP_A-4246654', 'SNP_A-2026459', 'SNP_A-8630776',
                    'SNP_A-8673051', 'SNP_A-1961190', 'SNP_A-8327571']
    snps_no_missing = ['SNP_A-4250282', 'SNP_A-4290848', 'SNP_A-8426563', 'SNP_A-1899510', 'SNP_A-2060328',
                       'SNP_A-1802035', 'SNP_A-4275453', 'SNP_A-8609627']

    # plot violin graphs of normalized confidences
    violin_plot({db.get_rs_id(chro, id): get_zscores(id, confidences) for id in snps_missing},
                {db.get_rs_id(chro, id): get_zscores(id, confidences) for id in snps_no_missing})

    # plot violin graph of non normalized confidences
    violin_plot({db.get_rs_id(chro, id): confidences['confidences'].get(id) for id in snps_missing},
                {db.get_rs_id(chro, id): confidences['confidences'].get(id) for id in snps_no_missing})

    plot_fmiss_vs_threshold(chro, snp_thresholds, db)


if __name__ == "__main__":
    print 'Started:', datetime.datetime.now().isoformat()

    birdseed_base = '/home/victor/Escritorio/matesanz2015'
    birdseed2_base = '/home/victor/Escritorio/Genotipado_Alternativo/data/birdseed_out'

    tfam = Tfam(os.path.join(birdseed_base, 'matesanz2015.tfam'))
    db = SnpDatabase()
    db.load_from_birdseed(birdseed_base, '8090939')
    intensities = Intensities.load_birdseed_summary_intensities(
        os.path.join(birdseed2_base, 'birdseed-dev.summary.txt.gz'))

    print 'Finished:', datetime.datetime.now().isoformat()
