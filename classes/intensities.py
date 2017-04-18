import os
import gzip
import numpy as np


class Intensities(object):
    @staticmethod
    def load_birdseed_summary_intensities(path, limit=-1):
        sep_char = '\t'
        f = gzip.open(path, 'r')

        intensities = {}  # dict with key=probe_id and value=np.array of intensities
        num_read = 0
        for line in f:
            if line.startswith('#'):
                continue  # skip comment
            elif line.startswith('probeset_id'):
                subjects = map(lambda s: os.path.splitext(s)[0], line.split(sep_char)[1:])
            else:
                toks = line.split(sep_char)
                intensities[toks[0]] = np.array(map(lambda i: float(i), toks[1:]))
                num_read += 1
            if num_read == limit:
                break

        f.close()
        return {'subjects': subjects, 'intensities': intensities}
