import gzip
import re
import os


class SnpDatabase(object):
    F_PROBEID = 0
    F_BASES = 6  # 6 or 7
    F_RSID = 8
    F_POS = 9

    def __init__(self):
        self.snp_map = {}  # rs_id to position
        self.snp_map_pos = {}  # position to rs_id
        self.probe_id_map = {}  # probe_id to rs_id

    def load_from_birdseed(self, path, subject):
        p = re.compile(subject + '\.birdseed-v2.(\w+)\.txt\.gz')
        files = filter(lambda name: p.match(name), os.listdir(path))

        for i in range(1, 23):
            chro = str(i)
            self.snp_map[chro] = {}  # rs_id to position
            self.snp_map_pos[chro] = {}  # position to rs_id
            self.probe_id_map[chro] = {}  # probe_id to rs_id
            file_name = filter(lambda name: p.match(name).groups(0)[0] == chro, files)[0]
            self.__load_from_birdseed_out(chro, os.path.join(path, file_name))

    def __load_from_birdseed_out(self, chro, file_path):
        with (gzip.open(file_path, 'r') if file_path.endswith('.gz') else open(file_path, 'r')) as f:
            for line in f:
                if not line.startswith('#') and not line.startswith('Probe'):
                    toks = line.split('\t')
                    rs_id = toks[self.F_RSID]
                    probe_id = toks[self.F_PROBEID]
                    pos = int(toks[self.F_POS])

                    if rs_id not in self.snp_map[chro]:
                        self.snp_map[chro][rs_id] = pos
                        self.snp_map_pos[chro][pos] = rs_id
                        self.probe_id_map[chro][probe_id] = rs_id

            f.close()

    def get_rs_id(self, probe_id):
        return self.probe_id_map.get(probe_id)

    def get_position(self, rs_id):
        return self.snp_map.get(rs_id)
