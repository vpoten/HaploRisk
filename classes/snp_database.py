import gzip
import re
import os


class SnpDatabase(object):
    F_PROBEID = 0
    F_BASES = 6  # 6 or 7
    F_RSID = 8
    F_POS = 9

    # UCSC dbSNP file fields
    F_UCSC_RSID = 4
    F_UCSC_POS = 2
    F_UCSC_CHR = 1

    def __init__(self):
        self.snp_map = {}  # rs_id to position
        self.snp_map_pos = {}  # position to rs_id
        self.probe_id_map = {}  # probe_id to rs_id
        self.probe_id_map_direct = {}  # probe_id to rs_id

    def load_from_ucsc_db_snp(self, file_path):
        valid_chrs = [str(i) for i in range(1, 23)] + ['X', 'Y']

        with (gzip.open(file_path, 'r') if file_path.endswith('.gz') else open(file_path, 'r')) as f:
            for line in f:
                toks = line.split('\t')
                rs_id = toks[self.F_UCSC_RSID]
                pos = int(toks[self.F_UCSC_POS])
                chro = toks[self.F_UCSC_CHR][3:]

                if chro not in valid_chrs:
                    continue

                if chro not in self.snp_map:
                    self.snp_map[chro] = {}
                    self.snp_map_pos[chro] = {}

                if rs_id not in self.snp_map[chro]:
                    self.snp_map[chro][rs_id] = {'id': rs_id, 'position': pos}
                    self.snp_map_pos[chro][pos] = rs_id

            f.close()

    def add_ucsc_db_snp_position(self, file_path, new_field_name):
        valid_chrs = self.snp_map.keys()

        with (gzip.open(file_path, 'r') if file_path.endswith('.gz') else open(file_path, 'r')) as f:
            for line in f:
                toks = line.split('\t')
                rs_id = toks[self.F_UCSC_RSID]
                pos = int(toks[self.F_UCSC_POS])
                chro = toks[self.F_UCSC_CHR][3:]

                if chro not in valid_chrs:
                    continue

                data = self.get_snp_data(chro, rs_id)

                if data is not None:
                    data[new_field_name] = pos

            f.close()

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
                        self.snp_map[chro][rs_id] = {'id': rs_id, 'position': pos}
                        self.snp_map_pos[chro][pos] = rs_id
                        self.probe_id_map[chro][probe_id] = rs_id
                        self.probe_id_map_direct[probe_id] = rs_id

            f.close()

    def load_missing_info(self, path, lmiss_regex='_(\d+)\.lmiss'):
        p = re.compile(lmiss_regex)
        files = filter(lambda name: p.match(name), os.listdir(path))

        for i in range(1, 23):
            chro = str(i)
            file_name = filter(lambda name: p.match(name).groups(0)[0] == chro, files)[0]
            self.__load_lmiss_out(chro, os.path.join(path, file_name))

    def __load_lmiss_out(self, chro, file_path):
        with (gzip.open(file_path, 'r') if file_path.endswith('.gz') else open(file_path, 'r')) as f:
            f.readline()  # skip header
            for line in f:
                toks = filter(None, line.split(' '))
                # fields: CHR SNP N_MISS N_GENO F_MISS
                rs_id = toks[1]
                data = self.get_snp_data(chro, rs_id)
                if data is not None:
                    data['lmiss'] = {'N_MISS': int(toks[2]), 'N_GENO': int(toks[3]), 'F_MISS': float(toks[4])}
                else:
                    print 'lmiss: %(snp)s in chr %(chr)s not found' % {'snp': rs_id, 'chr': chro}

            f.close()

    def load_from_map_files(self, path, ped_file_name_regex):
        p = re.compile(ped_file_name_regex + '\.map')
        map_files = filter(lambda name: p.match(name), os.listdir(path))

        for map_file in map_files:
            chro = p.match(map_file).groups()[0]
            self.snp_map[chro] = {}
            self.snp_map_pos[chro] = {}

            f = open(os.path.join(path, map_file), 'r')
            for line in f:
                toks = line.split('\t')
                rs_id = toks[1]
                pos = int(toks[3])
                if rs_id not in self.snp_map[chro]:
                    self.snp_map[chro][rs_id] = {'id': rs_id, 'position': pos}
                    self.snp_map_pos[chro][pos] = rs_id
            f.close()

    def get_rs_id(self, chro, probe_id):
        if chro is None:
            return self.probe_id_map_direct.get(probe_id)
        return self.probe_id_map[chro].get(probe_id)

    def get_snp_data(self, chro, rs_id):
        return self.snp_map[chro].get(rs_id)

    def get_chr_data(self, chro):
        return self.snp_map[chro].values()

    def get_position(self, chro, rs_id):
        return self.get_snp_data(chro, rs_id)['position']

    def get_snp_ids(self, chro):
        return self.snp_map[chro].keys()

    def get_probe_ids(self, chro):
        return self.probe_id_map[chro].keys()

    def print_stats(self):
        for chro in self.snp_map:
            print "Chr", chro, len(self.snp_map[chro]), 'snps'

    def get_snps_in_region(self, chro, min_pos, max_pos, pos_field='position'):
        chr_snps = self.snp_map[chro]
        return filter(lambda rs_id: max_pos > chr_snps[rs_id][pos_field] >= min_pos, chr_snps)

    def write_position(self, file_path, pos_field):
        f = open(file_path, 'w')
        for chro in self.snp_map:
            for snp in self.snp_map[chro].values():
                if snp.get(pos_field) is not None:
                    f.write(snp['id'] + '\t' + chro + '\t' + str(snp[pos_field]) + '\n')
        f.close()

    def read_position(self, file_path, pos_field):
        f = open(file_path, 'r')
        for line in f:
            toks = line.split('\t')
            rs_id = toks[0]
            chro = toks[1]
            pos = int(toks[2])
            data = self.get_snp_data(chro, rs_id)
            if data is not None:
                data[pos_field] = pos
        f.close()
