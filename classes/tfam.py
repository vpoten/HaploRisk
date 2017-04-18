import collections

Subject = collections.namedtuple('Subject', ['fid', 'id', 'father', 'mother', 'sex', 'phenotype'])


class Tfam(object):
    def __init__(self, tfam_file):
        """
        TFAM fields:
        1. Family ID ('FID')
        2. Within-family ID ('IID'; cannot be '0')
        3. Within-family ID of father ('0' if father isn't in dataset)
        4. Within-family ID of mother ('0' if mother isn't in dataset)
        5. Sex code ('1' = male, '2' = female, '0' = unknown)
        6. Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control) 
        """
        self.parents = []
        self.offspring = []
        self.subjects = {}

        f = open(tfam_file, 'r')

        for line in f:
            toks = line.rstrip('\n').split('\t')
            iid = toks[1]
            subject = Subject(fid=toks[0], id=iid, father=toks[2], mother=toks[3], sex=toks[4],
                              phenotype=toks[5])
            self.subjects[iid] = subject

            if subject.father == '0' and subject.mother == '0':
                self.parents.append(subject.id)
            else:
                self.offspring.append(subject.id)

        f.close()

    def get_subject(self, id):
        return self.subjects.get(id)

    def get_parents(self):
        return self.parents

    def get_offspring(self):
        return self.offspring

    def get_parents_index(self, subjects):
        return map(lambda pid: subjects.index(pid), self.parents)

    def get_offspring_index(self, subjects):
        return map(lambda oid: subjects.index(oid), self.offspring)
