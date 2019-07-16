import re
from os.path import join
import cPickle
import constants as const


class Thermodynamics:

    def __init__(self):

        self.data = dict()

    def read_data(self, filename):
        with open(filename, 'r') as f:
            entries = f.readlines()

        entries = [re.split(r'\t+', l.rstrip('\t')) for l in entries]
        entries = [e for e in entries if e]

        for entry in entries[1:]:
            a, b, c, d, e = entry
            self.data[a] = float(c), float(d), float(e)

        self.__set_reversibilities(-15)

    def dump_data(self):
        with open(join(const.DATAPATH_THERMODYNAMICS, 'equilibrator_data.pickle'), 'wb') as handle:
            cPickle.dump(self.data, handle, protocol=cPickle.HIGHEST_PROTOCOL)

    def read_dump(self):
        with open(join(const.DATAPATH_THERMODYNAMICS, 'equilibrator_data.pickle'), 'rb') as handle:
            self.data = cPickle.load(handle)

    def __set_reversibilities(self, threshold):

        for reaction_id, tdata in self.data.iteritems():

            dg0, dgu, dgm = tdata

            if abs(dgm) > 30 and dgu > 5 * abs(dgm):
                direction = 1

            elif -abs(threshold) < dgm < abs(threshold):
                direction = 0

            elif dgm < threshold:
                direction = 1

            elif dgm > abs(threshold):
                direction = -1

            else:
                direction = 1

            self.data[reaction_id] = dg0, dgu, dgm, direction

            if direction < 0:
                self.__modify_data_entry(reaction_id)

    def __modify_data_entry(self, keggid):

        entry = self.data[keggid]
        dg0, dgu, dgm, direction = entry
        self.data[keggid] = -dg0, dgu, -dgm, direction


if __name__ == '__main__':

    thermodynamics_data_file = "dGs_050619.txt"
    thermodynamics = Thermodynamics()
    thermodynamics.read_data(join(const.DATAPATH_THERMODYNAMICS, thermodynamics_data_file))
    thermodynamics.dump_data()

