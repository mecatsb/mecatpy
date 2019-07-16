import re
import helper


class Compound:
    def __init__(self, entry):
        self.molweight = float('nan')
        self.comment = " "
        self.formula = " "
        self.reactions = ()
        self.atomtypes = dict()
        self.heterocycles = set()
        self.functional_groups = set()
        self.pathways = list()
        self.number_of_rings = 0
        self.brite = ''
        self.remark = ()
        self.parse(entry)
        if not hasattr(self, 'id'):
            raise ValueError
        if not hasattr(self, 'names'):
            self.names = [self.id]

    def parse(self, entry):
        fieldkeys = ['ENTRY', 'NAME', 'FORMULA', 'EXACT_MASS', 'MOL_WEIGHT', 'SEQUENCE', 'REMARK', 'COMMENT',
                     'REACTION', 'PATHWAY', 'BRITE', 'ENZYME', 'PATHWAY', 'OTHER DBS', 'LINKDB', 'MODULE',
                     'STRUCTURE', 'ATOM', 'BOND', 'DBLINKS']
        pattern = '|'.join(fieldkeys)

        data = [x.strip() for x in re.split(pattern, entry)[1:]]
        keys = re.findall(pattern, entry)
        for key, val in zip(keys, data):
            if key == 'ENTRY':
                self.id = re.search('C\\d{5}', val).group(0)
            elif key == 'NAME':
                self.names = [x.strip('; ') for x in val.split('\n')]
            elif key == 'FORMULA':
                self.formula = val.strip()
            elif key == 'COMMENT':
                self.comment = val.strip()
            elif key == 'REACTION':
                self.reactions = re.findall('R\d{5}', val)
            elif key == 'MOL_WEIGHT':
                self.molweight = helper.f(val.strip())
            elif key == 'PATHWAY':
                self.pathways = re.findall('map(\d{5})', val)
            elif key == 'ATOM':
                self.__parse_atomtypes(val.strip())
            elif key == 'BRITE':
                self.brite = re.findall('(br\d{5}|ko\d{5})', val)
            elif key == 'REMARK':
                self.remark = re.findall('[DG]\\d{5}', val)

    def id_is_valid(self):
        return re.match('C\d{5}', self.id) is not None

    def is_generic(self):
        if not hasattr(self, 'comment'):
            return False

        if 'generic' in self.comment or 'Generic' in self.comment:
            return True
        else:
            return hasattr(self, 'formula') and 'R' in self.formula

    def has_reaction(self):
        return hasattr(self, 'reactions') and len(self.reactions) > 0

    def __parse_atomtypes(self,at):
        heterocycles = set(["C7x", "N1x", "N1y", "N2x", "N2y", "N4x", "N4y", "N5x", "N5y", "O2x", "O7x", "S2x", "S3x"])
        functional_groups = set(["6a", "C8x", "C8y", "N4x", "N4y", "X", "N1a", "N1b", "N1c", "N1d", "O2a", "O2b", "O2c",
                                 "O2x", "C7a", "C7x"])
        atom_types = [x[1] for x in [d.strip().split() for d in at.split('\n')[1:]]]
        self.number_of_rings = sum('x' in s for s in atom_types) + sum('y' in s for s in atom_types)
        self.atomtypes = dict((x, atom_types.count(x)) for x in set(atom_types))
        self.heterocycles = set(atom_types) & heterocycles
        self.functional_groups = set(atom_types) & functional_groups

