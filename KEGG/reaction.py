import re
import math
import helper
import copy
import itertools


class Reaction:
    #counter = 0 #debug Zwecke
    def __init__(self, entry):

        self.rclass = ()
        self.enzymes = ()
        self.comment = " "
        self.pathways = ()
        self.parse(entry)
        if not hasattr(self, 'id') or not hasattr(self, 'equation'):
            raise ValueError
        if not hasattr(self, 'names'):
            self.names = [self.id]

    def parse(self, entry):
        fieldkeys = ['ENTRY', 'NAME', 'DEFINITION', 'EQUATION', 'REMARK', 'COMMENT', 'RPAIR', 'RCLASS', 'ENZYME',
                     'PATHWAY', 'MODULE', 'ORTHOLOGY', 'LINKDB']
        pattern = '|'.join(fieldkeys)
        data = [x.strip() for x in re.split(pattern, entry)[1:]]
        keys = re.findall(pattern, entry)

        for key, val in zip(keys, data):
            if key == 'ENTRY':
                self.id = re.search('R\d{5}', val).group(0)
            elif key == 'NAME':
                self.names = [x.strip('; ') for x in val.split('\n')]
            elif key == 'RCLASS':
                self.rclass = re.findall('(C\d{5})_(C\d{5})', val)
            elif key == 'ENZYME':
                self.enzymes = re.findall('\d[.]\d+[.]\d+[.]\d+', val)
            elif key == 'EQUATION':
                self.equation = re.sub('///', '', val).strip()
            elif key == 'DEFINITION':
                self.definition = val.strip()
            elif key == 'COMMENT':
                self.comment = val.strip()
            elif key == 'PATHWAY':
                self.pathways = re.findall('rn(\d{5})', val)

        if not hasattr(self, 'equation'):
            raise ValueError

        self.__parse_equation(self.equation)

    def rpairs(self, mode):
        if mode == 'rclass':
            return self.rclass
        if mode == 'all':
            return list(itertools.product(self.__substrates, self.__products))
        if mode == 'both':
            if self.rclass:
                return self.rclass
            else:
                return list(itertools.product(self.__substrates, self.__products))

    def has_valid_stoichiometry(self):
        if any(math.isnan(v) or abs(v) < 0.5 for k, v in self.reactants.iteritems()):
            return False
        return True

    def has_valid_reactants(self):
        if any(re.match('G', k) for k, v in self.reactants.iteritems()):
            return False
        return True

    def is_g(self):
        return re.search('generic', self.comment) is not None

    def is_generic(self):
        return re.search('generic|incomplete|general', self.comment) is not None

    def is_incomplete(self):
        return re.search('incomplete', self.comment) is not None

    def is_general(self):
        return re.search('general', self.comment) is not None

    def is_non_enzymatic(self):
        return re.search('non-enzymatic', self.comment) is not None

    def is_product(self, compound_id):
        if compound_id in self.reactants:
            return self.reactants[compound_id] > 0
        else:
            return False

    def is_substrate(self, compound_id):
        if compound_id in self.reactants:
            return self.reactants[compound_id] < 0
        else:
            return False

    def is_substrate_product_distinct(self):
        return len(self.__reactants) == len(set(self.__reactants))

    def substrates(self):
        return [c for c, val in self.reactants.iteritems() if val < 0]

    def products(self):
        return [c for c, val in self.reactants.iteritems() if val > 0]

    def has_rclass(self):
        return (len(self.rclass)) > 0

    def has_enzyme(self):
        return (len(self.enzymes)) > 0

    def reverse(self):
        # TODO: equation und definition auch umdrehen
        reverse = copy.deepcopy(self)
        reverse.id = "-" + self.id
        for k in self.reactants:
            reverse.reactants[k] = -reverse.reactants[k]
        return reverse

    def __parse_equation(self, equation):
        tokens = ['-->', '<==>', '<=>']
        pattern = '|'.join(tokens)
        left, right = [x.strip() for x in re.split(pattern, equation)[0:]]

        self.reactants = dict()
        self.__reactants = list()
        self.__substrates = list()
        self.__products = list()
        for m in re.finditer(r'((?:[0-9]*|\([^)]*\)|[n,m]))\s*((?:C|G)[0-9]{5})', left):
            self.reactants[m.group(2)] = -1.0 if m.group(1) is '' else -helper.f(m.group(1))
            self.__reactants.append(m.group(2))
            self.__substrates.append(m.group(2))
        for m in re.finditer(r'((?:[0-9]*|\([^)]*\)|[n,m]))\s*((?:C|G)[0-9]{5})', right):
            stoich = 1.0 if m.group(1) is '' else helper.f(m.group(1))
            self.reactants[m.group(2)] = self.reactants[m.group(2)] + stoich if m.group(2) in self.reactants else stoich
            self.__reactants.append(m.group(2))
            self.__products.append(m.group(2))

    def __cmp__(self, other):
        return 1 if self.id>other.id else 0 if self.id == other.id else -1

    def __eq__(self, other):
        return other.id==self.id

    def has_changes(self, other):
        id = other.id != self.id
        equation = other.equation != self.equation
        definition = other.definition != self.definition
        comment = other.comment != self.comment
        rclass = other.rclass != self.rclass
        enzymes = other.enzymes != self.enzymes

        return id or rclass or equation or definition or comment or enzymes

    def has_changed_equation(self, other):
        return  other.equation != self.equation

    def has_changed_definition(self, other):
        return  other.definition != self.definition

    def has_changed_comment(self, other):
        return  other.comment != self.comment

    def has_changed_rclass(self, other):
        return  other.rclass != self.rclass

    def has_changed_enzymes(self, other):
        return  other.enzymes != self.enzymes