import re


class Enzyme:
    def __init__(self, entry):
        self.comment = " "
        self.organisms = set()
        self.parse(entry)
        if not hasattr(self, 'id'):
            raise ValueError
        if not hasattr(self, 'names'):
            self.names = [self.id]

    def parse(self, entry):
        fieldkeys = ['ENTRY', 'NAME', 'CLASS', 'SYSNAME', 'SUBSTRATE', 'PRODUCT', 'COMMENT', 'GENES', 'HISTORY',
                     'PATHWAY', 'ORTHOLOGY', 'REFERENCE', 'DBLINKS', 'REACTION', 'ALL_REAC']
        pattern = '|'.join(fieldkeys)

        data = [x.strip() for x in re.split(pattern, entry)[1:]]
        keys = re.findall(pattern, entry)
        for key, val in zip(keys, data):
            if key == 'ENTRY':
                self.id = re.search('EC\s*\d[.]\d+[.]\d+[.]\d+', val).group(0)
            elif key == 'NAME':
                self.names = [x.strip('; ') for x in val.split('\n')]
            elif key == 'COMMENT':
                self.comment = val.strip()
            elif key == 'GENES':
                self.organisms = set(re.findall('(\w{3,4})[:]', val))

    def is_in_organism(self, organism):
        return organism in self.organisms

    def has_organism(self):
        return (len(self.organisms)) > 0
