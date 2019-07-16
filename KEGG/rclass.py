import re


class RClass:
    def __init__(self, entry):
        self.parse(entry)
        if not hasattr(self, 'id'):
            raise ValueError

    def parse(self, entry):
        fieldkeys = ['ENTRY', 'RPAIR', 'REACTION', 'ENZYME', 'PATHWAY']
        pattern = '|'.join(fieldkeys)

        data = [x.strip() for x in re.split(pattern, entry)[1:]]
        keys = re.findall(pattern, entry)
        for key, val in zip(keys, data):
            if key == 'ENTRY':
                self.id = re.search('RC\d{5}', val).group(0)
            elif key == 'RPAIR':
                self.rpair = re.findall('(C\d{5})_(C\d{5})', val)
            elif key == 'ENZYME':
                self.enzyme = re.findall('(\d[.]\d+[.]\d+[.]\d+)', val)
            elif key == 'REACTION':
                self.reaction = re.findall('R\d{5}', val)
            elif key == 'PATHWAY':
                self.pathway = re.findall('(rn\d{5})\s*(.*)', val)
