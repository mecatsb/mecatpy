import re


class Pathway:
    def __init__(self, entry):
        self.parse(entry)
        if not hasattr(self, 'id'):
            raise ValueError

    def parse(self, entry):
        splitted = entry.split('\t')
        self.id = re.search('\d{5}', splitted[0]).group(0)
        self.name = splitted[1]
