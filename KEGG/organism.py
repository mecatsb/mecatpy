class Organism:
    def __init__(self, entry):
        self.comment = " "
        self.parse(entry)
        if not hasattr(self, 'id'):
            raise ValueError

    def parse(self, entry):
        splitted = entry.split('\t')
        self.id = splitted[0]
        self.name = splitted[1]
        self.definition = splitted[2]
        self.taxonomy = splitted[3]
