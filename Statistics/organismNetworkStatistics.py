from os.path import join
import constants as const
from pylab import *
import cPickle
import re
from scipy import stats


class OrganismNetworkStatistics:

    def __init__(self):

        with open(join(const.MECAT_BASE, "organisms.txt"), 'r') as f:
            organisms = [x.strip() for x in f.readlines()]

        hosts = """cge
        eco
        vna
        ppun
        mxa
        sce
        spo
        cgb
        mpe
        """.split()

        self.organisms = {k: organisms[k] for k in hosts}
        self.organisms = organisms

        with open(join(const.DATAPATH_KEGG, "KEGG_reaction.pickle"), 'rb') as handle:
            self.reactions = load(handle)

        with open(join(const.DATAPATH_KEGG, "KEGG_compound.pickle"), 'rb') as handle:
            self.compounds = load(handle)

        with open(join(const.DATAPATH_KEGG, "KEGG_enzyme.pickle"), 'rb') as handle:
            self.enzymes = load(handle)

        self.organism_reactions = []
        self.organism_reactions_model = []
        self.reactions_in_organism = dict()
        self.__reactions_per_organism()

    def __reaction_is_in_host(self, reaction, host):
        for e in reaction.enzymes:
            enzyme = self.enzymes['EC ' + e]
            if enzyme.is_in_organism(host):
                return True
        return False

    def __reaction_index(self, rid):
        if 'rIndex' not in self.__dict__:
            self.rIndex = {}
            for i, r in enumerate(self.reactions):
                self.rIndex[r] = i
        return self.rIndex[rid]

    def __reactions_per_organism(self):
        reactions = self.reactions.values()

        for o in self.organisms:
            self.reactions_in_organism[o] = [r for r in reactions if self.__reaction_is_in_host(r, o.upper())]
            self.organism_reactions.append((len(self.reactions_in_organism[o]), o))

        self.organism_reactions = sorted(self.organism_reactions, reverse=True)

    def reactions_in_model_per_organism(self):
        reactions_per_organism_kegg = dict((org, c) for c, org in self.organism_reactions)
        reactions_per_organism_model = dict()

        for o in self.organisms:
            if reactions_per_organism_kegg[o] > 0:

                reacs = [r for r in self.reactions_in_organism[o] if r.has_valid_stoichiometry()
                             and r.has_valid_reactants() and not r.is_generic() and r.has_rclass()
                             and not self.reaction_has_generic_compound(r.id)]
                reactions_per_organism_model[o] = len(reacs)

            else:
                reactions_per_organism_model[o] = 0

            reactions = filter(None, [r if self.reaction_is_in_host(r.id, o.upper())
                                      else None for k, r in self.reactions.iteritems()])

            has_valid_reactants = [r for r in reactions if r.has_valid_reactants()]
            valid_stoichiometry = [r for r in reactions if r.has_valid_stoichiometry()]
            generic_reactions = [r for r in reactions if r.is_generic()]
            reaction_has_generic_compound = [r for r in reactions if self.reaction_has_generic_compound(r.id)]

            reaction_ids = set([r.id for r in reactions])
            valid_reactant_ids = set([r.id for r in has_valid_reactants])
            valid_stoichiometry_ids = set([r.id for r in valid_stoichiometry])
            generic_reaction_ids = set([r.id for r in generic_reactions]) | \
                                   set([r.id for r in reaction_has_generic_compound])

            not_generic = reaction_ids - generic_reaction_ids
            valid_reactants = not_generic & valid_reactant_ids
            valid_stoichiometry = valid_reactants & valid_stoichiometry_ids

            assert sorted(valid_reactants) == sorted(not_generic)

            self.organism_reactions_model.append((reactions_per_organism_kegg[o], reactions_per_organism_model[o], len(valid_reactants), len(valid_stoichiometry), len(not_generic), o))
            self.organism_reactions_model = sorted(self.organism_reactions_model, key=lambda x: x[1], reverse=True) #key=lambda x: x[4]

        self.write_reactions_per_organism_kegg_model('reactions_organism_kegg_model.txt')

        kegg = np.array([float(x) for x, _, _, _, _, _ in self.organism_reactions_model])
        model = np.array([float(x) for _, x, _, _, _, _ in self.organism_reactions_model])

        with open(join(const.DATAPATH_KEGG, "kegg.pickle"), 'wb') as handle:
            cPickle.dump(kegg, handle, protocol=cPickle.HIGHEST_PROTOCOL)

        with open(join(const.DATAPATH_KEGG, "model.pickle"), 'wb') as handle:
            cPickle.dump(model, handle, protocol=cPickle.HIGHEST_PROTOCOL)

        with open(join(const.DATAPATH_STATISTICS_KEGG,'reactions_per_organism_kegg.pickle'), 'wb') as handle:
            cPickle.dump(reactions_per_organism_kegg, handle, protocol=cPickle.HIGHEST_PROTOCOL)

        with open(join(const.DATAPATH_STATISTICS_KEGG,'reactions_per_organism_model.pickle'), 'wb') as handle:
            cPickle.dump(reactions_per_organism_model, handle, protocol=cPickle.HIGHEST_PROTOCOL)

    def write_reactions_per_organism(self, filename):
        with open(join(const.DATAPATH_STATISTICS_KEGG, filename + 'KEGG' + '.txt'), 'w') as f:
            for c, o in self.organism_reactions:
                f.write(o + '\t' + self.organisms[o].definition + '\t' + str(c) + '\n')

    def write_reactions_per_organism_kegg_model(self, filename):
        with open(join(const.DATAPATH_STATISTICS_KEGG, filename), 'w') as f:
            for entry in self.organism_reactions_model:
                kegg = entry[0]
                model = entry[1]
                o = entry[5]
                f.write(o + '\t' + self.organisms[o].definition + '\t' + str(kegg) + '\t' + str(model) + '\t'
                        + str(kegg-model) + '\n')

    def plot_statistics(self, kegg, model):
        figure()
        ax = plt.axes()

        ls = 'None'
        ax.plot(kegg, marker='.', linestyle=ls, label='KEGG')
        ax.plot(model, marker='.', linestyle=ls, label='model')

        ax.set_xlabel("organism")
        ax.set_ylabel("number of reactions")

        slope, intercept, r_value, p_value, std_err = stats.linregress(kegg, model)
        line = slope * kegg + intercept

        ax2 = plt.axes([0.55, 0.55, 0.3, 0.3])
        ax2.plot(kegg, model, 'g.', linestyle='None', label='organism')
        ax2.plot(kegg, line, 'k-', label='linear regression')
        ax2.text(0.95, 0.01, '$r^2 =$ ' + str(round(r_value, 2)), verticalalignment='bottom',
                 horizontalalignment='right', transform=ax2.transAxes)

        ax2.set_xlabel("reactions in KEGG")
        ax2.set_ylabel("reactions in the model")
        ax2.legend(loc='upper left')
        ax.legend(loc='lower left')

        savefig(join(const.DATAPATH_STATISTICS_KEGG, "reactions per organism.pdf"))

    def reaction_statistics(self, host, filename):

        host = host.upper()
        if host == 'KEGG':
            reactions = self.reactions.values()
        else:
            reactions = filter(None, [r if self.reaction_is_in_host(r.id, host) else None for k, r in self.reactions.iteritems()])

        has_valid_reactants = [r for r in reactions if r.has_valid_reactants()]
        valid_stoichiometry = [r for r in reactions if r.has_valid_stoichiometry()]
        generic_reactions = [r for r in reactions if r.is_generic()]
        reactions_with_rclasses = [r for r in reactions if r.has_rclass()]
        reaction_has_generic_compound = [r for r in reactions if self.reaction_has_generic_compound(r.id)]
        reaction_has_enzyme = [r for r in reactions if r.has_enzyme()]

        all_reaction_ids = set([r.id for r in reactions])
        valid_reactant_ids = set([r.id for r in has_valid_reactants])
        valid_stoichiometry_ids = set([r.id for r in valid_stoichiometry])
        generic_reaction_ids = set([r.id for r in generic_reactions]) | set([r.id for r in reaction_has_generic_compound])
        rclass_ids = set([r.id for r in reactions_with_rclasses])
        enzyme_ids = set([r.id for r in reaction_has_enzyme])

        reactions = len(all_reaction_ids)
        generic = len(generic_reactions)
        rclass = len(reactions_with_rclasses)
        stoich = len(valid_stoichiometry)
        valid_reactants = len(has_valid_reactants)
        generic_reactants = len(reaction_has_generic_compound)

        valid_stoichiometry = valid_reactant_ids & valid_stoichiometry_ids
        not_generic = valid_stoichiometry - generic_reaction_ids
        has_rclass = not_generic & rclass_ids
        has_enzyme = has_rclass & enzyme_ids

        reactions_with_organisms = set()

        for rid in has_enzyme:
            reaction = self.reactions[rid]
            for e in reaction.enzymes:
                enzyme = self.enzymes['EC ' + e]
                if enzyme.has_organism():
                    reactions_with_organisms.add(rid)


        numvalid = len(valid_reactant_ids)
        numstoich = len(valid_stoichiometry)
        numnotgeneric = len(not_generic)
        numrclass = len(has_rclass)
        numenzyme = len(has_enzyme)

        pvalid = round(float(numvalid) / reactions*100,0)
        pstoich = round(float(numstoich)/ numvalid*100,0)
        pnotgeneric = round(float(numnotgeneric) / numstoich*100,0)
        prclass = round(float(numrclass) / numnotgeneric*100,0)
        penzyme = round(float(numenzyme)/ numrclass*100,0)
        porganism = round(float(len(reactions_with_organisms))/ numenzyme*100,0)

        wvalid = 300
        wstoich = pvalid * wvalid / 100
        wnotgeneric = pstoich * wstoich / 100
        wrclass = pnotgeneric * wnotgeneric / 100
        wenzyme = prclass * wrclass / 100
        worganism = penzyme * wenzyme / 100

        with open(join(const.DATAPATH_KEGG, filename + '.txt'), 'w') as f:
            f.write('total & ' + str(reactions) + ' \\\\' + '\n')
            f.write('with valid reactants & ' + str(valid_reactants) + ' \\\\' + '\n')

            f.write('with rclasses & ' + str(rclass) + ' \\\\' + '\n')
            f.write('with valid stoichiometry & ' + str(stoich) + ' \\\\' + '\n')
            f.write('generic & ' + str(generic) + ' \\\\' + '\n')
            f.write('with generic reactants & ' + str(generic_reactants) + ' \\\\' + '\n')

            f.write('\n' + 'fuer Abbildung' + '\n')
            f.write('valid reac: ' + str(numvalid) + '\t' + str(pvalid) +  '\n')
            f.write('valid stoich: ' + str(numstoich) + '\t' + str(pstoich)  + '\t' + str(wstoich) + '\t' + '\n')
            f.write('not generic: ' + str(numnotgeneric) + '\t' + str(pnotgeneric) + '\t' +  str(wnotgeneric) + '\t' + '\n')
            f.write('rclass: ' + str(numrclass) + '\t' + str(prclass) + '\t' +  str(wrclass) + '\t' +  '\n')
            f.write('enzyme: ' + str(len(has_enzyme)) + '\t' + str(penzyme) + '\t' +  str(wenzyme) + '\t' + '\n')
            f.write('organism: ' + str(len(reactions_with_organisms)) + '\t' + str(porganism) + '\t' +  str(worganism) + '\t' + '\n')

    def reaction_has_generic_compound(self, reaction_id):
        rid = re.findall('R\d{5}', reaction_id)[0]
        reaction = self.reactions[rid]

        for c, v in reaction.reactants.iteritems():
            try:
                compound = self.compounds[c]
            except KeyError:
                return True

            if compound.is_generic():
                return True

        return False

    def reaction_is_in_host(self, reaction_id, host):
        reaction = self.reactions[reaction_id]

        for e in reaction.enzymes:
            enzyme = self.enzymes['EC ' + e]
            if enzyme.is_in_organism(host):
                return True
        return False


if __name__ == '__main__':

    filename = "OrganismNetworkStatistics.pickle"

    statistics = OrganismNetworkStatistics()
    statistics.reaction_statistics('kegg', 'current_reactions_KEGG')
    statistics.reactions_in_model_per_organism()

    with open(join(const.DATAPATH_KEGG, "kegg.pickle"), 'rb') as handle:
       kegg = load(handle)

    with open(join(const.DATAPATH_KEGG, "model.pickle"), 'rb') as handle:
        model = load(handle)

    statistics.plot_statistics(kegg, model)
    statistics.write_reactions_per_organism('reactions_per_organism')

