import cPickle
import re
import constants as const
from os.path import join
from scipy import io as sio
import graph_tool.all as gt
from collections import defaultdict
from KEGG.compound import Compound as kCompound
import numpy as np
import traceback
import sys
from scipy import sparse


class Model:
    """
    Model class that represents a metabolic network
    """

    def __init__(self, reactions, metabolites, stoich_matrix, name):
        self.reactions = reactions
        self.metabolites = metabolites
        self.stoich_matrix = stoich_matrix
        self.name = name


class Modelbuilder:

    def __init__(self, name, path):

        with open(join(const.DATAPATH_KEGG, "KEGG_compound.pickle"), 'rb') as handle:
            self.compounds = cPickle.load(handle)

        with open(join(const.DATAPATH_KEGG, "KEGG_reaction.pickle"), 'rb') as handle:
            self.reactions = cPickle.load(handle)

        with open(join(const.DATAPATH_KEGG, "KEGG_enzyme.pickle"), 'rb') as handle:
            self.enzymes = cPickle.load(handle)

        with open(join(const.DATAPATH_THERMODYNAMICS, "equilibrator_data.pickle"), 'rb') as handle:
            self.thermodynamics = cPickle.load(handle)

        self.reactiongraph = gt.Graph()
        self.arcgraph = gt.Graph()  # rpairs, no cofactors/inorganics
        self.arcgraph_rclass_all = gt.Graph()  # rpairs
        self.arcgraph_rpairs_all = gt.Graph()  # all substrate-product pairs
        self.arcgraph_rpairs_filter = gt.Graph()  # all substrate-product pairs, no cofactors/inorganics

        self.host = ''  # graph property
        self.filtered_compounds_weight = []
        self.filtered_compounds_carbons = []
        self.filtered_compounds_no_carbons = []
        self.compounds_no_arcs = []
        self.basis_compound_list = []  # vertex property
        self.cofactor_list = []  # vertex property
        self.start_compound_list = []  # vertex property
        self.long_start_list = []
        self.generic_compound_list = []  # vertex property
        self.reversible_reactions = []
        self.rpairs = set()
        self.reaction_rpair = dict()
        self.reaction_reversibilities = dict()
        self.rpair_reaction = defaultdict(list)  # edge property
        self.feasible_reaction_list = set()
        self.all_target_compounds = set()
        self.targets_that_are_startmetabolites = set()
        self.path = path
        self.name = name
        self.all_reactions = []

    def buildKEGG(self, host, mode):
        self.host = host.upper()
        with_thermodynamics = True
        reacs = self.reactions.values()
        filter_reactions = True
        if mode == 'extended':
            filter_reactions = False
        if mode == 'reversible':
            with_thermodynamics = False

        self.__populate(reacs, filter_reactions, with_thermodynamics)

    def build_host(self, host, mode):
        self.host = host.upper()
        reacs = filter(None, [r if self.reaction_is_in_host(r.id) else None for k, r in self.reactions.iteritems()])

        filter_reactions = True
        if mode == 'extended':
            filter_reactions = False
        self.__populate(reacs, filter_reactions, self.host)

    def build(self, reactionlist, with_thermodynamics=False):
        reacs = [self.reactions[r[1:]].reverse() if r.startswith('-') else self.reactions[r] for r in reactionlist]
        self.__populate(reacs, filterreacs=False, with_thermodynamics=with_thermodynamics)

    def build_from_arcs(self, arcsfile, reactionfile):
        with open(arcsfile, 'r') as f:
            reacids = set([re.search('\d+', line).group(0).strip() for line in f if line.strip() != ''])

        with open(reactionfile, 'r') as f:
            rkeggids = [re.search('R\d{5}', line).group(0).strip() for line in f if line.strip() != '']

        reactions = [rkeggids[int(rid)] + '\n' for rid in reacids]
        reacs = set([self.reactions[re.search('R\d{5}', r).group(0).strip()] for r in reactions])

        self.__populate(reacs, True, 'complete')

    def build_arcs(self, mode, filter, filename):

        rpairs = set()
        reaction_rpair = dict()
        rpair_reaction = defaultdict(list)

        with open(join(self.path, filename + self.model.name + '.txt'), 'w') as f, \
                open(join(self.path, filename + self.model.name + '_names.txt'), 'w') as f2:
            try:
                model = self.model
                cofactors = set([c.id for c in self.cofactor_list])

                for reaction in model.reactions:
                    if reaction.is_substrate_product_distinct():
                        ri = self.__reaction_index(reaction.id, const.FOR_MATLAB)
                        reaction_rpair[reaction.id] = []
                        if re.search('-', reaction.id) is None:
                            kreaction = self.reactions[reaction.id]
                            for rclass in kreaction.rpairs(mode):
                                rpair = rclass
                                if kreaction.is_substrate(rclass[1]):
                                    rpair = rclass[::-1]

                                if filter & bool(set(rpair) & cofactors):
                                    continue

                                rpairs.add(rpair)
                                reaction_rpair[reaction.id].append(rpair)
                                m1 = self.__compound_index(rpair[0], const.FOR_MATLAB)
                                m2 = self.__compound_index(rpair[1], const.FOR_MATLAB)
                                f.write(str(m1) + '\t' + str(m2) + '\t' + str(ri) + '\t' + '1' + '\n')
                                f2.write(rpair[0] + '\t' + rpair[1] + '\t' + reaction.id + '\n')

                        else:
                            kreaction = self.reactions[re.search('R\d{5}', reaction.id).group(0)]
                            for rclass in kreaction.rpairs(mode):
                                rpair = rclass
                                if kreaction.is_substrate(rclass[0]):
                                    rpair = rclass[::-1]

                                if bool(set(rpair) & cofactors):
                                    continue
                                rpairs.add(rpair)
                                reaction_rpair[reaction.id].append(rpair)
                                m1 = self.__compound_index(rpair[0], const.FOR_MATLAB)
                                m2 = self.__compound_index(rpair[1], const.FOR_MATLAB)
                                f.write(str(m1) + '\t' + str(m2) + '\t' + str(ri) + '\t' + '1' + '\n')
                                f2.write(rpair[0] + '\t' + rpair[1] + '\t' + reaction.id + '\n')

                for r, rps in reaction_rpair.iteritems():
                    for rpair in rps:
                        rpair_reaction[rpair].append(r)

            except AttributeError as e:
                print(e)
                print(traceback.format_exc())
                print(sys.exc_info()[0])
                print('Error in arc creation.')
                return

        return rpairs, reaction_rpair, rpair_reaction

    def pickle_model(self, filename):
        with open(join(self.path, filename + ".pickle"), 'wb') as handle:
            cPickle.dump(self, handle, protocol=cPickle.HIGHEST_PROTOCOL)

    def save_arcgraph(self, filename):
        self.arcgraph.save(join(self.path, filename + "arcs.xml.gz"))

    def save_reactiongraph(self, filename):
        self.reactiongraph.save(join(self.path, filename + "reactions.graphml"))

    def write_metabolites(self, filename):
        with open(join(self.path, filename + self.model.name + '.txt'), 'w') as f1, \
                open(join(self.path, filename + '_names' + self.model.name + '.txt'), 'w') as f2:

            for id, metabolite in enumerate(self.model.metabolites, start=1):
                f1.write(metabolite.id + '\n')
                f2.write(metabolite.id + '\t' + metabolite.names[0] + '\t' + str(id) + '\n')

    def write_generic_compounds(self, filename):
        with open(join(self.path, filename + self.model.name + '.txt'), 'w') as f1, \
                open(join(self.path, filename + '_names' + self.model.name + '.txt'), 'w') as f2:
            for metabolite in self.generic_compound_list:
                m = self.__compound_index(metabolite.id, const.FOR_MATLAB)
                f1.write(str(m) + '\n')
                f2.write(metabolite.id + '\t' + metabolite.names[0] + '\n')

    def write_reactions(self, filename):
        with open(join(self.path, filename + self.model.name + '.txt'), 'w') as f1, \
                open(join(self.path, filename + '_names' + self.model.name + '.txt'), 'w') as f2, \
                open(join(self.path, filename + '_equations' + self.model.name + '.txt'), 'w') as f3:
            for reaction in self.model.reactions:
                f1.write(reaction.id + '\n')
                f2.write(reaction.id + '\t' + reaction.names[0] + '\n')
                f3.write(reaction.id + '\t' + reaction.equation + '\n')

        with open(join(self.path, filename + '_unique' + self.model.name + '.txt'), 'w') as f:
            for reaction in self.unique_reactions:
                f.write(reaction.id + '\t' + reaction.definition + '\n')

    def write_reversible_reactions(self, filename):
        with open(join(self.path, filename + self.model.name + '.txt'), 'w') as f1, \
                open(join(self.path, filename + '_names' + self.model.name + '.txt'), 'w') as f2:
            for reactions in self.reversible_reactions:
                index1 = self.__reaction_index(reactions[0], const.FOR_MATLAB)
                index2 = self.__reaction_index(reactions[1], const.FOR_MATLAB)
                f1.write(str(index1) + '\t' + str(index2) + '\n')
                f2.write(reactions[0] + '\t' + reactions[1] + '\n')

    def write_non_enzymatic_reactions(self, filename):
        with open(join(self.path, filename + self.model.name + '.txt'), 'w') as f1, \
                open(join(self.path, filename + '_names' + self.model.name + '.txt'), 'w') as f2:

            for reaction in self.model.reactions:
                r = self.reactions[re.search('R\d{5}', reaction.id).group(0)]
                if r.is_non_enzymatic():
                    f1.write(reaction.id + '\n')
                    f2.write(reaction.names[0] + '\n')

    def write_cofactors(self, filename):
        with open(join(self.path, filename + self.model.name + '.txt'), 'w') as f, open(join(self.path, filename
                                                                    + '_names' + self.model.name + '.txt'), 'w') as f2:
            for metabolite in self.cofactor_list:
                m = self.__compound_index(metabolite.id, const.FOR_MATLAB)
                f.write(str(m) + '\n')
                f2.write(metabolite.id + '\t' + metabolite.names[0] + '\t' + str(m) + '\n')

    def stoichiometric_matrix(self, filename):
        self.model.stoich_matrix.eliminate_zeros()
        sio.savemat(join(self.path, filename + self.model.name), {'S': self.model.stoich_matrix}, do_compression="True")

    def write_cytoscape_datatable(self):
        with open(self.model.name + '_data', 'w') as f:
            for reaction in self.model.reactions:
                f.write("%s : %s\n" % (reaction.id, "Blue"))

    def print_model(self):
        print("Reactions")
        print("---------")
        for x in self.model.reactions:
            print("%s : %s" % (x.id, x.name))

        print("Metabolites")
        print("-----------")
        for x in self.model.metabolites:
            print('%9s : %s' % (x.id, x.name))

    def filter_compounds_weight(self, minw, maxw):
        self.filtered_compounds_weight = \
            [compound for compound in self.model.metabolites if (minw <= compound.molweight <= maxw)]

    def filter_compounds_carbons(self, minc, maxc):
        for compound in self.model.metabolites:
            if hasattr(compound, 'formula'):
                numc = (re.findall('C((?![a-z])\d*)', compound.formula))
                if len(numc) > 1:
                    # more than one block of C in formula
                    continue
                if len(numc) == 0:
                    self.filtered_compounds_no_carbons.append(compound)
                    continue
                try:
                    n = int(numc[0])
                    if minc <= n <= maxc:
                        self.filtered_compounds_carbons.append(compound)
                except ValueError:
                    self.filtered_compounds_carbons.append(compound)
                    continue
            # no else case, because formula should always be set by the compound parser

    def basis_compounds(self, filename):
        id_list = self.__read_list(filename)
        basis_compounds = [compound for compound in self.model.metabolites if compound.id in id_list]
        self.basis_compound_list = basis_compounds

    def read_cofactors(self, filename):
        id_list = self.__read_list(filename)
        cofactors = [compound for compound in self.model.metabolites if compound.id in id_list]
        self.cofactor_list = cofactors

    def feasible_reactions(self):
        compounds = set([c.id for c in self.start_compound_list]) | set([c.id for c in self.cofactor_list]) | \
                    set([c.id for c in self.basis_compound_list]) | set([c.id for c in self.generic_compound_list])

        num_new_compounds = len(compounds)
        model_reactions = list(self.model.reactions)
        reactions = set()

        while num_new_compounds > 0:
            num_new_compounds = 0
            for r in model_reactions[:]:
                subs = set(r.substrates())
                if subs <= compounds:
                    reactions.add(r.id)
                    model_reactions.remove(r)
                    products = set(r.products())
                    num_new_compounds += len(products - compounds)
                    compounds = compounds.union(products)

        self.feasible_reaction_list = reactions

        reacs = [self.reactions[r[1:]].reverse() if r.startswith('-') else self.reactions[r] for r in reactions]

        with open(join(self.path, 'feasible_reactions.txt'), 'w') as f:
            for r in reacs:
                f.write(r.id + '\t' + r.equation + '\n')

        with open(join(self.path, 'substrates_products.txt'), 'w') as f:
            for compound in compounds:
                f.write(compound + '\t' + self.compounds[compound].names[0] + '\n')

    def __start_compounds(self):
        scomps = list(set(self.filtered_compounds_weight) - set(self.filtered_compounds_no_carbons) -
                      set(self.cofactor_list) - set(self.generic_compound_list) - set(self.filtered_compounds_carbons))

        ids = [comp.id for comp in scomps]
        start = [compound for k, compound in self.compounds.iteritems() if compound.id in ids]
        return start

    def write_start_compounds(self, filename):
        with open(join(self.path, filename + self.model.name + '.txt'), 'w') as f1, \
                open(join(self.path, filename + '_names' + self.model.name + '.txt'), 'w') as f2:
            for metabolite in self.start_compound_list:
                m = self.__compound_index(metabolite.id, const.FOR_MATLAB)
                f1.write(str(m) + '\n')
                f2.write(metabolite.id + '\t' + metabolite.names[0] + '\n')

    def write_basis_compounds(self, filename):
        with open(join(self.path, filename + self.model.name + '.txt'), 'w') as f:
            for compound in self.basis_compound_list:
                if compound in self.model.metabolites:
                    m = self.__compound_index(compound.id, const.FOR_MATLAB)
                    f.write(str(m) + '\n')

    def build_targets(self, filename):
        targets = set(self.model.metabolites) - set(self.basis_compound_list) - set(self.generic_compound_list) \
                      - set(self.cofactor_list) - set(self.start_compound_list)

        targets_with_start = set(self.model.metabolites) - set(self.basis_compound_list) \
                                    - set(self.generic_compound_list) - set(self.cofactor_list)

        producible = set([self.compounds[x[1]] for x in self.rpairs])
        # all targets in model that are potentially producible
        self.all_target_compounds = targets & producible
        self.targets_that_are_startmetabolites = (targets_with_start - targets) & producible

        with open(join(self.path, 'all_targets.txt'), 'w') as f:
            for metabolite in self.all_target_compounds:
                f.write(metabolite.id + '\t' + metabolite.names[0] + '\n')

        with open(join(self.path, filename + '.txt'), 'w') as f1,\
                open(join(self.path, filename + '_names' + self.model.name + '.txt'), 'w') as f2:
            for metabolite in self.all_target_compounds:
                f1.write(metabolite.id + '\n')
                f2.write(metabolite.id + '\t' + metabolite.names[0] + '\n')

    def write_thermodynamics(self, filename):
        sio.savemat(join(self.path, filename + self.model.name), {'thermodynamics': self.thermodynamics})

    def write_rpairs(self, filename):
        with open(join(self.path, filename + self.model.name + '.txt'), 'w') as f:
            for rp in self.rpairs:
                f.write(rp[0] + '\t' + rp[1] + '\n')

    def compound_is_in_arc(self, compound_id):
        tmp = [item for item in self.rpairs if compound_id in item]
        return len(tmp) > 0

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

    def reaction_is_in_host(self, reaction_id):
        reaction = self.reactions[reaction_id]

        for e in reaction.enzymes:
            enzyme = self.enzymes['EC ' + e]
            if enzyme.is_in_organism(self.host):
                return True
        return False

    def build_metabolite_pool(self, basis=False):
        start = self.__start_compounds()
        if basis:
            assert len(self.basis_compound_list) > 0
            self.start_compound_list = self.basis_compound_list
            self.long_start_list = start

        else:
            self.start_compound_list = start
        self.__generic_compounds()

    def write_metabolite_lists(self):

        """TODO: check for empty lists and print a warning"""
        self.write_metabolites('compounds')
        self.write_cofactors('cofactors')
        self.write_generic_compounds('generic_compounds')
        self.write_non_enzymatic_reactions('non_enzymatic_reactions')
        self.write_start_compounds('start_compounds')
        self.write_basis_compounds('terminal_compounds')

    def __populate(self, reacs, filterreacs, with_thermodynamics=True):

        metabsset = set()
        reactions = []

        self.all_reactions = reacs

        if filterreacs:
            rea = [r for r in reacs if r.has_valid_stoichiometry() and r.has_valid_reactants()
                   and not r.is_generic() and r.has_rclass() and not self.reaction_has_generic_compound(r.id)]
        else:
            rea = [r for r in reacs if r.has_valid_stoichiometry() and r.has_valid_reactants() and not r.is_general()
                   and not self.reaction_has_generic_compound(r.id)]

        self.unique_reactions = rea #for the statistics: reactions after filtering without reverse

        for r in rea:
            for mId, _ in r.reactants.iteritems():
                metabsset.add(mId)

            reverse_reaction = r.reverse()

            if with_thermodynamics:
                if r.id in self.thermodynamics:
                    dg0, dgu, dgm, direction = self.thermodynamics[r.id]

                    if direction == -1:
                        reactions.append(reverse_reaction)

                    elif direction == 1:
                        reactions.append(r)

                    elif direction == 0:
                        reactions.append(r)
                        reactions.append(reverse_reaction)
                        self.reversible_reactions.append((r.id, reverse_reaction.id))

                    self.reaction_reversibilities[r.id] = direction

                else:
                    reactions.append(r)
                    self.reaction_reversibilities[r.id] = 1

            else:
                ## all irreversible
                reactions.append(r)
                self.reaction_reversibilities[r.id] = 1

                ## all reversible
                # reactions.append(r)
                # reactions.append(reverse_reaction)
                # self.reversible_reactions.append((r.id, reverse_reaction.id))
                # self.reaction_reversibilities[r.id] = 0



        stoch_mat = np.zeros((len(metabsset), len(reactions)), dtype='double')
        metabs = sorted(metabsset)
        c_index = {}
        for i, c in enumerate(metabs):
            c_index[c] = i
        for i, r in enumerate(reactions):
            for mId, stoich in r.reactants.iteritems():
                stoch_mat[c_index[mId], i] = stoich

        try:
            self.model = \
                Model(reactions=reactions,
                      metabolites=[self.compounds[m] if m in self.compounds else kCompound("ENTRY" + m) for m in metabs],
                      stoich_matrix=sparse.csr_matrix(stoch_mat, dtype='double'), name=self.host)

        except AttributeError:
            self.model = Model([], [], None, name=self.host)

    def __filter_compounds_arcs(self):
        for metabolite in self.model.metabolites:
            try:
                if not self.compound_is_in_arc(metabolite.id):
                    self.compounds_no_arcs.append(self.compounds[metabolite.id])
            except KeyError:
                continue

    def __generic_compounds(self):
        generic = [compound for compound in self.model.metabolites if compound.is_generic()]
        self.generic_compound_list = generic

    def compounds_atom_types(self):
        with open(join(self.path, 'compounds_atomtypes.txt'), 'w') as f:
            for metabolite in self.model.metabolites:
                f.write('{}\t{}\t{}\n'.format(metabolite.id, metabolite.names[0], len(metabolite.atomtypes)))

    # def build_reaction_only_graph(self, filename):
    #     g = gt.Graph()
    #
    #     cofactors = [s.id for s in self.cofactor_list]
    #     generic_compound = [s.id for s in self.generic_compound_list]
    #
    #     edges = []
    #     for reaction1 in self.model.reactions:
    #         products = reaction1.products()
    #         for reaction2 in self.model.reactions:
    #             substrates = reaction2.substrates()
    #
    #             inter = set(products) & set(substrates)
    #             inter2 = inter - set(cofactors) - set(generic_compound)
    #             if bool(inter2):
    #                 edges.append((reaction1.id, reaction2.id))
    #
    #     ids = g.add_edge_list(edges, hashed=True, string_vals=True)
    #
    #     g.vertex_properties["ids"] = ids
    #     g.vertex_properties["names"] = g.new_vertex_property("string")
    #     g.vertex_properties["color"] = g.new_vertex_property("string")
    #     g.vertex_properties["size"] = g.new_vertex_property("int")
    #     g.edge_properties["arrows"] = g.new_edge_property("string")
    #
    #     for v in g.vertices():
    #         v_id = g.vp.ids[v]
    #         g.vp.names[v] = v_id
    #         g.vp.size[v] = 10
    #
    #     for e in g.edges():
    #         g.ep.arrows[e] = "none"
    #
    #     deg = g.degree_property_map("in")
    #     deg.a = 4 * (deg.a * 0.5 + 0.4)+10
    #     ebet = gt.betweenness(g)[1]
    #     ebet.a /= ebet.a.max() / 10.
    #     eorder = ebet.copy()
    #     eorder.a *= -1
    #
    #     root = gt.find_vertex(g, g.vertex_properties["names"], "R00289")[0]
    #     pos = gt.radial_tree_layout(g, g.vertex_index[root])
    #     control = g.new_edge_property("vector<double>")
    #     for e in g.edges():
    #         d = math.sqrt(sum((pos[e.source()].a - pos[e.target()].a) ** 2)) / 5
    #         control[e] = [0.3, d, 0.7, d]
    #     gt.graph_draw(g, pos=pos, vertex_size=deg, vertex_fill_color=deg,edge_pen_width=ebet, vorder=deg, eorder = eorder,
    #     edge_control_points = control,
    #     vertex_text=g.vp.names,
    #     output_size=(1000, 1000),
    #     output = filename)
    #
    #     #root = gt.find_vertex(g, g.vertex_properties["names"], "R00289")[0]
    #     #pos = gt.radial_tree_layout(g, g.vertex_index[root])
    #     #pos = gt.arf_layout(g, max_iter=10000, dt=1e-4)
    #     #pos = gt.sfdp_layout(g, multilevel=True)
    #
    #     #gt.graph_draw(g, pos,
    #      #             vertex_size=g.vp.size,
    #       #            vertex_text=g.vp.names,
    #        #           edge_pen_width=0.5,
    #         #          edge_end_marker=g.ep.arrows, output=filename, fit_view=True, output_size=(10000, 10000))

    def build_result_graph(self, reactionlist, product, filename):
        g = gt.Graph()

        cofactors = [s.id for s in self.cofactor_list]
        generic_compound = [s.id for s in self.generic_compound_list]

        edges = []
        for reaction in self.model.reactions:
            for c, val in reaction.reactants.iteritems():
                if c not in cofactors and c not in generic_compound:
                    if val < 0:
                        edges.append((c, reaction.id))
                    else:
                        edges.append((reaction.id, c))

        ids = g.add_edge_list(edges, hashed=True, string_vals=True)

        g.vertex_properties["ids"] = ids
        g.vertex_properties["names"] = g.new_vertex_property("string")
        g.vertex_properties["color"] = g.new_vertex_property("string")
        g.vertex_properties["size"] = g.new_vertex_property("int")
        g.edge_properties["arrows"] = g.new_edge_property("string")

        for v in g.vertices():
            v_id = g.vp.ids[v]
            g.vp.names[v] = v_id
            g.vp.color[v] = "red" if v_id in reactionlist else "white"
            g.vp.size[v] = 10 if v_id in reactionlist else 3

        for e in g.edges():
            g.ep.arrows[e] = "none"

        root = gt.find_vertex(g, g.vertex_properties["names"], product)[0]
        pos = gt.radial_tree_layout(g, g.vertex_index[root])

        gt.graph_draw(g, pos,
                      vertex_size=g.vp.size,
                      vertex_text=g.vp.names,
                      vertex_fill_color=g.vp.color,
                      edge_pen_width=0.5,
                      edge_end_marker=g.ep.arrows, output=filename, fit_view=True, output_size=(10000, 10000))

    def build_result_graph_in_whole(self, reactionlist, filename):
        g = gt.Graph()

        reactions = {r.id: r for r in self.model.reactions}

        edges = []
        for reaction in self.model.reactions:
            for c, val in reaction.reactants.iteritems():
                if val < 0:
                    edges.append((c, reaction.id))
                else:
                    edges.append((reaction.id, c))

        ids = g.add_edge_list(edges, hashed=True, string_vals=True)

        g.vertex_properties["ids"] = ids
        g.vertex_properties["names"] = g.new_vertex_property("string")
        g.vertex_properties["color"] = g.new_vertex_property("string")
        g.vertex_properties["size"] = g.new_vertex_property("int")
        g.edge_properties["arrows"] = g.new_edge_property("string")

        for v in g.vertices():
            v_id = g.vp.ids[v]
            g.vp.names[v] = v_id
            g.vp.color[v] = "red" if v_id in reactionlist else "lightgrey" if v_id in reactionlist else "white"
            g.vp.size[v] = 10 if v_id in reactions else 3

        for e in g.edges():
            g.ep.arrows[e] = "none"

        pos = gt.sfdp_layout(g, multilevel=True)

        gt.graph_draw(g, pos,
                      vertex_size=g.vp.size,
                      vertex_fill_color=g.vp.color,
                      edge_pen_width=0.5,
                      edge_end_marker=g.ep.arrows, output=filename, fit_view=True, output_size=(10000, 10000))

    def build_reaction_graph(self):
        g = self.reactiongraph
        reactions = {r.id: r for r in self.model.reactions}
        metabolites = {m.id: m for m in self.model.metabolites}

        edges = []
        for reaction in self.model.reactions:
            for c, val in reaction.reactants.iteritems():
                if val < 0:
                    edges.append((c, reaction.id))
                else:
                    edges.append((reaction.id, c))

        ids = g.add_edge_list(edges, hashed=True, string_vals=True)

        # VERTEX PROPERTIES

        # assign compound ids to vertices
        g.vertex_properties["ids"] = ids
        g.vertex_properties["names"] = g.new_vertex_property("string")
        g.vertex_properties["shapes"] = g.new_vertex_property("string")
        g.vertex_properties["size"] = g.new_vertex_property("int")
        g.vertex_properties["color"] = g.new_vertex_property("string")

        cofactors = [s.id for s in self.cofactor_list]
        excluded = [s.id for s in self.excluded_compounds]
        start_compounds = [s.id for s in self.start_compound_list]

        for v in g.vertices():
            v_id = g.vp.ids[v]
            g.vp.names[v] = metabolites[v_id].names[0] if v_id.startswith('C') else reactions[v_id].id
            g.vp.shapes[v] = "octagon" if v_id in cofactors else "circle"
            g.vp.size[v] = 10 if v_id in cofactors or v_id in excluded else 20
            g.vp.color[v] = "lightgrey" if v_id in excluded else "lightgrey" if v_id in cofactors else \
                "yellow" if v_id in start_compounds else "red" if v_id == const.PRODUCT else "floralwhite"

        g.edge_properties["arrows"] = g.new_edge_property("string")

        for e in g.edges():
            # source = e.source()
            target = e.target()
            eid = g.vp.ids[target]
            if eid.startswith('-') | eid.startswith('R'):
                g.ep.arrows[e] = "none"
                g.vp.shapes[target] = "pie"
                g.vp.color[target] = "royalblue"
            else:
                g.ep.arrows[e] = "arrow"

    def build_arc_graph(self, rpairs):
        g = gt.Graph()
        metabolites = {m.id: m for m in self.model.metabolites}

        # v: compounds, e: arcs
        cids = g.add_edge_list(rpairs, hashed=True, string_vals=True)

        # VERTEX PROPERTIES

        # assign compound ids to vertices
        g.vertex_properties["compound_ids"] = cids
        g.vertex_properties["name"] = g.new_vertex_property("string")
        g.vertex_properties["is_start_compound"] = g.new_vertex_property("bool")
        g.vertex_properties["is_basis_compound"] = g.new_vertex_property("bool")
        g.vertex_properties["is_generic_compound"] = g.new_vertex_property("bool")
        g.vertex_properties["is_cofactor"] = g.new_vertex_property("bool")
        g.vertex_properties["molweight"] = g.new_vertex_property("double")
        g.vertex_properties["color"] = g.new_vertex_property("string")

        start_compounds = [s.id for s in self.start_compound_list]
        cofactors = [s.id for s in self.cofactor_list]
        generic_compound = [s.id for s in self.generic_compound_list]

        for v in g.vertices():
            v_id = g.vp.compound_ids[v]
            g.vp.is_start_compound[v] = v_id in start_compounds
            g.vp.is_basis_compound[v] = v_id in self.basis_compound_list
            g.vp.is_cofactor[v] = v_id in cofactors
            g.vp.is_generic_compound[v] = v_id in generic_compound
            g.vp.molweight[v] = self.compounds[v_id].molweight
            g.vp.name[v] = metabolites[v_id].names[0]
            g.vp.color[v] = "lightgrey" if v_id in cofactors else "yellow" if v_id in start_compounds else "floralwhite"

        # EDGE PROPERTIES

        # assign rpairs to edges
        g.edge_properties["rpairs"] = g.new_edge_property("string")

        # assign reactions containing the respective rpair
        g.edge_properties["reaction_ids"] = g.new_edge_property("object")

        for e in g.edges():
            source_id = g.vp.compound_ids[e.source()]
            target_id = g.vp.compound_ids[e.target()]
            rpair = (source_id, target_id)
            reaction_ids = self.rpair_reaction[rpair]
            g.ep.reaction_ids[e] = reaction_ids
            g.ep.rpairs[e] = rpair

        g.graph_properties["host"] = g.new_graph_property("string")  # according to the documentation this is necessary
        g.gp.host = self.host
        return g

    @staticmethod
    def read_pathways_for_color(filename):
        pathway_list = []
        with open(join(const.DATAPATH_KEGG, filename), 'r') as f:
            for line in f:
                pathway_list.append(re.search('(?:(\d{5}))', line).group(0))
        return pathway_list

    def write_input_file(self):
        with open(join(const.MECAT_BASE, 'input_{}.txt'.format(self.name)), 'w') as f:
            f.write('DATA_PATH\t{}\n'.format(self.name))
            f.write('GENERIC_METABOLITES\tgeneric_compounds{}.txt\n'.format(self.host))
            f.write('COFACTORS\tcofactors{}.txt\n'.format(self.host))
            f.write('COFACTOR_NAMES\tcofactors_names{}.txt\n'.format(self.host))
            f.write('START_METABOLITES\tstart_compounds{}.txt\n'.format(self.host))
            f.write('BASIS_METABOLITES\tterminal_compounds{}.txt\n'.format(self.host))
            f.write('START_METABOLITE_NAMES\tstart_compounds_names{}.txt\n'.format(self.host))
            f.write('COMPOUNDS\tcompounds{}.txt\n'.format(self.host))
            f.write('REACTIONS\treactions{}.txt\n'.format(self.host))
            f.write('COMPOUND_NAMES\tcompounds_names{}.txt\n'.format(self.host))
            f.write('REACTION_NAMES\treactions_names{}.txt\n'.format(self.host))
            f.write('ARCS\t{}{}.txt\n'.format('RxR_kegg', self.host))
            f.write('ARC_NAMES\t{}{}_names.txt\n'.format('RxR_kegg', self.host))
            f.write('STOICHIOMETRIC_MATRIX\tS{}.mat\n'.format(self.host))
            f.write('REVERSIBLE_REACTIONS\treversible_reactions{}.txt\n'.format(self.host))
            f.write('KEGG_DATA_PATH\t{}\n'.format(join('..', 'KEGG_data', '030619')))
            f.write('REACTION_MAP\t{}\n'.format('KEGG_REACTION.mat'))
            f.write('COMPOUND_MAP\t{}\n'.format('KEGG_COMPOUND.mat'))
            f.write('ENZYME_MAP\t{}\n'.format('KEGG_ENZYME.mat'))
            f.write('PATHWAY_MAP\t{}\n'.format('KEGG_PATHWAY.mat'))
            f.write('THERMODYNAMICS_MAP\tthermodynamics{}.mat\n'.format(self.host))
            f.write('NON_ENZYMATIC_REACTIONS\t{}{}.txt\n'.format('non_enzymatic_reactions', self.host))
            f.write('REVERSIBILITIES_MAP\t{}\n'.format('reversibilities_map.mat'))
            f.write('HOST\t{}\n'.format(self.host))

    @staticmethod
    def __read_list(filename):
        id_list = []
        with open(filename, 'r') as f:
            for line in f:
                id_list.append(re.search('(?:C|R)[0-9]{5}', line).group(0))
        return id_list

    def __compound_index(self, cid, for_matlab):
        if 'cIndex' not in self.__dict__:
            self.cIndex = {}
            for i, c in enumerate(self.model.metabolites):
                self.cIndex[c.id] = i
        index = self.cIndex[cid]
        return index + 1 if for_matlab else index

    def __reaction_index(self, rid, for_matlab):
        if 'rIndex' not in self.__dict__:
            self.rIndex = {}
            for i, r in enumerate(self.model.reactions):
                self.rIndex[r.id] = i
        index = self.rIndex[rid]
        return index + 1 if for_matlab else index
