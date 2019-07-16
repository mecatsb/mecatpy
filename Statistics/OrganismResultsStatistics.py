import itertools
from cPickle import load
import os
from os.path import join
import constants as const
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import cPickle
import graph_tool.all as gt
import numpy as np
import helper
import re
import pandas as pd
from matplotlib_venn import venn3_unweighted, venn3_circles


class OrganismResultStatistics:

    def __init__(self):

        with open(join(const.DATAPATH_KEGG, "KEGG_organism.pickle"), 'rb') as handle:
            self.organisms = load(handle)

        with open(join(const.DATAPATH_KEGG, "KEGG_reaction.pickle"), 'rb') as handle:
            self.reactions = load(handle)

        with open(join(const.DATAPATH_KEGG, "KEGG_compound.pickle"), 'rb') as handle:
            self.compounds = load(handle)

        with open(join(const.DATAPATH_KEGG, "KEGG_enzyme.pickle"), 'rb') as handle:
            self.enzymes = load(handle)

        with open(join(const.DATAPATH_KEGG, "KEGG_pathway.pickle"), 'rb') as handle:
            self.pathways = load(handle)

        with open(join(const.DATAPATH_KEGG, "KEGG_atom_types.pickle"), 'rb') as handle:
            self.atom_types = load(handle)

        with open (join(const.MECAT_BASE, 'organisms.txt'), 'r') as f:
            organisms = [x.strip().split('\n')[0] for x in f.readlines()]

        self.suffix = ['basis']
        self.o = self.organisms
        self.organisms = list()
        for o in organisms:
            for s in self.suffix:
                self.organisms.append(o+s)

        self.results = dict()
        self.result_counts = dict()
        self.target_directories = dict()
        self.targets = dict()
        self.targets_that_are_startmetabolites = dict()  # long target list for the basis models
        self.found = dict()
        self.not_found = dict()
        self.really_not_found = dict()
        self.really_really_not_found = dict()
        self.bfs_reachable = dict()
        self.bfs_not_reachable = dict()
        self.builder = dict()
        self.predicted_feasibility = dict()
        self.targets_in_satellites = dict()
        self.__load_builder()

    def organism_statistics(self):
        targets = list()
        target_lists = dict()

        for o in self.organisms:
            with open(join(const.MECAT_BASE, 'Targets '+o, 'targets_overview'), 'r') as handle:
                content = handle.readlines()
            all_targets = [x.id for x in self.builder[o].all_target_compounds]

            content = [x for x in (x.strip().split('\t') for x in content) if x[0] in all_targets]
            res = [(int(c[2]), c[1], c[0]) for c in content]
            self.results[o] = res
            self.target_directories[o] = [x[0] for x in os.walk(join(const.MECAT_BASE, 'Targets'+o)) if 'results' in x[0]]
            not_found = [id for c, _, id in res if c is 0]
            found = [id for c, _, id in res if c is 1]

            self.found[o] = found
            self.not_found[o] = not_found
            self.targets[o] = set([id for c, _, id in res])
            self.targets_that_are_startmetabolites[o] = [x.id for x in self.builder[o].targets_that_are_startmetabolites]

            # from ModelStatistics
            cnames = dict()
            with open(join(const.MECAT_BASE, o, "Statistik", "targets_in_satellites.txt"), 'r') as handle:
                content = handle.readlines()
            content = [x.strip().split('\t') for x in content]
            for c in content:
                if len(c) == 1:
                    cnames[c[0]] = c[0]
                else:
                    cnames[c[0]] = c[1]
            self.targets_in_satellites[o] = [c[0] for c in content]

            targets.append(set(found))
            target_lists[o] = found

            with open(join(const.MECAT_BASE, o, "Statistik", "compounds_components.txt"), 'r') as handle:
                compounds_arcs = handle.readlines()
            compounds_arcs = [x.strip().split('\t')[1] for x in compounds_arcs]

            nicht_gefunden = (set(not_found) & set(compounds_arcs)) - set(self.targets_in_satellites[o])

            self.__write_targets(nicht_gefunden, join(const.MECAT_BASE, 'Targets '+o, 'not_found.txt'), join(const.MECAT_BASE, 'Targets '+o, 'not_found_atom_types.txt'))
            self.__write_targets(found, join(const.MECAT_BASE, 'Targets ' + o, 'found.txt'), join(const.MECAT_BASE, 'Targets ' + o, 'found_atom_types.txt'))
            self.__write_found_but_not_predicted(o, found, join(const.MECAT_BASE, 'Targets ' + o, 'found_not_predicted.txt'))
            self.__write_predicted_but_not_found(o, found, join(const.MECAT_BASE, 'Targets ' + o, 'predicted_not_found.txt'))


            # annotate pathways
            for target in self.target_directories:

                resultdirectories = [join(const.MECAT_BASE, 'Targets ' + target, x) for x in os.listdir(join(const.MECAT_BASE, 'Targets '+ target)) if os.path.isdir(join(const.MECAT_BASE, 'Targets '+ target, x))]

                for target in resultdirectories:
                    allpathways = []
                    reactions = set()
                    reactionfiles = [join(target, x) for x in os.listdir(target) if x.endswith('reactions')]
                    if not reactionfiles:
                        continue
                    reactionfile = reactionfiles[0]
                    with open(reactionfile) as f:
                        for i, l in enumerate(f):
                            reactionids = l.strip().split(' ')
                            for r in reactionids:
                                rid = re.findall('R\d{5}', r)[0]
                                reactions.add(rid)

                    with open(join(target, "pathways.txt"), 'w') as f:
                        for r in reactions:
                            reaction = self.reactions[r]
                            pathways = reaction.pathways
                            for p in pathways:
                                allpathways.append(p)
                            pathway_names = [self.pathways[p].name for p in pathways]
                            pws = '\t'.join(pathway_names)
                            f.write(reaction.id + '\t' + reaction.names[0] + '\t' + pws + '\n')

        self.__write_common_targets(targets)
        self.__write_unique_targets(target_lists)

    def __load_builder(self):
        for organism in self.organisms:
            with open(join(const.MECAT_BASE, organism, "KEGG_model.pickle"), 'rb') as handle:
                builder = cPickle.load(handle)
                self.builder[organism] = builder

            feasible = set()
            with open(join(const.MECAT_BASE, organism, "Statistik", "target_is_feasible.txt"), 'rb') as handle:
                for i, l in enumerate(handle):
                    feasible.add(re.findall('C\d{5}', l)[0])

            self.predicted_feasibility[organism] = feasible

    # find and write targets that are found in all organisms under consideration
    def __write_common_targets(self, targetlists):
        common_targets = set.intersection(*targetlists)

        with open(join(const.MECAT_BASE, 'common_targets_in_organisms.txt'), 'w') as h:
            for target in common_targets:
                compound = self.compounds[target]
                h.write('{}\t{}\n'.format(target, compound.names[0]))

    def __write_unique_targets(self, targetlists):
        organisms = [o for o in targetlists if o != 'kegg']
        for organism, targetlist in (x for x in targetlists.iteritems() if x[0] != 'kegg'):
            found = self.found[organism]
            other_targets = [targetlists[o] for o in organisms if o is not organism]
            combined_targets = [item for sublist in other_targets for item in sublist]

            unique_targets = sorted((set(targetlist) - set(combined_targets)) & set(found))

            with open(join(const.MECAT_BASE, 'unique_targets_in_{}.txt'.format(organism)), 'w') as h:
                for target in unique_targets:
                    compound = self.compounds[target]
                    h.write('{}\t{}\n'.format(target, compound.names[0]))

    @staticmethod
    def __write_atomtype_statistics(functional_groups, heterocycles, file):
        with open(file, 'w') as handle:
            handle.write("functional groups:\n".format())
            for fg, count in functional_groups.iteritems():
                handle.write("{}\t{}\n".format(fg, str(count)))
            handle.write("\nheterocycles:\n".format())
            for hc, count in heterocycles.iteritems():
                handle.write("{}\t{}\n".format(hc, str(count)))

    def __write_targets(self, targetlist, targetfile, statisticsfile):
        functional_group_count = dict()
        heterocycle_count = dict()
        with open(targetfile, 'w') as handle :
            for id in sorted(targetlist):
                compound = self.compounds[id]
                heterocycles = compound.heterocycles
                functional_groups = compound.functional_groups

                for functional_group in functional_groups:
                    if functional_group in functional_group_count:
                        functional_group_count[functional_group] = functional_group_count[functional_group]+1
                    else:
                        functional_group_count[functional_group] = 1

                for heterocycle in heterocycles:
                    if heterocycle in heterocycle_count:
                        heterocycle_count[heterocycle] = heterocycle_count[heterocycle]+1
                    else:
                        heterocycle_count[heterocycle] = 1

                handle.write('{}\t{}\n'.format(id, compound.names[0]))
                self.__write_atomtype_statistics(functional_group_count, heterocycle_count, statisticsfile)

    def __write_found_but_not_predicted(self, o, targetlist, targetfile):
        found_not_predicted = set(targetlist) - self.predicted_feasibility[o]
        with open(targetfile, 'w') as handle:
            for id in found_not_predicted:
                compound = self.compounds[id]
                handle.write('{}\t{}\n'.format(id, compound.names[0]))

    def __write_predicted_but_not_found(self, o, targetlist, targetfile):
        predicted_not_found = self.predicted_feasibility[o] - set(targetlist)
        with open(targetfile, 'w') as handle:
            for id in predicted_not_found:
                compound = self.compounds[id]
                handle.write('{}\t{}\n'.format(id, compound.names[0]))

    def __dump_statistics(self, organisms, filename):

        found_not_predicted = {o:set(self.found[o]) -  self.predicted_feasibility[o] for o in organisms}
        kein_pfad = {o:set(self.not_found[o]) - set(self.really_not_found[o]) - set(self.targets_in_satellites[o]) for o in organisms}

        for k,v in kein_pfad.iteritems():
            with open(join(const.MECAT_BASE, k + ' kein_pfad.txt'), 'w') as f:
                for id in v:
                    compound = self.compounds[id]
                    f.write(id + '\t' +  compound.names[0] + '\n')

        combined = [(len(self.targets[o]),
                     len(found_not_predicted[o]),
                     len(self.found[o]) - len(found_not_predicted[o]),
                     len(self.not_found[o]) - len(self.really_not_found[o]) - len(kein_pfad[o]),
                     len(kein_pfad[o]),
                     len(self.really_not_found[o]) - len(self.really_really_not_found[o]),
                     len(self.really_really_not_found[o])) for o in organisms ]

        found_not_predicted = np.array([i[1] for i in combined],dtype=float)
        found = np.array([i[2] for i in combined],dtype=float)
        not_found_sat = np.array([i[3] for i in combined],dtype=float)
        not_found_main = np.array([i[4] for i in combined], dtype=float)
        not_found_feasibility = np.array([i[5] for i in combined],dtype=float)
        not_found_other = np.array([i[6] for i in combined],dtype=float)

        all2 = np.stack((found_not_predicted, found, not_found_sat, not_found_main, not_found_feasibility, not_found_other))
        legend_text = ('found and not predicted', 'found and predicted', 'not found (connectivity - satellite components)', 'not found (connectivity - components with start metabolites)', 'not found (feasibility)', 'not found (other)')
        organism_labels = list(map(lambda x: x.replace('basis', ''), organisms))
        all3 = pd.DataFrame(data=all2, index=legend_text, columns=organism_labels)
        all3.to_csv(filename, float_format='%d')

    def neighbor_graphs(self):
        for o in self.organisms:

            builder = self.builder[o]
            cofactors = set([c.id for c in builder.cofactor_list])

            g = self.builder[o].arcgraph
            model_reactions = {r.id: r for r in self.builder[o].model.reactions}

            has_start_metabolite_in_graph = list()

            for nf in self.not_found[o]:
                g.set_reversed(True)
                root = gt.find_vertex(g, g.vertex_properties["compound_ids"], nf)[0]

                helper.add_properties(g)
                g.vp.size[root] = 100
                old_color = g.vp.color[root]
                g.vp.color[root] = "red"
                g.vp.besucht[root] = True
                visitor = helper.NeighborVisitor(g)
                gt.bfs_search(g,root,visitor)

                g.set_reversed(False)
                g.set_vertex_filter(g.vp.besucht)
                g.vp.size.a = 90 * 0.9**g.vp.dist.a+10
                if any([v for v in g.vertices() if g.vp.color[v] == "yellow"]):
                    has_start_metabolite_in_graph.append(nf)

                # reset graph
                g.clear_filters()
                helper.remove_properties(g)
                g.vp.color[root] = old_color

            feasible_reactions = self.builder[o].feasible_reaction_list

            self.really_not_found[o] = has_start_metabolite_in_graph

            with open(join(const.MECAT_BASE, 'Targets '+o,'has_start_metabolite_in_graph.txt'),'w') as f1, \
                    open(join(const.MECAT_BASE, 'Targets ' + o, 'reactions_for_really_not_found.txt'), 'w') as f2:
                rrnf = list()

                substrates_not_found = set(has_start_metabolite_in_graph + self.not_found[o] + self.really_not_found[o])

                for c in has_start_metabolite_in_graph:
                    f1.write('{}\t{}\n'.format(c, self.compounds[c].names[0]))

                    reactions_for_target = self.__reactions_for_target(c, model_reactions)

                    fr = set(reactions_for_target) & feasible_reactions
                    if len(fr) > 0:
                        rrnf.append(c)

                    f2.write('{}\t{}:\n'.format(c, self.compounds[c].names[0]))
                    for reac_id in fr:
                        t = ''
                        # test if the substrates of the reaction are all cofactors
                        substrates = set(next(itertools.ifilter(lambda x:x.id == reac_id, builder.model.reactions)).substrates())
                        if substrates <= cofactors:
                            t += ' keine Kante (Cofactors)'

                        # test if the reaction has arcs that end with the target
                        rpairs = builder.reaction_rpair[reac_id]
                        if not any((x[1] == c) for x in rpairs):
                                t += ' keine Kante zu Target'

                        #  test if any substrate is also in list of not found
                        if substrates & substrates_not_found:
                                t += ' Substrat auch nicht gefunden'

                        fx = re.findall('R\d{5}', reac_id)[0]
                        f2.write('\t{}\t{}\t{}\n'.format(t, reac_id, self.reactions[fx].definition))

            self.really_really_not_found[o] = rrnf
        self.__dump_statistics(self.organisms, join(const.MECAT_BASE, 'targets_in_organisms.csv'))

    def __reactions_for_target(self, compound_id, model_reactions):

        reactions = set()
        reaction_ids = self.compounds[compound_id].reactions

        for r in reaction_ids:
            if r in model_reactions:
                reac = model_reactions[r]
                if compound_id in reac.products():
                    reactions.add(r)

            rev = '-'+r
            if rev in model_reactions:
                reac = model_reactions[rev]
                if compound_id in reac.products():
                    reactions.add(rev)
        return reactions

    def reachability_graphs(self):
        for o in self.organisms:
            g = self.builder[o].arcgraph
            targets = self.targets[o]

            root = g.add_vertex()
            for v in g.vertices():
                if g.vp.is_start_compound[v]:
                    g.add_edge(root, v)

                if g.vp.compound_ids[v] in targets:
                    g.vp.color[v] = 'green'

            helper.add_properties(g)
            g.vp.size[root] = 100
            old_color = g.vp.color[root]
            g.vp.color[root] = "red"
            g.vp.besucht[root] = True
            visitor = helper.NeighborVisitor(g)
            gt.bfs_search(g, root, visitor)

            pos = gt.radial_tree_layout(g, g.vertex_index[root], r=0.001)
            g.vp.besucht[root] = False
            g.set_vertex_filter(g.vp.besucht)
            g.vp.size.a = 50 * 0.9 ** g.vp.dist.a + 10

            gt.graph_draw(g, pos, vertex_size=g.vp.size,
                          vertex_text=g.vp.dist, vertex_text_position=-0.1, vertex_fill_color=g.vp.color,
                          vertex_text_color="black", edge_marker_size=10, edge_pen_width=0.5, output_size=(5000, 5000),
                          output=join(const.MECAT_BASE, 'Targets ' + o, 'reachability.pdf'))

            bfs_reachable = list()

            for v in g.vertices():
                if g.vp.compound_ids[v] in self.targets[o]:
                    bfs_reachable.append(g.vp.compound_ids[v])

            bfs_not_reachable = self.targets[o] - set(bfs_reachable)

            self.bfs_reachable[o] = bfs_reachable
            self.bfs_not_reachable[o] = bfs_not_reachable
            g.clear_filters()
            helper.remove_properties(g)
            g.vp.color[root] = old_color


def plot_statistics(filename):

    data = pd.read_csv(filename + '.csv')
    data.set_index('Unnamed: 0', inplace=True)
    # for percentage

    help = data.sum()
    found_not_predicted = (data.iloc[0] / help) * 100.
    found = (data.iloc[1]/ help) * 100.
    not_found_sat = (data.iloc[2] / help) * 100.
    not_found_mc = (data.iloc[3] / help) * 100.
    really_not_found = (data.iloc[4] / help) * 100.
    really_really_not_found = (data.iloc[5] / help) * 100.

    x = range(len(really_really_not_found))
    plt.figure()
    ax = plt.subplot(111)

    p0 = ax.bar(x, found_not_predicted, width=0.8)
    p1 = ax.bar(x, found, bottom=found_not_predicted, width=0.8)
    p2 = ax.bar(x, not_found_sat, bottom=[f + nf for f, nf in zip(found_not_predicted, found)], width=0.8)
    p3 = ax.bar(x, not_found_mc, bottom=[fnp + f + nf for f, nf, fnp in zip(found, not_found_sat, found_not_predicted)], width=0.8)
    p4 = ax.bar(x, really_not_found, bottom=[f + nf + rnf + fnp for f, nf, rnf, fnp in zip(found, not_found_sat, not_found_mc, found_not_predicted)], width=0.8)
    p5 = ax.bar(x, really_really_not_found, bottom=[f + nf + rnf + fnp + bla for f, nf, rnf, fnp, bla in zip(found, not_found_sat, not_found_mc, really_not_found, found_not_predicted)],  width=0.8)

    lgd = ax.legend((p0[0], p1[0], p2[0], p3[0], p4[0], p5[0]), data.index, loc="center left", bbox_to_anchor=(1, 0, 0.5, 1))
    plt.title('Targets in organisms')
    plt.xticks(x, data.columns, rotation=45, ha='right')
    plt.xlabel("model")
    plt.ylabel("percentage of targets")
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f %%'))

    plt.subplots_adjust(bottom=0.2)
    plt.savefig(filename + '.pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')


if __name__ == '__main__':

    statistics = OrganismResultStatistics()
    statistics.organism_statistics()
    statistics.neighbor_graphs()
    plot_statistics(join(const.MECAT_BASE, 'targets_in_organisms'))

