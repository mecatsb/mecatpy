# coding=utf-8
from __future__ import division, absolute_import, print_function
import pandas as pd
from pylab import *
from KEGG.modelbuilder import *
import csv
import graph_tool.all as gt
import os
import helper

from scipy import stats


class Statistics:

    def __init__(self, name, builder, path_to_model):
        self.builder = builder
        self.model_reactions = builder.model.reactions
        self.model_compounds = builder.model.metabolites
        self.rpairs = builder.rpairs
        self.arcgraph = builder.arcgraph
        self.start_compounds = [c.id for c in self.builder.start_compound_list]
        self.reaction_pathways = dict()
        self.bfs_reachable = list()
        self.bfs_not_reachable = list()
        self.targets = set([t.id for t in builder.all_target_compounds])
        self.name = name
        self.model_path = path_to_model

        self.statistics_path = join(path_to_model, "Statistik")
        if not os.path.exists(self.statistics_path):
            os.makedirs(self.statistics_path)

        with open(join(const.DATAPATH_KEGG, "KEGG_reaction.pickle"), 'rb') as handle:
            self.reactions = load(handle)

        with open(join(const.DATAPATH_KEGG, "KEGG_compound.pickle"), 'rb') as handle:
            self.compounds = load(handle)

        with open(join(const.DATAPATH_KEGG, "KEGG_pathway.pickle"), 'rb') as handle:
            self.pathways = cPickle.load(handle)

    def print_metabolite_statistics(self, filename):

        compounds = len(self.model_compounds)
        generic = len(self.builder.generic_compound_list)
        start = len(self.builder.start_compound_list)
        basis = len(self.builder.basis_compound_list)
        cofactors = len(self.builder.cofactor_list)
        pool = cofactors + basis + start
        external = compounds - pool

        with open(join(self.statistics_path, filename + self.builder.model.name + '.txt'), 'w') as f:
            f.write('metabolites M:' + '\t' + str(compounds) + '\n')
            f.write('external metabolites :' + '\t' + str(external) + '\n')
            f.write('\t' + 'generic metabolites:' + '\t' + str(generic) + '\n')
            f.write('metabolite pool Em:' + '\t' + str(pool) + '\n')
            f.write('\t' + 'start metabolites S:' + '\t' + str(start) + '\n')
            f.write('\t\t' + 'basis metabolites B:' + '\t' + str(basis) + '\n')
            f.write('\t' + 'cofactors:' + '\t' + str(cofactors) + '\n')

    def print_reaction_statistics(self, filename):
        reactions = len(self.model_reactions)
        unique_reactions = len(self.builder.unique_reactions)
        reversible_reactions = len(self.builder.reversible_reactions)

        with open(join(self.statistics_path, filename + self.builder.model.name + '.txt'), 'w') as f:
            f.write('reactions R:' + '\t' + str(reactions) + '\n')
            f.write('unique reactions: ' + '\t' + str(unique_reactions) + '\n')
            f.write('reversible reactions: ' + '\t' + str(reversible_reactions) + '\n')

    def draw_graph(self):
        g = self.arcgraph
        pos = gt.sfdp_layout(self.arcgraph)
        gt.graph_draw(g, pos, output_size=(10000, 10000), vertex_fill_color=g.vp.color, edge_pen_width=10, output=join('/mnt', 'g', 'LisaDaten',
                                                                                         'Paper2', 'figures' ,'arcgraph' + self.name + '.pdf'))

    def topology(self):
        g = self.arcgraph.copy()
        components, chist = gt.label_components(g, directed=False)  # directed = False because True would look for strongly connected components
        self.__plot_component_hist(chist, 'componenthist')
        start_components = set()
        number_compounds_in_start_components = 0
        for c in self.start_compounds:
            for v in gt.find_vertex(g,g.vp.compound_ids,c):
                start_components.add(components[v])

        cg = gt.Graph()
        cg.vertex_properties["size"] = cg.new_vertex_property("int", val=10)
        for c in start_components:
            v = cg.add_vertex()
            cg.vp.size[v] = chist[c]
            number_compounds_in_start_components += chist[c]

        satellites = set()

        clustering_coefficient = gt.global_clustering(g)
        with open(join(self.statistics_path, "clustering_coefficient.txt"), 'w') as f:
            f.write(str(clustering_coefficient[0]) + '\t' + str(clustering_coefficient[1]) + '\n')

        with open(join(self.statistics_path, "compounds_components.txt"), 'w') as f, \
                open(join(self.statistics_path, "component_hist.txt"), 'w') as f2:

            for componentid, elem in enumerate(chist):
                u = gt.GraphView(g, vfilt=components.a == componentid)
                u = gt.Graph(u, prune=True)

                f2.write(str(componentid+1) + '\t' + str(elem) + '\n')

                for v in u.vertices():
                    f.write(str(componentid+1) + '\t' + u.vp.compound_ids[v] + '\t' + u.vp.name[v] + '\n')

                    if componentid not in start_components:
                        satellites.add(u.vp.compound_ids[v])

                # gt.graph_draw(u, output=join(self.statistics_path, "component{i}.pdf".format(i=componentid)))

        targets_in_main_component = self.targets - satellites
        targets_in_satellites = self.targets & satellites

        with open(join(self.statistics_path, "targets_in_main_component.txt"), 'w') as f:
            for c in targets_in_main_component:
                compound = self.builder.compounds[c]
                f.write(c + '\t' + compound.names[0] + '\n')

        with open(join(self.statistics_path, "targets_in_satellites.txt"), 'w') as f:
            for c in targets_in_satellites:
                compound = self.builder.compounds[c]
                f.write(c + '\t' + compound.names[0] + '\n')

        with open(join(self.statistics_path, "components_with_start_metabolites.txt"), 'w') as f:
            for cid in start_components:
                f.write(str(cid) + '\n')

        p = number_compounds_in_start_components/g.num_vertices()*100

        with open(join(const.MECAT_BASE, "component_table.txt"), 'a') as f:
            f.write(self.name + ' & ' + str(len(chist)) + ' & ' + str(np.amax(chist)) + ' & ' + str(len(start_components)) + ' & ' + str(int(number_compounds_in_start_components)) + ' & ' + str(int(round(p, 0))) + '\%' + '\\\\ \n')

        #largest = gt.label_largest_component(g, directed=False)
        #gt.graph_draw(g, vertex_fill_color=largest, output=join(self.statistics_path,"largest_component.pdf"))

        g.vertex_properties["start_components"] = g.new_vertex_property("string", val='white')

        for v in g.vertices():
            if components[v] in start_components:
                g.vp.start_components[v] = 'red'
            else:
                g.vp.start_components[v] = 'blue'

        gt.graph_draw(g, vertex_fill_color=g.vp.start_components, output=join('/mnt', 'g', 'LisaDaten', 'Paper2', 'figures' ,'arcgraph' + self.name + '.pdf'))

    def degree_distribution(self, g, name):
        total_hist = gt.vertex_hist(g, "total", float_count=False)
        self.__plot_degree(total_hist, self.name + ' ' + name + ' ' + "totaldegdist.pdf", "total node degree")

        # in_hist = gt.vertex_hist(g, "in", float_count=False)
        # self.__plot_degree(in_hist, self.name + ' ' + name + ' ' + "indegdistloglog.pdf", "in node degree")
        #
        # out_hist = gt.vertex_hist(g, "out", float_count=False)
        # self.__plot_degree(out_hist, self.name + ' ' + name + ' ' + "outdegdistloglog.pdf", "out node degree")

        [atot, stdm] = gt.vertex_average(g, "total")
        stdtot = stdm * np.sqrt(g.num_vertices())
        return atot, stdtot, stdm

    def vertex_degrees(self, g, n, name ):
        alldegrees = [(v.in_degree()+v.out_degree(), v) for v in g.vertices()]
        alldegrees_ids = {g.vp.compound_ids[v] : v.in_degree()+v.out_degree() for v in g.vertices() }
        sortedall = sorted(alldegrees, reverse=True)
        self.__print_degree(sortedall, name + "degrees")

        s = pd.Series(alldegrees_ids)
        topn = pd.DataFrame(s.groupby(s, sort=False)).nlargest(n, 0).apply({0: lambda x: x, 1: lambda x: list(x.index)})
        hubs = set()

        with open(join(self.statistics_path, 'hubs' + self.name + ' ' + name + ' ' + '.txt'), 'w') as f, \
                open(join(self.model_path, 'hubs' + self.name + ' ' + name + ' ' + '.txt'), 'w') as f2:
            for row in topn.itertuples(index=False):
                for id in getattr(row, '_1'):
                    hubs.add(id)
                    met = self.compounds[id]
                    f.write(str(getattr(row, '_0')) + '\t' + id + '\t' + met.names[0] + '\n')
                    f2.write(id + '\t' + met.names[0] + '\n')

        return hubs

        # indegrees = [(v.in_degree(), v) for v in g.vertices()]
        # sortedin = sorted(indegrees, reverse=True)
        # self.__print_degree(sortedin, "indegrees")
        #
        # outdegrees = [(v.out_degree(), v) for v in g.vertices()]
        # sortedout = sorted(outdegrees, reverse=True)
        # self.__print_degree(sortedout, "outdgrees")

    def reactions_in_pathways(self):
        with open(join(self.statistics_path, "reactions_pathways.txt"), 'w') as f:
            for r in self.model_reactions:
                if '-' not in r.id:
                    pathways = r.pathways
                    pathway_names = [self.pathways[p].name for p in pathways]
                    pws = '\t'.join(pathway_names)
                    f.write(r.id + '\t' + pws + '\n')
                    self.reaction_pathways[r.id] = pathway_names

    def target_reachability(self):
        # part 1: arcs
        g = self.arcgraph.copy()
        g.vertex_properties["size"] = g.new_vertex_property("int", val=40)
        g.vertex_properties["dist"] = g.new_vertex_property("int", val=0)
        g.vertex_properties["besucht"] = g.new_vertex_property("bool", False)

        root = g.add_vertex()
        for v in g.vertices():
            if g.vp.is_start_compound[v]:
                g.add_edge(root, v)
            if g.vp.compound_ids[v] in self.targets:
                g.vp.color[v] = 'green'

        helper.add_properties(g)
        g.vp.size[root] = 100
        g.vp.color[root] = "red"
        g.vp.besucht[root] = True
        visitor = helper.NeighborVisitor(g)
        gt.bfs_search(g, root, visitor)

        pos = gt.radial_tree_layout(g, g.vertex_index[root], r=0.001)
        g.vp.besucht[root] = False
        g.set_vertex_filter(g.vp.besucht)
        g.vp.size.a = 50 * 0.9 ** g.vp.dist.a + 10

        bfs_reachable = list()

        for v in g.vertices():
            if g.vp.compound_ids[v] in self.targets:
                bfs_reachable.append(g.vp.compound_ids[v])

        bfs_not_reachable = self.targets - set(bfs_reachable)

        self.bfs_reachable = sorted(bfs_reachable)
        self.bfs_not_reachable = sorted(bfs_not_reachable)

        assert len(bfs_reachable) + len(bfs_not_reachable) == len(self.targets)

        with open(join(self.statistics_path, "bfs_reachable.txt"), 'w') as f:
            for id in self.bfs_reachable:
                c = self.compounds[id]
                f.write("{}\t{}\n".format(id, c.names[0]))

        with open(join(self.statistics_path, "bfs_not_reachable.txt"), 'w') as f:
            for id in self.bfs_not_reachable:
                c = self.compounds[id]
                f.write("{}\t{}\n".format(id, c.names[0]))

        g.set_vertex_filter(None)

        # part 2: reactions
        model_reactions = {r.id: r for r in self.model_reactions}
        feasible_reactions = self.builder.feasible_reaction_list
        feasible_targets = set()
        not_feasible_targets = set()

        predicted_targets = set()

        for c in self.targets:
            reactions_for_target = self.__reactions_for_target(c, model_reactions)
            fr = set(reactions_for_target) & feasible_reactions
            if len(fr) > 0:
                predicted_targets.add(c)
                if c in self.bfs_reachable:
                    feasible_targets.add(c)
            else:
                if c in self.bfs_reachable:
                    not_feasible_targets.add(c)

        with open(join(self.statistics_path, "target_is_feasible.txt"), 'w') as f:
            for id in sorted(feasible_targets):
                c = self.compounds[id]
                f.write("{}\t{}\n".format(id, c.names[0]))

        numt = len(self.targets)        # all targets
        numbfsf = len(bfs_reachable)    # BFS reachable
        numf = len(feasible_targets)    # feasible targets
        nump = len(predicted_targets)   # predicted targets

        p = numf/numt*100

        with open(join(const.MECAT_BASE, 'target_table.txt'), 'a') as f:
            f.write(self.name + ' & ' + str(numt) + ' & ' + str(numbfsf) + ' & ' + str(nump) + ' & ' + str(numf) +  ' & ' + str(int(round(p, 0))) + '\%'  + '\\\\ \n')

    def __reactions_for_target(self, compound_id, model_reactions):
        reactions = set()
        reaction_ids = self.compounds[compound_id].reactions

        for r in reaction_ids:
            forw = r
            if forw in model_reactions:
                reac = model_reactions[forw]
                if compound_id in reac.products():
                    reactions.add(forw)

            rev = '-'+r
            if rev in model_reactions:
                reac = model_reactions[rev]
                if compound_id in reac.products():
                    reactions.add(rev)
        return reactions

    def __print_degree(self, data, filename):
        g = self.arcgraph
        with open(join(self.statistics_path, filename + self.name + '.txt'), 'w') as f:
            for deg, v in data:
                f.write(str(deg) + '\t&\t' + g.vp.compound_ids[v] + '\t&\t' + g.vp.name[v] + '\\'+'\\' + '\n')

    def __plot_degree(self, data, filename, type):
        figure()
        bar(data[1][:-1], data[0])
        gca().set_yscale("symlog")

        subplots_adjust(left=0.2, bottom=0.2)
        xlabel("node degree")
        ylabel("log count")
        title(type + ' of model ' + self.name)

        tight_layout()
        savefig(join(self.statistics_path, filename))

    def __plot_component_hist(self, data, filename):
        plt.subplots()

        hist, bin_edges = np.histogram(data, bins = range(int(max(data))+1))
        non_zero = np.nonzero(hist)
        hist = hist[non_zero]
        bin_edges = bin_edges[non_zero]

        x_ticks = [str(int(edge)) for edge in bin_edges]
        indices = np.arange(len(bin_edges))

        plt.bar(indices, hist, 0.75, align='center')
        plt.xticks(indices, x_ticks)

        xlabel("component size")
        ylabel("number components")
        title('model ' + self.name)

        savefig(join('/mnt', 'g', 'LisaDaten', 'Paper2', 'figures', self.name + filename + '.pdf'))

    def model_properties(self):

        num_reactions_kegg = str(len(self.builder.all_reactions))
        num_reactions_reversible = str(len(self.builder.reversible_reactions))
        num_reactions_model = str(len(self.builder.model.reactions) - len(self.builder.reversible_reactions))
        num_feasible_reactions = str(len(self.builder.feasible_reaction_list))
        num_metabolites = str(len(self.builder.model.metabolites))
        num_metabolites_basis = str(len(self.builder.basis_compound_list))

        with open(join(const.MECAT_BASE, 'model_table.txt'), 'a') as f:
            f.write(self.name + '& &' + '{}/{}/{}&{}&{}/{}'.format(num_reactions_kegg, num_reactions_model,
                                                              num_reactions_reversible, num_feasible_reactions,
                                                              num_metabolites, num_metabolites_basis) + '\\\\ \n')


def topn_hubs(g, topn, name, f):
    num_vertices = g.num_vertices()
    num_arcs = g.num_edges()
    hubs = statistics.vertex_degrees(g, topn, name)
    (average_degree, std, stdmean) = statistics.degree_distribution(g, name)
    f.write(organism + ' & ' + str(num_vertices) + ' & ' + str(num_arcs) + ' & ' +
            str(round(average_degree, 2)) + ' (' + str(round(std, 2)) + ', ' + str(round(stdmean, 2)) + ') \\\\ \n')
    return hubs


if __name__ == '__main__':

    with open(join(const.MECAT_BASE, "organisms.txt"), 'r') as f:
        organisms = [x.strip() for x in f.readlines()]

    suffix = ['basis']

    hubs = []
    hubsfilter = []
    hubsall = []
    hubsallfilter = []

    with open(join(const.MECAT_BASE, "graphproperties_1_arcgraph.txt"), 'w') as f1, \
        open(join(const.MECAT_BASE, "graphproperties_2_arcgraphfilter.txt"), 'w') as f2, \
        open(join(const.MECAT_BASE, "graphproperties_4_arcgraphall.txt"), 'w') as f3, \
        open(join(const.MECAT_BASE, "graphproperties_3_arcgraphallfilter.txt"), 'w') as f4:
        f1.write('model\t&\tnumber of metabolites\t&number of arcs \\\\ \n')

        for organism in organisms:
            for s in suffix:
                path = join(const.MECAT_BASE, organism + s)
                with open(join(path, "KEGG_model.pickle"), 'rb') as handle:
                    builder = cPickle.load(handle)

                # builder.build_reaction_graph()
                # builder.save_arcgraph(const.MODEL)
                # builder.save_reactiongraph(const.MODEL)
                # builder.draw_graph(builder_kegg.reactiongraph, filename=join(const.DATAPATH_MODEL, "reactiongraph.pdf"))

                statistics = Statistics(organism, builder, path)
                statistics.print_metabolite_statistics('metabolites')
                statistics.print_reaction_statistics('reactions')
                #statistics.draw_graph()
                statistics.topology()
                #statistics.reactions_in_pathways()
                statistics.target_reachability()
                statistics.model_properties()

                # hubs in the network
                topn = 10

                hubs.append(topn_hubs(builder.arcgraph, topn, 'rclass_filter', f1)) # rclass no cofactors/inorganics, 1
                hubs_rclass = topn_hubs(builder.arcgraph_rclass_all, topn, 'rclass_all', f2) # rclass, 2
                hubsfilter.append(hubs_rclass)
                hubs_rpairs = topn_hubs(builder.arcgraph_rpairs_all, topn, 'rpairs_all', f3) # all substrate-product pairs, 4
                hubsall.append(hubs_rpairs)
                hubsallfilter.append(topn_hubs(builder.arcgraph_rpairs_filter, topn, 'rpairs_filter', f4)) # all substrate-product pairs, no cofactors/inorganics, 3

                cofactors = set([c.id for c in builder.cofactor_list])
                cofactors_in_topn_hubsrpair = cofactors & set(hubs_rpairs)
                cofactors_in_topn_hubsrclass = cofactors & set(hubs_rclass)

                with open(join(path, "Statistik", "cofactors_in_top" + str(topn) + ".txt"), 'w') as f:
                    f.write("rpairs:" + '\t' + str(len(cofactors_in_topn_hubsrpair)) + '\t' + str(len(hubs_rpairs)) + '\n' )
                    f.write("rclass:" + '\t' + str(len(cofactors_in_topn_hubsrclass)) + '\t' + str(len(hubs_rclass)) + '\n')

        common_hubs = set.intersection(*hubs)
        common_hubsfilter = set.intersection(*hubsfilter)
        common_hubsall = set.intersection(*hubsall)
        common_hubsallfilter = set.intersection(*hubsallfilter)

        with open(join(const.MECAT_BASE, "common_hubs.txt"), 'w') as f:
            f.write("rpairs without cofactors/inorganics: " +'\t'.join(common_hubs) + '\n')
            f.write("rpairs : " + '\t'.join(common_hubsfilter) + '\n')
            f.write("all substrate-product pairs: " + '\t'.join(common_hubsall) + '\n')
            f.write("all substrate-product pairs without cofactors/inorganics: "+ '\t'.join(common_hubsallfilter) + '\n')
