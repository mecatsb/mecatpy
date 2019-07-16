from KEGG.modelbuilder import Modelbuilder
from os import makedirs
from os.path import join
from os.path import exists
import constants as const
import cPickle

with open(join(const.MECAT_BASE, "organisms.txt"), 'r') as f:
    organisms = [x.strip() for x in f.readlines()]

suffix = ['basis']
mode = 'standard'  # 'standard': 'normal' filter, 'extended': no general reactions, 'reversible': every reaction is reversible

for organism in organisms:
    for s in suffix:
        organism_model_path = join(const.MECAT_BASE, organism+s)
        if not exists(organism_model_path):
            makedirs(organism_model_path)

        # build model with all reaction in host
        builder = Modelbuilder(organism+s, organism_model_path)

        if organism == 'kegg':
            builder.buildKEGG('kegg', mode)

        else:
            builder.build_host(organism, mode)

        builder.write_thermodynamics('thermodynamics')
        builder.write_reactions('reactions')
        builder.write_reversible_reactions('reversible_reactions')

        # metabolites
        builder.read_cofactors(join(const.DATAPATH_KEGG, 'cofactors_inorganic.txt'))
        builder.basis_compounds(join(const.DATAPATH_KEGG, 'possible_terminal_metabolites.txt'))

        builder.filter_compounds_weight(0, 300)
        builder.filter_compounds_carbons(0, 1)
        builder.stoichiometric_matrix('S')

        builder.build_metabolite_pool(s == 'basis')
        builder.write_metabolite_lists()

        # arcs for pathfinding (only rpairs and no cofactors/inorganics)
        if mode == 'extended':
            rpairs, reaction_rpair, rpair_reaction = builder.build_arcs('both', True, 'RxR_kegg')
        else:
            rpairs, reaction_rpair, rpair_reaction = builder.build_arcs('rclass', True, 'RxR_kegg')
        builder.rpairs = rpairs
        builder.reaction_rpair = reaction_rpair
        builder.rpair_reaction = rpair_reaction
        builder.arcgraph = builder.build_arc_graph(builder.rpairs)

        # rpairs
        rpairs, reaction_rpair, rpair_reaction = builder.build_arcs('rclass', False, 'RxR_kegg_filter')
        builder.arcgraph_rclass_all = builder.build_arc_graph(rpairs)

        # all arcs
        rpairs, reaction_rpair, rpair_reaction = builder.build_arcs('all', False, 'RxR_all')
        builder.arcgraph_rpairs_all = builder.build_arc_graph(rpairs)

        # all arcs, no cofactors/inorganics
        rpairs, reaction_rpair, rpair_reaction = builder.build_arcs('all', True, 'RxR_all_filter')
        builder.arcgraph_rpairs_filter = builder.build_arc_graph(rpairs)

        builder.feasible_reactions()

        builder.compounds_atom_types()
        builder.build_targets('targets_organism')

        with open(join(organism_model_path, "KEGG_model.pickle"), 'wb') as handle:
            cPickle.dump(builder, handle, protocol=cPickle.HIGHEST_PROTOCOL)

        builder.write_input_file()
