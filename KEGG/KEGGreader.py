from KEGG.reaction import Reaction
from KEGG.compound import Compound
from KEGG.rclass import RClass
from KEGG.enzyme import Enzyme
from KEGG.organism import Organism
from KEGG.pathway import Pathway
from os.path import join
import cPickle
import constants as const
import re
from scipy import io as sio


def read_kegg_reaction_dump(dump):
    f = open(dump)
    reactions = f.read()
    reactions = reactions.split('///')
    reactions = [Reaction(r) for r in reactions if r.strip() != '']

    reactiondict = dict()
    with open(join(const.DATAPATH_KEGG, "KEGG_reaction_ids.txt"), 'wb') as handle:
        for r in reactions:
            reactiondict[r.id] = r
            handle.write(r.id + '\t' + r.equation + '\n')

        __reactiondictToMATLABMap(reactions)

    with open(join(const.DATAPATH_KEGG, "KEGG_reaction.pickle"), 'wb') as handle:
        cPickle.dump(reactiondict, handle, protocol=cPickle.HIGHEST_PROTOCOL)


def read_kegg_compound_dump(dump):
    f = open(dump)
    compounds = f.read()
    compounds = compounds.split('///')
    compounds = [Compound(c) for c in compounds if c.strip() != '']

    compounddict = dict()
    has_brite = 0
    for c in compounds:
        compounddict[c.id] = c

        if c.brite:
            has_brite += 1

    __compounddictToMATLABMap(compounds)

    with open(join(const.DATAPATH_KEGG, "KEGG_compound.pickle"), 'wb') as handle:
        cPickle.dump(compounddict, handle, protocol=cPickle.HIGHEST_PROTOCOL)


def read_kegg_enzyme_dump(dump):
    f = open(dump)
    enzymes = f.read()
    enzymes = enzymes.split('///')
    enzymes = [Enzyme(e) for e in enzymes if e.strip() != '']

    enzymedict = dict()
    for e in enzymes:
        enzymedict[e.id] = e

    __enzymedictToMATLABMap(enzymes)

    with open(join(const.DATAPATH_KEGG, "KEGG_enzyme.pickle"), 'wb') as handle:
        cPickle.dump(enzymedict, handle, protocol=cPickle.HIGHEST_PROTOCOL)


def read_kegg_organism_dump(dump):

    with open(dump) as f:
        organisms = f.readlines()
    organisms = [x.strip() for x in organisms]
    organisms = [Organism(o) for o in organisms if o.strip() != '']

    organismdict = dict()
    for o in organisms:
        organismdict[o.name] = o

    with open(join(const.DATAPATH_KEGG, "KEGG_organism.pickle"), 'wb') as handle:
        cPickle.dump(organismdict, handle, protocol=cPickle.HIGHEST_PROTOCOL)


def read_kegg_pathway_dump(dump):

    with open(dump) as f:
        pathways = f.readlines()
    pathways = [x.strip() for x in pathways]
    pathways = [Pathway(p) for p in pathways if p.strip() != '']

    pathwaydict = dict()
    for p in pathways:
        pathwaydict[p.id] = p

    __pathwaydictToMATLABMap(pathways)

    with open(join(const.DATAPATH_KEGG, "KEGG_pathway.pickle"), 'wb') as handle:
        cPickle.dump(pathwaydict, handle, protocol=cPickle.HIGHEST_PROTOCOL)


def read_kegg_rclass_dump(dump):
    f = open(dump)
    rclasses = f.read()
    rclasses = rclasses.split('///')
    rcs = [RClass(rc) for rc in rclasses if rc.strip() != '']

    rclassdict = dict()
    for rc in rcs:
        rclassdict[rc.id] = rc

    with open(join(const.DATAPATH_KEGG, "KEGG_rclass.pickle"), 'wb') as handle:
        cPickle.dump(rclassdict, handle, protocol=cPickle.HIGHEST_PROTOCOL)


def read_kegg_pathway_hierarchy_dump(dump):

    f = open(dump)
    entry = f.read()

    pathway_hierarchy_dict = dict()
    pathway_dict = dict()

    for l in re.split('^A', entry, flags=re.MULTILINE)[1:]:
        A = re.search('>(.*)<', l).group(0).strip('><')
        category_dict = dict()
        for l2 in re.split('^B', l, flags=re.MULTILINE)[1:]:
            Bs = re.split('^C', l2, flags=re.MULTILINE)
            B = Bs[0].strip('\t\n\r ')
            sub_category_dict = dict()
            for l3 in Bs[1:]:
                C = re.findall('(\d{5})\s*(\S[^\n]*)', l3)
                sub_category_dict[C[0][0]] = C[0][1]
                pathway_dict[C[0][0]] = (A, B, C[0][1])
            category_dict[B] = sub_category_dict
        pathway_hierarchy_dict[A] = category_dict

    with open(join(const.DATAPATH_KEGG, "KEGG_pathway_hierachy.pickle"), 'wb') as handle:
        cPickle.dump(pathway_hierarchy_dict, handle, protocol=cPickle.HIGHEST_PROTOCOL)

    with open(join(const.DATAPATH_KEGG, "KEGG_pathways.pickle"), 'wb') as handle:
        cPickle.dump(pathway_dict, handle, protocol=cPickle.HIGHEST_PROTOCOL)


def read_kegg_atom_type_dump(dump):
    atom_type_dict = dict()
    with open(dump, 'r') as handle:
        next(handle)
        functional_group = ''
        for i, l in enumerate(handle):
            lp = l.strip().split(',')
            if lp[1] is not '':
                functional_group = lp[1]
            if lp[2] is not '':
                k = lp[2]
                atom_type_dict[k] = (functional_group, lp[3], lp[4])

    with open(join(const.DATAPATH_KEGG, "KEGG_atom_types.pickle"), 'wb') as handle:
        cPickle.dump(atom_type_dict, handle, protocol=cPickle.HIGHEST_PROTOCOL)


def __reactiondictToMATLABMap(reactions):
    reaction_dict = dict()

    for r in reactions:
        data = (r.id, r.names, r.equation, r.enzymes, r.reactants, r.rpairs)
        reaction_dict[r.id] = data

    sio.savemat(join(const.DATAPATH_KEGG, "KEGG_REACTION_dict"), {'reaction_dict': reaction_dict})


def __compounddictToMATLABMap(compounds):
    compound_dict = dict()

    for c in compounds:
        data = (c.names, c.formula, c.reactions, c.id, c.comment, c.molweight)
        compound_dict[c.id] = data

    sio.savemat(join(const.DATAPATH_KEGG, "KEGG_COMPOUND_dict"), {'compound_dict': compound_dict})


def __enzymedictToMATLABMap(enzymes):
    enzyme_dict = dict()

    for e in enzymes:
        organisms = list(e.organisms)
        data = (e.names, e.id, organisms)
        enzyme_dict[e.id] = data

    sio.savemat(join(const.DATAPATH_KEGG, "KEGG_ENZYME_dict"), {'enzyme_dict': enzyme_dict})


def __pathwaydictToMATLABMap(pathways):
    pathway_dict = dict()

    for p in pathways:
        data = (p.id, p.name)
        pathway_dict['p' + p.id] = data

    sio.savemat(join(const.DATAPATH_KEGG, "KEGG_PATHWAY_dict"), {'pathway_dict': pathway_dict})


if __name__ == '__main__':

    read_kegg_reaction_dump(join(const.DATAPATH_KEGG, "KEGG_reaction_raw.txt"))
    read_kegg_compound_dump(join(const.DATAPATH_KEGG, "KEGG_compound_raw.txt"))
    read_kegg_rclass_dump(join(const.DATAPATH_KEGG, "KEGG_rclass_raw.txt"))
    read_kegg_enzyme_dump(join(const.DATAPATH_KEGG, "KEGG_enzyme_raw.txt"))
    read_kegg_organism_dump(join(const.DATAPATH_KEGG, "KEGG_organism_raw.txt"))
    read_kegg_pathway_dump(join(const.DATAPATH_KEGG, "KEGG_pathway_raw.txt"))
