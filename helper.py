import graph_tool.all as gt


def f(value):
    try:
        val = float(value)
        return val
    except ValueError:
        return float('NaN')


def add_properties(g):
    g.vertex_properties["size"] = g.new_vertex_property("int", val=40)
    g.vertex_properties["dist"] = g.new_vertex_property("int", val=0)
    g.vertex_properties["besucht"] = g.new_vertex_property("bool", False)


def remove_properties(g):
    del g.vertex_properties["dist"]
    del g.vertex_properties["size"]
    del g.vertex_properties["besucht"]


class NeighborVisitor(gt.BFSVisitor):

    def __init__(self, g):
        self.g = g

    def discover_vertex(self, u):
        self.g.vp.besucht[u] = True

    def examine_vertex(self, u):
        pass

    def tree_edge(self, e):
        self.g.vp.dist[e.target()] = self.g.vp.dist[e.source()] + 1
