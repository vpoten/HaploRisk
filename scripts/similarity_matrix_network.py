import json
import os
import re
import itertools
import math

import networkx as nx
import matplotlib.pyplot as plt


def __dist_func(val):
    return max(1.0 - val, 1e-6)


def create_graph_from_go_similarity(sim_input_json_file, _label=None):
    data = json.load(open(sim_input_json_file))

    # get feature names and sort them
    names = set(map(lambda d: d['product1']['name'], data))
    names.update(set(map(lambda d: d['product2']['name'], data)))

    if _label is None:
        _label = __default_node_label

    names = map(lambda n: _label(n), names)
    names = sorted(names)

    G = nx.Graph()
    G.add_nodes_from(names)

    # add edges to the graph (distance < 1)
    edges = map(
        lambda d: (_label(d['product1']['name']), _label(d['product2']['name']), __dist_func(d['similarity'])),
        data)
    edges = filter(lambda t: t[2] < 1.0 and t[0] != t[1], edges)

    G.add_weighted_edges_from(edges)

    return G


def draw_graph(G):
    options = {
        'with_labels': True,
        'font_weight': 'bold',
        'font_color': 'grey',
        'node_size': 3000,
        'edge_color': range(len(G.edges)),
        'node_color': range(len(G.nodes)),
        'edge_cmap': plt.cm.Blues,
        'width': 2,
        'cmap': plt.cm.Blues
    }

    nx.draw(G, pos=nx.spring_layout(G), **options)
    plt.show()


def __default_node_label(name):
    return name


def __node_label_enrichr(name):
    """
    node label formatter for enrichr studies
    :param name:
    :return:
    """
    return name.split(' ')[0]


def print_connected_components(G):
    comp = list(nx.connected_components(G))
    isolated = filter(lambda cmp: len(cmp) == 1, comp)
    print "%i connected components" % len(comp)
    print "%i isolated components" % len(isolated)


def analyze_graph_connectivity(graphs):
    print ''

    for name in graphs:
        print '.... Analyzing connectivity %s' % name
        print_connected_components(graphs[name])

    isolateds = map(lambda name: sorted(nx.isolates(graphs[name])), graphs.keys())
    same_isolateds = all(map(lambda t: isolateds[t[0]] == isolateds[t[1]],
                             itertools.combinations(range(0, len(isolateds)), 2)))
    if same_isolateds is True:
        print '.... Analyzing connectivity\nAll similarity files produce the same set of isolated nodes'


def _not_isolated_nodes(G):
    return filter(lambda n: not nx.is_isolate(G, n), G.nodes)


def analyze_cliques(graphs):
    print ''

    for name in graphs:
        print '.... Analyzing cliques %s' % name
        cliques = list(nx.find_cliques(graphs[name]))
        print "%i maximal cliques found" % len(cliques)

        for i, clique in enumerate(cliques):
            print '[%i] %i nodes' % (i + 1, len(clique))


if __name__ == "__main__":
    base_path = '/home/victor/Escritorio/Genotipado_Alternativo/colocalizacion/go_distances'

    file_pattern = 'go_similarities_(\d+)_BP\.json'
    prog = re.compile(file_pattern)

    sim_files = [f for f in os.listdir(base_path) if os.path.isfile(os.path.join(base_path, f)) and prog.match(f)]

    graphs = {sim_file: create_graph_from_go_similarity(os.path.join(base_path, sim_file),
                                                        _label=__node_label_enrichr) for sim_file in sim_files}

    analyze_graph_connectivity(graphs)

    graphs2 = {name: graphs[name].subgraph(_not_isolated_nodes(graphs[name])).copy() for name in graphs.keys()}

    analyze_cliques(graphs2)

    draw_graph(graphs2['go_similarities_010_BP.json'])
