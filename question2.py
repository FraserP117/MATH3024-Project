import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import random
from itertools import chain

# import sys
#
# print(f"old recursion limit: {sys.getrecursionlimit()}")
# sys.setrecursionlimit(7000)
# print(f"new recursion limit: {sys.getrecursionlimit()}")



assortativities = []
probabilities = []
transitivities = []
degree_centralities = []
eigenvector_centralities = []
percolation_centralities = []
average_degrees = []
diameters = []

WS_assortativities = []
WS_edge_rewire_probabilities = []
WS_transitivities = []
WS_degree_centralities = []
WS_eigenvector_centralities = []
WS_percolation_centralities = []
WS_average_degrees = []
WS_diameters = []


'''
QUESTIONS

1. do the number of edges in the assortative swaps need to remain constant?
2. How to use the ruth-hurwitz criterion for Q1
3. don't need to plot agianst average degree for Q2 a and b?

'''

def average_degree(G):
    return G.number_of_edges() / 300

def average_degree_centrality(G):
    return sum(nx.degree_centrality(G).values())/len(nx.degree_centrality(G).values())

def average_eigenvector_centrality(G):
    return sum(nx.eigenvector_centrality(G, max_iter = 1000, tol = 1000).values())/len(nx.eigenvector_centrality(G, max_iter = 1000, tol = 1000).values())

def average_percolation_centrality(G):
    nx.set_node_attributes(G, 0.1, 'percolation')
    percolation_cent_dict = nx.percolation_centrality(
        G = G,
        attribute='percolation',
    )
    return sum(nx.percolation_centrality(G).values())/len(nx.percolation_centrality(G).values())

# compute all relevant metrics for 100 random graphs with harmoically
# decreasing edge connection probabilities
# def generate_ER_graph_metrics():
#     for i in range(1, 101):
#         G = nx.erdos_renyi_graph(300, 1/i)
#         probabilities.append(1/i)
#         average_degrees.append(round(average_degree(G)))
#         assortativities.append(nx.degree_assortativity_coefficient(G))
#         transitivities.append(nx.transitivity(G))
#         degree_centralities.append(average_degree_centrality(G))
#         eigenvector_centralities.append(average_eigenvector_centrality(G))
#         percolation_centralities.append(average_percolation_centrality(G))
#
#         if nx.is_connected(G): # lcc
#             diameters.append(nx.diameter(G))
#         else:
#             diameters.append(np.Infinity)

def generate_ER_graph_metrics():
    for i in range(1, 101):
        G = nx.erdos_renyi_graph(300, 1/i)

        Gcc = sorted(nx.connected_components(G), key = len, reverse = True)
        G0 = G.subgraph(Gcc[0])

        probabilities.append(1/i)
        average_degrees.append(round(average_degree(G0)))
        assortativities.append(nx.degree_assortativity_coefficient(G0))
        transitivities.append(nx.transitivity(G0))
        degree_centralities.append(average_degree_centrality(G0))
        eigenvector_centralities.append(average_eigenvector_centrality(G0))
        percolation_centralities.append(average_percolation_centrality(G0))
        diameters.append(nx.diameter(G0))

# compute all relevant metrics for 100 random graphs with harmoically
# decreasing edge connection probabilities
def generate_WS_graph_metrics():
    for i in range(1, 101):

        G = nx.watts_strogatz_graph(300, 4, 1/i)

        Gcc = sorted(nx.connected_components(G), key = len, reverse = True)
        G0 = G.subgraph(Gcc[0])

        WS_edge_rewire_probabilities.append(1/i)
        WS_average_degrees.append(round(average_degree(G0)))
        WS_assortativities.append(nx.degree_assortativity_coefficient(G0))
        WS_transitivities.append(nx.transitivity(G0))
        WS_degree_centralities.append(average_degree_centrality(G0))
        WS_eigenvector_centralities.append(average_eigenvector_centrality(G0))
        WS_percolation_centralities.append(average_percolation_centrality(G0))

def assortative_rewiring_random(G):
    changed_assortativities = []
    number_swaps = []
    G_prime = G.copy()
    iter = 0
    for i in range(15000):
        A = nx.degree_assortativity_coefficient(G)
        G_prime = G.copy()
        nx.double_edge_swap(G_prime)
        if nx.degree_assortativity_coefficient(G_prime) > nx.degree_assortativity_coefficient(G):
            G = G_prime
            changed_assortativities.append(nx.degree_assortativity_coefficient(G_prime))
            number_swaps.append(iter)
            iter += 1

    plt.scatter(number_swaps, changed_assortativities, s = 1.5, color = 'darkviolet')
    plt.title("BA Graph Assortativity vs number of 2-OPT swaps")
    plt.xlabel("number 2-OPT swaps")
    plt.ylabel("assortativity")
    plt.show()

    return G

def disassortative_rewiring_random(G):
    changed_assortativities = []
    number_swaps = []
    G_prime = G.copy()
    iter = 0
    for i in range(15000):
        A = nx.degree_assortativity_coefficient(G)
        G_prime = G.copy()
        nx.double_edge_swap(G_prime)
        if nx.degree_assortativity_coefficient(G_prime) < nx.degree_assortativity_coefficient(G):
            G = G_prime
            changed_assortativities.append(nx.degree_assortativity_coefficient(G_prime))
            number_swaps.append(iter)
            iter += 1

    plt.scatter(number_swaps, changed_assortativities, s = 1.5, color = 'darkviolet')
    plt.title("BA Graph Assortativity vs number of 2-OPT swaps")
    plt.xlabel("number 2-OPT swaps")
    plt.ylabel("assortativity")
    plt.show()

    return G


'''
If new_edge has any vertices in common with edge_list, return true
else return False.
'''
def any_common_ele(edge_list, new_edge):
    full = [item for t in edge_list for item in t]
    for i in new_edge:
        if i in full:
            return True
    return False

'''
select 2 edges from G_prime: (x, y) (u, v)
find the degrees of x, y, u, v
delete the edges (x, y), (u, v)
connect the highest degree vertex with the second hightes degree vertex
and connect the lowest degree vertex with the second lowest degree vertex

must check if the selected edges are unique
must check if proposed edges don't already exist

selected edges: (u, v), (w, x)
d_0 = |k_u - k_v| + |k_w - k_x|
d_1 = |k_u - k_w| + |k_v - k_x|
d_2 = |k_u - k_x| + |k_v - k_w|

if d_0 > min(d_1, d_2):
    take the min
    add edge for max assortativity

if d_0 < max(d_1, d_2):
    take the max
    add edge for min assortativity
'''
def assortative_rewiring(G):
    changed_assortativities = []
    number_swaps = []
    iter = 0
    for i in range(15000):

        G_prime = G.copy()

        # randomly select 2 edges and remove these from the graph
        all_edges = [e for e in G_prime.edges]
        selected_edges = []
        while len(selected_edges) != 2:

            random_edge = random.choice(all_edges)

            if not any_common_ele(selected_edges, random_edge):
                selected_edges.append(random_edge)

            print(f"selected_edges: {selected_edges}")

        d_0 = np.abs(G_prime.degree[selected_edges[0][0]] - G_prime.degree[selected_edges[0][1]]) + np.abs(G_prime.degree[selected_edges[1][0]] - G_prime.degree[selected_edges[1][1]]) # d_0 = |k_u - k_v| + |k_w - k_x|
        d_1 = np.abs(G_prime.degree[selected_edges[0][0]] - G_prime.degree[selected_edges[1][0]]) + np.abs(G_prime.degree[selected_edges[0][1]] - G_prime.degree[selected_edges[1][1]]) # d_1 = |k_u - k_w| + |k_v - k_x|
        d_2 = np.abs(G_prime.degree[selected_edges[0][0]] - G_prime.degree[selected_edges[1][1]]) + np.abs(G_prime.degree[selected_edges[0][1]] - G_prime.degree[selected_edges[1][0]]) # d_2 = |k_u - k_x| + |k_v - k_w|
        min_d = min(d_1, d_2)

        if d_0 > min_d:
            if min_d == d_1:
                # add edges between selected edges for maximal assortativity:
                if G_prime.has_edge(selected_edges[0][0], selected_edges[1][0]):
                    pass
                else:
                    print(f"edge added: {(selected_edges[0][0], selected_edges[1][0])}")
                    G_prime.add_edge(selected_edges[0][0], selected_edges[1][0])

                if G_prime.has_edge(selected_edges[0][1], selected_edges[1][1]):
                    pass
                else:
                    print(f"edge added: {(selected_edges[0][1], selected_edges[1][1])}")
                    G_prime.add_edge(selected_edges[0][1], selected_edges[1][1])

                # remove the selected edges from the graph:
                for edge in selected_edges:
                    G_prime.remove_edge(*edge)

            elif min_d == d_2:
                if G_prime.has_edge(selected_edges[0][0], selected_edges[1][1]):
                    pass
                else:
                    # add edges between selected edges for maximal assortativity:
                    print(f"edge added: {(selected_edges[0][0], selected_edges[1][1])}")
                    G_prime.add_edge(selected_edges[0][0], selected_edges[1][1])

                if G_prime.has_edge(selected_edges[0][1], selected_edges[1][0]):
                    pass
                else:
                    print(f"edge added: {(selected_edges[0][1], selected_edges[1][0])}")
                    G_prime.add_edge(selected_edges[0][1], selected_edges[1][0])

                # remove the selected edges from the graph:
                for edge in selected_edges:
                    G_prime.remove_edge(*edge)

        G = G_prime

        # increment the iterations:
        iter += 1

    return G


def disassortative_rewiring(G):
    changed_assortativities = []
    number_swaps = []
    iter = 0
    for i in range(15000):

        G_prime = G.copy()

        # randomly select 2 edges and remove these from the graph
        all_edges = [e for e in G_prime.edges]
        selected_edges = []
        while len(selected_edges) != 2:

            random_edge = random.choice(all_edges)

            if not any_common_ele(selected_edges, random_edge):
                selected_edges.append(random_edge)

            print(f"selected_edges: {selected_edges}")

        d_0 = np.abs(G_prime.degree[selected_edges[0][0]] - G_prime.degree[selected_edges[0][1]]) + np.abs(G_prime.degree[selected_edges[1][0]] - G_prime.degree[selected_edges[1][1]]) # d_0 = |k_u - k_v| + |k_w - k_x|
        d_1 = np.abs(G_prime.degree[selected_edges[0][0]] - G_prime.degree[selected_edges[1][0]]) + np.abs(G_prime.degree[selected_edges[0][1]] - G_prime.degree[selected_edges[1][1]]) # d_1 = |k_u - k_w| + |k_v - k_x|
        d_2 = np.abs(G_prime.degree[selected_edges[0][0]] - G_prime.degree[selected_edges[1][1]]) + np.abs(G_prime.degree[selected_edges[0][1]] - G_prime.degree[selected_edges[1][0]]) # d_2 = |k_u - k_x| + |k_v - k_w|
        max_d = max(d_1, d_2)

        if d_0 < max_d:
            if max_d == d_1:
                # add edges between selected edges for maximal assortativity:
                if G_prime.has_edge(selected_edges[0][0], selected_edges[1][0]):
                    pass
                else:
                    print(f"edge added: {(selected_edges[0][0], selected_edges[1][0])}")
                    G_prime.add_edge(selected_edges[0][0], selected_edges[1][0])

                if G_prime.has_edge(selected_edges[0][1], selected_edges[1][1]):
                    pass
                else:
                    print(f"edge added: {(selected_edges[0][1], selected_edges[1][1])}")
                    G_prime.add_edge(selected_edges[0][1], selected_edges[1][1])

                # remove the selected edges from the graph:
                for edge in selected_edges:
                    G_prime.remove_edge(*edge)

            elif max_d == d_2:
                if G_prime.has_edge(selected_edges[0][0], selected_edges[1][1]):
                    pass
                else:
                    # add edges between selected edges for maximal assortativity:
                    print(f"edge added: {(selected_edges[0][0], selected_edges[1][1])}")
                    G_prime.add_edge(selected_edges[0][0], selected_edges[1][1])

                if G_prime.has_edge(selected_edges[0][1], selected_edges[1][0]):
                    pass
                else:
                    print(f"edge added: {(selected_edges[0][1], selected_edges[1][0])}")
                    G_prime.add_edge(selected_edges[0][1], selected_edges[1][0])

                # remove the selected edges from the graph:
                for edge in selected_edges:
                    G_prime.remove_edge(*edge)

        G = G_prime

        # increment the iterations:
        iter += 1

    return G

# # compute all relevant metrics for 100 random graphs with harmoically
# # decreasing edge connection probabilities
# def generate_BA_graph_metrics():
#     for i in range(1, 10):
#         G = nx.barabasi_albert_graph(300, 6)
#         WS_edge_rewire_probabilities.append(1/i)
#         WS_average_degrees.append(round(average_degree(G)))
#         WS_assortativities.append(nx.degree_assortativity_coefficient(G))
#         WS_transitivities.append(nx.transitivity(G))
#         WS_degree_centralities.append(average_degree_centrality(G))
#         WS_eigenvector_centralities.append(average_eigenvector_centrality(G))
#         WS_percolation_centralities.append(average_percolation_centrality(G))
#
#         if nx.is_connected(G):
#             diameters.append(nx.diameter(G))
#         else:
#             diameters.append(np.Infinity)

# ''' Assortativity '''
# def show_assortativities(probabilities, average_degrees, assortativities):
#     # plot assortativity against p
#     plt.scatter(probabilities, assortativities, s = 2, color = 'blue')
#     plt.title("assortativity against p")
#     plt.xlabel("probability of edge connection")
#     plt.ylabel("assortativity")
#     plt.show()
#
#     # plot assortativity against average degree
#     plt.scatter(average_degrees, assortativities, s = 2, color = 'blue')
#     plt.title("assortativity against average degree")
#     plt.xlabel("probability of edge connection")
#     plt.ylabel("average degree")
#     plt.show()

def create_sandpile(num_grains):
    sandpile = {}
    # create the sandpile:
    # for node in range(0, 401):
    for node in range(0, num_grains):
        sandpile[node] = random.randint(0, 3)

    return sandpile

def is_saturated(sandpile, threshold):
    for node, num_grains in sandpile.items():
        if sandpile[node] < threshold:
            return False

    return True

def add_sand(num_grains_to_add, sandpile):
    for node, num_grains in sandpile.items():
        if num_grains == 3:
            # num_grains += num_grains_to_add
            sandpile[node] = num_grains + num_grains_to_add
            print(f"{num_grains} grains on node: {node}")
            return sandpile
    return sandpile
'''
return True iff at least 1 critical node in the lattice.
node is critical if its number of grains > threshold and at least
one of its neighbors has less than threshold grains.
'''
def any_critical(sandpile, threshold):
    num_critical = 0
    for node, num_grains in sandpile.items():
        if num_grains >= threshold:
        # if num_grains > threshold:
            neighbors = [n for n in lattice.neighbors(node)]
            for neighbor in neighbors:
                if sandpile[neighbor] < threshold:
                    num_critical += 1
    if num_critical == 0:
        return False
    else:
        return True

def is_critical(node, sandpile, threshold):
    if sandpile[node] >= threshold:
        neighbors = [n for n in lattice.neighbors(node)]
        for neighbor in neighbors:
            if sandpile[neighbor] < threshold:
                return True
    return False

'''
Avalanche size = the number of toppling events required to return the sandpile
to a restful state
'''
def sandpile_sim(threshold, lattice, sandpile):
    num_avalanches = 0
    while any_critical(sandpile, threshold):
        for node, num_grains in sandpile.items():
            if is_critical(node, sandpile, threshold):
                num_avalanches += 1
                sandpile[node] -= 4
                print(f"\n\ncritical node: {node}, num_grains: {num_grains}")
                neighbors = [n for n in lattice.neighbors(node)]
                print(f"neighbors: {neighbors}\n\n")
                for neighbor in neighbors:
                    sandpile[neighbor] += 1
    # print(f"sandpile: {sandpile.values()}")
    return num_avalanches


''' Assortativity '''
def show_assortativities(probabilities, average_degrees, assortativities, is_ER):
    # plot assortativity against p
    plt.scatter(probabilities, assortativities, s = 2, color = 'blue')
    plt.title("assortativity against p")
    plt.xlabel("probability of edge connection")
    plt.ylabel("assortativity")
    plt.show()

    if is_ER:
        # plot assortativity against average degree
        plt.scatter(average_degrees, assortativities, s = 2, color = 'blue')
        plt.title("assortativity against average degree")
        plt.xlabel("average degree")
        plt.ylabel("assortativity")
        plt.show()

''' Transitivity '''
def show_transitivities(probabilities, average_degrees, transitivities, is_ER):
    # plot transitivity against p
    plt.scatter(probabilities, transitivities, s = 2, color = 'red')
    plt.title("transitivity against p")
    plt.xlabel("probability of edge connection")
    plt.ylabel("transitivity")
    plt.show()

    if is_ER:
        # plot transitivity against p
        plt.scatter(average_degrees, transitivities, s = 2, color = 'red')
        plt.title("transitivity against average degree")
        plt.xlabel("average degree")
        plt.ylabel("transitivity")
        plt.show()

''' Degree Centrality '''
def show_degree_centralities(probabilities, average_degrees, degree_centralities, is_ER):
    # plot degree_centrality against average degree
    plt.scatter(probabilities, degree_centralities, s = 2, color = 'green')
    plt.title("degree centrality against p")
    plt.xlabel("probability of edge connection")
    plt.ylabel("degree centrality")
    plt.show()

    if is_ER:
        # plot degree_centrality against average degree
        plt.scatter(average_degrees, degree_centralities, s = 2, color = 'green')
        plt.title("degree centrality against average degree")
        plt.xlabel("average degree")
        plt.ylabel("degree centrality")
        plt.show()

''' Eigenvector Centrality '''
def show_eigenvector_centralities(probabilities, average_degrees, eigenvector_centralities, is_ER):
    # plot eigenvector_centralities against p
    plt.scatter(probabilities, eigenvector_centralities, s = 2, color = 'purple')
    plt.title("eigenvector centrality against p")
    plt.xlabel("probability of edge connection")
    plt.ylabel("eigenvector centrality")
    plt.show()

    if is_ER:
        # plot eigenvector_centralities against average degree
        plt.scatter(average_degrees, eigenvector_centralities, s = 2, color = 'purple')
        plt.title("eigenvector centrality against average degree")
        plt.xlabel("average degree")
        plt.ylabel("eigenvector centrality")
        plt.show()

''' percolation Centrality '''
def show_percolation_centralities(probabilities, average_degrees, percolation_centralities, is_ER):
    # plot eigenvector_centralities against p
    plt.scatter(probabilities, percolation_centralities, s = 2, color = 'orange')
    plt.title("percolation centrality against p")
    plt.xlabel("probability of edge connection")
    plt.ylabel("percolation centrality")
    plt.show()

    if is_ER:
        # plot percolation_centralities against average degree
        plt.scatter(average_degrees, percolation_centralities, s = 2, color = 'orange')
        plt.title("percolation centrality against average degree")
        plt.xlabel("average degree")
        plt.ylabel("percolation centrality")
        plt.show()

''' Diameter'''
def show_diameters(probabilities, average_degrees, diameters, is_ER):
    # plot diameters against p
    plt.scatter(probabilities, diameters, s = 1.5, color = 'lime')
    plt.title("graph diameter against p")
    plt.xlabel("probability of edge connection")
    plt.ylabel("diameter")
    plt.show()

    if is_ER:
        # plot diameters against average degree
        plt.scatter(average_degrees, diameters, s = 1.5, color = 'lime')
        plt.title("graph diameter against average degree")
        plt.xlabel("average degree")
        plt.ylabel("diameter")
        plt.show()


if __name__ == '__main__':

    # # question 2.a
    # generate_ER_graph_metrics()
    # show_degree_centralities(probabilities, average_degrees, degree_centralities, True)
    # show_assortativities(probabilities, average_degrees, assortativities, True)
    # show_transitivities(probabilities, average_degrees, transitivities, True)
    # show_eigenvector_centralities(probabilities, average_degrees, eigenvector_centralities, True)
    # show_percolation_centralities(probabilities, average_degrees, percolation_centralities, True)
    # show_diameters(probabilities, average_degrees, diameters, True)

    # # question 2.b
    # generate_WS_graph_metrics()
    # show_degree_centralities(WS_edge_rewire_probabilities, WS_average_degrees, WS_degree_centralities, False)
    # show_assortativities(WS_edge_rewire_probabilities, WS_average_degrees, WS_assortativities, False)
    # show_transitivities(WS_edge_rewire_probabilities, WS_average_degrees, WS_transitivities, False)
    # show_eigenvector_centralities(WS_edge_rewire_probabilities, WS_average_degrees, WS_eigenvector_centralities, False)
    # show_percolation_centralities(WS_edge_rewire_probabilities, WS_average_degrees, WS_percolation_centralities, False)
    # show_diameters(WS_edge_rewire_probabilities, WS_average_degrees, WS_diameters, False)

    # BA Graph G_0
    # G_0 = nx.barabasi_albert_graph(300, 6)
    G_0 = nx.barabasi_albert_graph(300, 6)

    # G_a = assortative_rewiring(G_0)
    G_a = disassortative_rewiring(G_0)

    # assortativity of G_0:
    A_1 = nx.degree_assortativity_coefficient(G_0)

    # Assortativity of G_a:
    A_2 = nx.degree_assortativity_coefficient(G_a)


    # # The comparatively assortative graph:
    # G_a = assortative_rewiring_random(G_0)
    #
    # # The comparatively disassortative graph:
    # G_b = disassortative_rewiring_random(G_0)
    #
    # # Assortativity of G_a:
    # A_2 = nx.degree_assortativity_coefficient(G_a)
    #

    nx.draw(G_0, node_size = 1.0, width = 0.5, node_color = "cyan", edge_color = "darkviolet")
    plt.show()
    nx.draw(G_a, node_size = 1.0, width = 0.5, node_color = "cyan", edge_color = "darkviolet")
    plt.show()
    #
    # print(f"\n\nA_1: {A_1}")
    # print(f"A_2 {A_2}\n\n")

    # lattice = nx.random_regular_graph(4, 400)
    # nx.draw(lattice, node_size = 2.0, width = 0.75, node_color = "lime", edge_color = "green")
    # plt.show()

    # threshold = 4
    # lattice = nx.random_regular_graph(4, 400)
    # sandpile = create_sandpile()
    # critical_sandpile = add_sand(2, sandpile)
    # num_av = sandpile_sim(4, lattice, sandpile)

    '''
    threshold = 4
    num_grains = 400
    lattice = nx.random_regular_graph(4, num_grains)
    sandpile = create_sandpile(num_grains)
    critical_sandpile = add_sand(20, sandpile)
    num_av = sandpile_sim(4, lattice, sandpile)

    print(f"\n\nnumber of avalanches: {num_av}")
    '''


    # nx.draw(lattice, node_size = 2.0, width = 0.75, node_color = "lime", edge_color = "green")
    # plt.show()
