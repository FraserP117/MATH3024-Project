import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import random


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

def assortative_swaps_random(G):
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

# def assortative_swaps(G):
#     changed_assortativities = []
#     number_swaps = []
#     G_prime = G.copy()
#     iter = 0
#     for i in range(15000):
#         A = nx.degree_assortativity_coefficient(G)
#         G_prime = G.copy()
#         '''
#         select 2 edges from G_prime: (x, y) (u, v)
#         find the degrees of x, y, u, v
#         delete the edges (x, y), (u, v)
#         connect the highest degree vertex with the second hightes degree vertex
#         and connect the lowest degree vertex with the second lowest degree vertex
#         '''
#
#         # randomly select 2 edges and remove these from the graph
#         all_edges = [e for e in G.edges]
#         edges_to_swap = []
#         while len(edges_to_swap) < 2:
#             random_edge = random.choice(all_edges)
#             if not random_edge in edges_to_swap:
#                 edges_to_swap.append(random_edge)
#
#         # create a dictionary where key = edge and val = degree:
#         vertex_degrees = {}
#         for edge in edges_to_swap:
#             for vertex in edge:
#                 vertex_degrees[vertex] = G.degree[vertex]
#
#         # sort the dictionary of vertex degrees by degree - max frst
#         vertex_degrees = sorted(vertex_degrees.items(), key = lambda x: x[1], reverse=True)
#         print(vertex_degrees)
#
#         # remove the selected edges from the graph:
#         for edge in edges_to_swap:
#             G.remove_edge(edge)
#
#         # add edges between selected edges for maximal assortativity:
#         for vertex, degree in vertex_degrees.items():
#
#     # plt.scatter(number_swaps, changed_assortativities, s = 1.5, color = 'darkviolet')
#     # plt.title("BA Graph Assortativity vs number of 2-OPT swaps")
#     # plt.xlabel("number 2-OPT swaps")
#     # plt.ylabel("assortativity")
#     # plt.show()
#
#     return G

# '''
# select 2 edges from G_prime: (x, y) (u, v)
# find the degrees of x, y, u, v
# delete the edges (x, y), (u, v)
# connect the highest degree vertex with the second hightes degree vertex
# and connect the lowest degree vertex with the second lowest degree vertex
# '''
# def assortative_swaps(G):
#     changed_assortativities = []
#     number_swaps = []
#     iter = 0
#     for i in range(150):
#
#         # randomly select 2 edges and remove these from the graph
#         all_edges = [e for e in G.edges]
#         edges_to_swap = []
#         while len(edges_to_swap) != 2:
#             random_edge = random.choice(all_edges)
#             if not random_edge in edges_to_swap:
#                 edges_to_swap.append(random_edge)
#
#         print(f"\n\nedges_to_swap: {edges_to_swap}")
#
#         # create a dictionary where key = edge and val = degree:
#         vertex_degrees = {}
#         for edge in edges_to_swap:
#             for vertex in edge:
#                 vertex_degrees[vertex] = G.degree[vertex]
#
#         # sort the dictionary of vertex degrees by degree - max frst
#         vertex_degrees = sorted(vertex_degrees.items(), key = lambda x: x[1], reverse = True)
#         print(f"vertex_degrees: {vertex_degrees}")
#
#         # add edges between selected edges for maximal assortativity:
#         G.add_edge(vertex_degrees[0][0], vertex_degrees[1][0])
#         print(f"edge added: {(vertex_degrees[0][0], vertex_degrees[1][0])}")
#         G.add_edge(vertex_degrees[2][0], vertex_degrees[3][0])
#         print(f"edge added: {vertex_degrees[2][0], vertex_degrees[3][0]}\n\n")
#
#         # remove the selected edges from the graph:
#         for edge in edges_to_swap:
#             G.remove_edge(*edge)
#
#         # increment the iterations:
#         iter += 1
#
#     plt.scatter(number_swaps, changed_assortativities, s = 1.5, color = 'darkviolet')
#     plt.title("BA Graph Assortativity vs number of 2-OPT swaps")
#     plt.xlabel("number 2-OPT swaps")
#     plt.ylabel("assortativity")
#     plt.show()
#
#     return G

# '''
# select 2 edges from G_prime: (x, y) (u, v)
# find the degrees of x, y, u, v
# delete the edges (x, y), (u, v)
# connect the highest degree vertex with the second hightes degree vertex
# and connect the lowest degree vertex with the second lowest degree vertex
# '''
# def assortative_swaps(G):
#     changed_assortativities = []
#     number_swaps = []
#     iter = 0
#     for i in range(150):
#
#         # randomly select 2 edges and remove these from the graph
#         all_edges = [e for e in G.edges]
#         edges_to_swap = []
#         while len(edges_to_swap) != 2:
#             random_edge = random.choice(all_edges)
#             if not random_edge in edges_to_swap:
#                 edges_to_swap.append(random_edge)
#
#         print(f"\n\nedges_to_swap: {edges_to_swap}")
#
#         # create a dictionary where key = edge and val = degree:
#         vertex_degrees = {} # use a list instead?? possible key overlap????????????????????????????????????????????????????????????????
#         for edge in edges_to_swap:
#             for vertex in edge:
#                 vertex_degrees[vertex] = G.degree[vertex]
#
#         # sort the dictionary of vertex degrees by degree - max frst
#         vertex_degrees = sorted(vertex_degrees.items(), key = lambda x: x[1], reverse = True)
#         print(f"vertex_degrees: {vertex_degrees}")
#
#         # add edges between selected edges for maximal assortativity:
#         G.add_edge(vertex_degrees[0][0], vertex_degrees[1][0])
#         print(f"edge added: {(vertex_degrees[0][0], vertex_degrees[1][0])}")
#         G.add_edge(vertex_degrees[2][0], vertex_degrees[3][0])
#         print(f"edge added: {vertex_degrees[2][0], vertex_degrees[3][0]}\n\n")
#
#         # remove the selected edges from the graph:
#         for edge in edges_to_swap:
#             G.remove_edge(*edge)
#
#         # increment the iterations:
#         iter += 1
#
#     plt.scatter(number_swaps, changed_assortativities, s = 1.5, color = 'darkviolet')
#     plt.title("BA Graph Assortativity vs number of 2-OPT swaps")
#     plt.xlabel("number 2-OPT swaps")
#     plt.ylabel("assortativity")
#     plt.show()
#
#     return G

# '''
# select 2 edges from G_prime: (x, y) (u, v)
# find the degrees of x, y, u, v
# delete the edges (x, y), (u, v)
# connect the highest degree vertex with the second hightes degree vertex
# and connect the lowest degree vertex with the second lowest degree vertex
# '''
# def assortative_swaps(G):
#     changed_assortativities = []
#     number_swaps = []
#     iter = 0
#     for i in range(1500):
#
#         # randomly select 2 edges and remove these from the graph
#         all_edges = [e for e in G.edges]
#         edges_to_swap = []
#         while len(edges_to_swap) != 2:
#             random_edge = random.choice(all_edges)
#             if not random_edge in edges_to_swap:
#                 edges_to_swap.append(random_edge)
#
#         print(f"\n\nedges_to_swap: {edges_to_swap}")
#
#         # create a dictionary where key = edge and val = degree:
#         vertex_degrees = [] # use a list instead?? possible key overlap????????????????????????????????????????????????????????????????
#         for edge in edges_to_swap:
#             for vertex in edge:
#                 vertex_degrees.append((vertex, G.degree[vertex]))
#
#         # sort the dictionary of vertex degrees by degree - max frst
#         vertex_degrees = sorted(vertex_degrees, key = lambda x: x[1], reverse = True)
#         print(f"vertex_degrees: {vertex_degrees}")
#
#         # add edges between selected edges for maximal assortativity:
#         G.add_edge(vertex_degrees[0][0], vertex_degrees[1][0])
#         print(f"edge added: {(vertex_degrees[0][0], vertex_degrees[1][0])}")
#         G.add_edge(vertex_degrees[2][0], vertex_degrees[3][0])
#         print(f"edge added: {vertex_degrees[2][0], vertex_degrees[3][0]}\n\n")
#
#         # remove the selected edges from the graph:
#         for edge in edges_to_swap:
#             G.remove_edge(*edge)
#
#         # increment the iterations:
#         iter += 1
#
#     # plt.scatter(number_swaps, changed_assortativities, s = 1.5, color = 'darkviolet')
#     # plt.title("BA Graph Assortativity vs number of 2-OPT swaps")
#     # plt.xlabel("number 2-OPT swaps")
#     # plt.ylabel("assortativity")
#     # plt.show()
#
#     return G


# '''
# select 2 edges from G_prime: (x, y) (u, v)
# find the degrees of x, y, u, v
# delete the edges (x, y), (u, v)
# connect the highest degree vertex with the second hightes degree vertex
# and connect the lowest degree vertex with the second lowest degree vertex
# '''
# def assortative_swaps(G):
#     changed_assortativities = []
#     number_swaps = []
#     iter = 0
#     for i in range(15000):
#
#         G_prime = G.copy()
#
#         '''
#         must check if the selected edges are unique
#         must check if proposed edges don't already exist
#
#         selected edges: (u, v), (w, x)
#         d_0 = |k_u - k_v| + |k_w - k_x|
#         d_1 = |k_u - k_w| + |k_v - k_x|
#         d_2 = |k_u - k_x| + |k_v - k_w|
#
#         if d_0 > min(d_1, d_2):
#             take the min
#             add edge for max assortativity
#
#         if d_0 < max(d_1, d_2):
#             take the max
#             add edge for min assortativity
#         '''
#
#         # randomly select 2 edges and remove these from the graph
#         all_edges = [e for e in G_prime.edges]
#         print(f"\n\nlen(all_edges): {len(all_edges)}")
#         edges_to_swap = []
#         while len(edges_to_swap) != 2:
#             random_edge = random.choice(all_edges)
#             if not random_edge in edges_to_swap:
#                 edges_to_swap.append(random_edge)
#                 # must check if all nodes in edge list are unique
#
#         print(f"edges_to_swap: {edges_to_swap}")
#
#         # create a dictionary where key = edge and val = degree:
#         vertex_degrees = [] # use a list instead?? possible key overlap????????????????????????????????????????????????????????????????
#         for edge in edges_to_swap:
#             for vertex in edge:
#                 vertex_degrees.append((vertex, G_prime.degree[vertex]))
#
#         # sort the dictionary of vertex degrees by degree - max frst
#         vertex_degrees = sorted(vertex_degrees, key = lambda x: x[1], reverse = True)
#         print(f"vertex_degrees: {vertex_degrees}")
#         print(f"vertex_degrees[0]: {vertex_degrees[0]}")
#         print(f"vertex_degrees[1]: {vertex_degrees[1]}")
#         print(f"vertex_degrees[2]: {vertex_degrees[2]}")
#         print(f"vertex_degrees[3]: {vertex_degrees[3]}\n\n")
#
#         # add edges between selected edges for maximal assortativity:
#         G_prime.add_edge(vertex_degrees[0][0], vertex_degrees[1][0])
#         # print(f"edge added: {(vertex_degrees[0][0], vertex_degrees[1][0])}")
#         G_prime.add_edge(vertex_degrees[2][0], vertex_degrees[3][0])
#         # print(f"edge added: {vertex_degrees[2][0], vertex_degrees[3][0]}\n\n")
#
#         # remove the selected edges from the graph:
#         for edge in edges_to_swap:
#             G_prime.remove_edge(*edge)
#
#         G = G_prime
#
#         # increment the iterations:
#         iter += 1
#
#     # plt.scatter(number_swaps, changed_assortativities, s = 1.5, color = 'darkviolet')
#     # plt.title("BA Graph Assortativity vs number of 2-OPT swaps")
#     # plt.xlabel("number 2-OPT swaps")
#     # plt.ylabel("assortativity")
#     # plt.show()
#
#     return G

'''
select 2 edges from G_prime: (x, y) (u, v)
find the degrees of x, y, u, v
delete the edges (x, y), (u, v)
connect the highest degree vertex with the second hightes degree vertex
and connect the lowest degree vertex with the second lowest degree vertex
'''
def assortative_swaps(G):
    changed_assortativities = []
    number_swaps = []
    iter = 0
    for i in range(15000):

        G_prime = G.copy()

        '''
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

        # randomly select 2 edges and remove these from the graph
        all_edges = [e for e in G_prime.edges]
        print(f"\n\nlen(all_edges): {len(all_edges)}")
        edges_to_swap = []
        while len(edges_to_swap) != 2:
            random_edge = random.choice(all_edges)
            if not random_edge in edges_to_swap and not [v for v in random_edge] in [v for v in all_edges]:
                edges_to_swap.append(random_edge)
                # must check if all nodes in edge list are unique do the max/min checks here

        print(f"edges_to_swap: {edges_to_swap}")

        # # create a dictionary where key = edge and val = degree:
        # vertex_degrees = [] # use a list instead?? possible key overlap????????????????????????????????????????????????????????????????
        # for edge in edges_to_swap:
        #     for vertex in edge:
        #         vertex_degrees.append((vertex, G_prime.degree[vertex]))
        #
        # # sort the dictionary of vertex degrees by degree - max frst
        # vertex_degrees = sorted(vertex_degrees, key = lambda x: x[1], reverse = True)
        # print(f"vertex_degrees: {vertex_degrees}")
        # print(f"vertex_degrees[0]: {vertex_degrees[0]}")
        # print(f"vertex_degrees[1]: {vertex_degrees[1]}")
        # print(f"vertex_degrees[2]: {vertex_degrees[2]}")
        # print(f"vertex_degrees[3]: {vertex_degrees[3]}\n\n")

        # edges_to_swap: [(u, v), (w, x)]
        d_0 = np.abs(edges_to_swap[0][0] - edges_to_swap[0][1]) + np.abs(edges_to_swap[1][0] - edges_to_swap[1][1]) # d_0 = |k_u - k_v| + |k_w - k_x|
        d_1 = np.abs(edges_to_swap[0][0] - edges_to_swap[1][0]) + np.abs(edges_to_swap[0][1] - edges_to_swap[1][1]) # d_1 = |k_u - k_w| + |k_v - k_x|
        d_2 = np.abs(edges_to_swap[0][0] - edges_to_swap[1][1]) + np.abs(edges_to_swap[0][1] - edges_to_swap[1][0]) # d_2 = |k_u - k_x| + |k_v - k_w|
        min_d = min(d_1, d_2)

        if d_0 > min_d:
            if min_d == d_1:
                # add edges between selected edges for maximal assortativity:
                G_prime.add_edge(edges_to_swap[0][0], edges_to_swap[1][0])
                # print(f"edge added: {(vertex_degrees[0][0], vertex_degrees[1][0])}")
                G_prime.add_edge(edges_to_swap[0][1], edges_to_swap[1][1])

                # remove the selected edges from the graph:
                for edge in edges_to_swap:
                    G_prime.remove_edge(*edge)

            elif min_d == d_2:
                # add edges between selected edges for maximal assortativity:
                G_prime.add_edge(edges_to_swap[0][0], edges_to_swap[1][1])
                # print(f"edge added: {(vertex_degrees[0][0], vertex_degrees[1][0])}")
                G_prime.add_edge(edges_to_swap[0][1], edges_to_swap[1][0])

                # remove the selected edges from the graph:
                for edge in edges_to_swap:
                    G_prime.remove_edge(*edge)

            else:
                print("--------------------------- NONE SWAPPED ---------------------------")

        # # add edges between selected edges for maximal assortativity:
        # G_prime.add_edge(vertex_degrees[0][0], vertex_degrees[1][0])
        # # print(f"edge added: {(vertex_degrees[0][0], vertex_degrees[1][0])}")
        # G_prime.add_edge(vertex_degrees[2][0], vertex_degrees[3][0])
        # # print(f"edge added: {vertex_degrees[2][0], vertex_degrees[3][0]}\n\n")

        # # remove the selected edges from the graph:
        # for edge in edges_to_swap:
        #     G_prime.remove_edge(*edge)

        G = G_prime

        # increment the iterations:
        iter += 1

    # plt.scatter(number_swaps, changed_assortativities, s = 1.5, color = 'darkviolet')
    # plt.title("BA Graph Assortativity vs number of 2-OPT swaps")
    # plt.xlabel("number 2-OPT swaps")
    # plt.ylabel("assortativity")
    # plt.show()

    return G


def disassortative_swaps_random(G):
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

# 2.d
def sandpile(end):
    # end must be >= 4
    lattice = nx.random_regular_graph(4, 400)
    sandpile = {}

    # create the sandpile:
    for node in range(0, 401):
        sandpile[node] = [grain for grain in range(0, end + 1)]


''' Assortativity '''
def show_assortativities(probabilities, average_degrees, assortativities):
    # plot assortativity against p
    plt.scatter(probabilities, assortativities, s = 2, color = 'blue')
    plt.title("assortativity against p")
    plt.xlabel("probability of edge connection")
    plt.ylabel("assortativity")
    plt.show()

    # plot assortativity against average degree
    plt.scatter(average_degrees, assortativities, s = 2, color = 'blue')
    plt.title("assortativity against average degree")
    plt.xlabel("average degree")
    plt.ylabel("assortativity")
    plt.show()

''' Transitivity '''
def show_transitivities(probabilities, average_degrees, transitivities):
    # plot transitivity against p
    plt.scatter(probabilities, transitivities, s = 2, color = 'red')
    plt.title("transitivity against p")
    plt.xlabel("probability of edge connection")
    plt.ylabel("transitivity")
    plt.show()

    # plot transitivity against p
    plt.scatter(average_degrees, transitivities, s = 2, color = 'red')
    plt.title("transitivity against average degree")
    plt.xlabel("average degree")
    plt.ylabel("transitivity")
    plt.show()

''' Degree Centrality '''
def show_degree_centralities(probabilities, average_degrees, degree_centralities):
    # plot degree_centrality against average degree
    plt.scatter(probabilities, degree_centralities, s = 2, color = 'green')
    plt.title("degree centrality against p")
    plt.xlabel("probability of edge connection")
    plt.ylabel("degree centrality")
    plt.show()

    # plot degree_centrality against average degree
    plt.scatter(average_degrees, degree_centralities, s = 2, color = 'green')
    plt.title("degree centrality against average degree")
    plt.xlabel("average degree")
    plt.ylabel("degree centrality")
    plt.show()

''' Eigenvector Centrality '''
def show_eigenvector_centralities(probabilities, average_degrees, eigenvector_centralities):
    # plot eigenvector_centralities against p
    plt.scatter(probabilities, eigenvector_centralities, s = 2, color = 'purple')
    plt.title("eigenvector centrality against p")
    plt.xlabel("probability of edge connection")
    plt.ylabel("eigenvector centrality")
    plt.show()

    # plot eigenvector_centralities against average degree
    plt.scatter(average_degrees, eigenvector_centralities, s = 2, color = 'purple')
    plt.title("eigenvector centrality against average degree")
    plt.xlabel("average degree")
    plt.ylabel("eigenvector centrality")
    plt.show()

''' percolation Centrality '''
def show_percolation_centralities(probabilities, average_degrees, percolation_centralities):
    # plot eigenvector_centralities against p
    plt.scatter(probabilities, percolation_centralities, s = 2, color = 'orange')
    plt.title("percolation centrality against p")
    plt.xlabel("probability of edge connection")
    plt.ylabel("percolation centrality")
    plt.show()

    # plot percolation_centralities against average degree
    plt.scatter(average_degrees, percolation_centralities, s = 2, color = 'orange')
    plt.title("percolation centrality against average degree")
    plt.xlabel("average degree")
    plt.ylabel("percolation centrality")
    plt.show()

''' Diameter'''
def show_diameters(probabilities, average_degrees, diameters):
    # plot diameters against p
    plt.scatter(probabilities, diameters, s = 1.5, color = 'lime')
    plt.title("graph diameter against p")
    plt.xlabel("probability of edge connection")
    plt.ylabel("diameter")
    plt.show()

    # plot diameters against average degree
    plt.scatter(average_degrees, diameters, s = 1.5, color = 'lime')
    plt.title("graph diameter against average degree")
    plt.xlabel("average degree")
    plt.ylabel("diameter")
    plt.show()


if __name__ == '__main__':

    # # question 2.a
    # generate_ER_graph_metrics()
    # show_degree_centralities(probabilities, average_degrees, degree_centralities)
    # show_assortativities(probabilities, average_degrees, assortativities)
    # show_transitivities(probabilities, average_degrees, transitivities)
    # show_eigenvector_centralities(probabilities, average_degrees, eigenvector_centralities)
    # show_percolation_centralities(probabilities, average_degrees, percolation_centralities)
    # show_diameters(probabilities, average_degrees, diameters)

    # # question 2.b
    # generate_WS_graph_metrics()
    # show_degree_centralities(WS_edge_rewire_probabilities, WS_average_degrees, WS_degree_centralities)
    # show_assortativities(WS_edge_rewire_probabilities, WS_average_degrees, WS_assortativities)
    # show_transitivities(WS_edge_rewire_probabilities, WS_average_degrees, WS_transitivities)
    # show_eigenvector_centralities(WS_edge_rewire_probabilities, WS_average_degrees, WS_eigenvector_centralities)
    # show_percolation_centralities(WS_edge_rewire_probabilities, WS_average_degrees, WS_percolation_centralities)
    # print(f"\n\nlen(WS_edge_rewire_probabilities): {len(WS_edge_rewire_probabilities)}")
    # print(f"len(WS_average_degrees): {len(WS_average_degrees)}\n\n")
    # show_diameters(WS_edge_rewire_probabilities, WS_average_degrees, WS_diameters)

    # BA Graph G_0
    # G_0 = nx.barabasi_albert_graph(300, 6)
    G_0 = nx.barabasi_albert_graph(300, 6)

    G_a = assortative_swaps(G_0)

    # assortativity of G_0:
    A_1 = nx.degree_assortativity_coefficient(G_0)

    # Assortativity of G_a:
    A_2 = nx.degree_assortativity_coefficient(G_a)


    # # The comparatively assortative graph:
    # G_a = assortative_swaps_random(G_0)

    # # The comparatively disassortative graph:
    # G_b = disassortative_swaps_random(G_0)
    #
    # # Assortativity of G_a:
    # A_2 = nx.degree_assortativity_coefficient(G_a)
    #
    nx.draw(G_0, node_size = 1.0, width = 0.5, node_color = "cyan", edge_color = "darkviolet")
    plt.show()
    nx.draw(G_a, node_size = 1.0, width = 0.5, node_color = "cyan", edge_color = "darkviolet")
    plt.show()

    print(f"\n\nA_1: {A_1}")
    print(f"A_2 {A_2}\n\n")

    # lattice = nx.random_regular_graph(4, 400)
    # nx.draw(lattice, node_size = 2.0, width = 1.0, node_color = "lime", edge_color = "green")
    # plt.show()
