import networkx as nx
import matplotlib.pyplot as plt
import numpy as np


assortativities = []
probabilities = []
transitivities = []
degree_centralities = []
eigenvector_centralities = []
percolation_centralities = []
average_degrees = []
diameters = []

WS_assortativities = []
WS_probabilities = []
WS_transitivities = []
WS_degree_centralities = []
WS_eigenvector_centralities = []
WS_percolation_centralities = []
WS_average_degrees = []
WS_diameters = []

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
def generate_ER_graph_metrics():
    for i in range(1, 101):
        G = nx.erdos_renyi_graph(300, 1/i)
        probabilities.append(1/i)
        average_degrees.append(round(average_degree(G)))
        assortativities.append(nx.degree_assortativity_coefficient(G))
        transitivities.append(nx.transitivity(G))
        degree_centralities.append(average_degree_centrality(G))
        eigenvector_centralities.append(average_eigenvector_centrality(G))
        percolation_centralities.append(average_percolation_centrality(G))

        if nx.is_connected(G):
            diameters.append(nx.diameter(G))
        else:
            diameters.append(np.Infinity)

# compute all relevant metrics for 100 random graphs with harmoically
# decreasing edge connection probabilities
def generate_WS_graph_metrics():
    for i in range(1, 10):
        G = nx.watts_strogatz_graph(300, 4, 1/i)
        WS_probabilities.append(1/i)
        WS_average_degrees.append(round(average_degree(G)))
        WS_assortativities.append(nx.degree_assortativity_coefficient(G))
        WS_transitivities.append(nx.transitivity(G))
        WS_degree_centralities.append(average_degree_centrality(G))
        WS_eigenvector_centralities.append(average_eigenvector_centrality(G))
        WS_percolation_centralities.append(average_percolation_centrality(G))

        if nx.is_connected(G):
            diameters.append(nx.diameter(G))
        else:
            diameters.append(np.Infinity)

changed = []
itrs = []
def max_assortative_swaps(G):
    epsilon = 0.1
    delta = np.Inf
    G_prime = G.copy()
    iter = 0
    for i in range(15000):
        A = nx.degree_assortativity_coefficient(G)
        G_prime = G.copy()
        nx.double_edge_swap(G_prime)
        if nx.degree_assortativity_coefficient(G_prime) > nx.degree_assortativity_coefficient(G):
            G = G_prime
            changed.append(nx.degree_assortativity_coefficient(G_prime))
            itrs.append(iter)
            iter += 1
        # else:
        #     changed.append(0)

        # itrs.append(i)

    return G

# # compute all relevant metrics for 100 random graphs with harmoically
# # decreasing edge connection probabilities
# def generate_BA_graph_metrics():
#     for i in range(1, 10):
#         G = nx.barabasi_albert_graph(300, 6)
#         WS_probabilities.append(1/i)
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
    plt.xlabel("probability of edge connection")
    plt.ylabel("average degree")
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
    plt.xlabel("probability of edge connection")
    plt.ylabel("assortativity")
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

    # question 2.b
    # generate_WS_graph_metrics()
    # show_degree_centralities(WS_probabilities, WS_average_degrees, WS_degree_centralities)
    # show_assortativities(WS_probabilities, WS_average_degrees, WS_assortativities)
    # show_transitivities(WS_probabilities, WS_degree_centralities, WS_assortativities)
    # show_eigenvector_centralities(WS_probabilities, WS_degree_centralities, WS_assortativities)
    # show_percolation_centralities(WS_probabilities, WS_degree_centralities, WS_assortativities)
    # show_diameters(WS_probabilities, WS_degree_centralities, WS_assortativities)

    # BA Graph G_0
    G_0 = nx.barabasi_albert_graph(300, 6)
    #
    # # assortativity of G_0:
    # A_1 = nx.degree_assortativity_coefficient(G_0)
    #
    # # The comparatively assortative graph:
    # G_1 = max_assortative_swaps(G_0)
    #
    # # Assortativity of G_1:
    # A_2 = nx.degree_assortativity_coefficient(G_1)
    #
    # print(f"\n\nA_1: {A_1}")
    # print(f"A_2 {A_2}\n\n")

    G_1 = max_assortative_swaps(G_0)

    plt.scatter(itrs, changed, s = 1.0, color = 'green')
    plt.show()

    # assortativity of G_0:
    A_1 = nx.degree_assortativity_coefficient(G_0)

    # Assortativity of G_1:
    A_2 = nx.degree_assortativity_coefficient(G_1)

    nx.draw(G_0, node_size = 1.0, width = 0.5, node_color = "lime", edge_color = "green")
    plt.show()
    nx.draw(G_1, node_size = 1.0, width = 0.5, node_color = "lime", edge_color = "green")
    plt.show()

    print(f"\n\nA_1: {A_1}")
    print(f"A_2 {A_2}\n\n")
