import networkx as nx
import matplotlib.pyplot as plt
import numpy as np


def average_degree(G):
    return G.number_of_edges() / 300

def average_degree_centrality(G):
    return sum(nx.degree_centrality(G).values())/len(nx.degree_centrality(G).values())

def average_eigenvector_centrality(G):
    return sum(nx.eigenvector_centrality(G).values())/len(nx.eigenvector_centrality(G).values())

def average_percolation_centrality(G):
    nx.set_node_attributes(G, 0.1, 'percolation')
    percolation_cent_dict = nx.percolation_centrality(
        G = G,
        attribute='percolation',
    )
    return sum(nx.percolation_centrality(G).values())/len(nx.percolation_centrality(G).values())

assortativities = []
probabilities = []
transitivities = []
degree_centralities = []
eigenvector_centralities = []
percolation_centralities = []
average_degrees = []
diameters = []

# compute all relevant metrics for 100 random graphs with harmoically
# decreasing edge connection probabilities
for i in range(1, 100):
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
        # diams[i] = [nx.diameter(G), False]
    else:
        diameters.append(np.Infinity)
        # diameters.append(diameters[-1])
        # diams[i] = [nx.diameter(G), True]


''' Assortativity '''
# # plot assortativity against p
# plt.plot(probabilities, assortativities)
# plt.title("assortativity against p")
# plt.xlabel("probability of edge connection")
# plt.ylabel("assortativity")
# plt.show()
#
# # plot assortativity against average degree
# plt.plot(average_degrees, assortativities)
# plt.title("assortativity against average degree")
# plt.xlabel("probability of edge connection")
# plt.ylabel("average degree")
# plt.show()

''' Transitivity '''
# # plot transitivity against p
# plt.plot(probabilities, transitivities)
# plt.title("transitivity against p")
# plt.xlabel("probability of edge connection")
# plt.ylabel("transitivity")
# plt.show()
#
# # plot transitivity against p
# plt.plot(average_degrees, transitivities)
# plt.title("transitivity against average degree")
# plt.xlabel("probability of edge connection")
# plt.ylabel("assortativity")
# plt.show()

''' Degree Centrality '''
# # plot degree_centrality against average degree
# plt.plot(probabilities, degree_centralities)
# plt.title("degree centrality against p")
# plt.xlabel("probability of edge connection")
# plt.ylabel("degree centrality")
# plt.show()
#
# # plot degree_centrality against average degree
# plt.plot(average_degrees, degree_centralities)
# plt.title("degree centrality against average degree")
# plt.xlabel("average degree")
# plt.ylabel("degree centrality")
# plt.show()

''' Eigenvector Centrality '''
# # plot eigenvector_centralities against p
# plt.plot(probabilities, eigenvector_centralities)
# plt.title("eigenvector centrality against p")
# plt.xlabel("probability of edge connection")
# plt.ylabel("eigenvector centrality")
# plt.show()
#
# # plot eigenvector_centralities against average degree
# plt.plot(average_degrees, eigenvector_centralities)
# plt.title("eigenvector centrality against average degree")
# plt.xlabel("average degree")
# plt.ylabel("eigenvector centrality")
# plt.show()

''' percolation Centrality '''
# plot eigenvector_centralities against p
plt.plot(probabilities, percolation_centralities)
plt.title("percolation centrality against p")
plt.xlabel("probability of edge connection")
plt.ylabel("percolation centrality")
plt.show()

# plot percolation_centralities against average degree
plt.plot(average_degrees, percolation_centralities)
plt.title("percolation centrality against average degree")
plt.xlabel("average degree")
plt.ylabel("percolation centrality")
plt.show()

''' Diameter'''
# # plot diameters against p
# plt.scatter(probabilities, diameters, s = 1.5, color = 'blue')
# plt.title("graph diameter against p")
# plt.xlabel("probability of edge connection")
# plt.ylabel("diameter")
# plt.show()
#
# # plot diameters against average degree
# plt.scatter(average_degrees, diameters, s = 1.5, color = 'red')
# plt.title("graph diameter against average degree")
# plt.xlabel("average degree")
# plt.ylabel("diameter")
# plt.show()




# # GENERATE EXEMPLAR P VALUE VERSION
# # 2.a - Erdos-Renyi Graph with 300 vertices and 0.5 probability edge conection
# assortativities = []
# probabilities = []
# transitivities = []
# degree_centralities = []
#
# # compute all relevant metrics for 100 random graphs with harmoically
# # decreasing edge connection probabilities
# G = nx.erdos_renyi_graph(300, 0.5)
# probabilities.append(1/i)
# # degrees = G.degree()
# # assortativities.append(nx.degree_assortativity_coefficient(G))
# # transitivities.append(nx.transitivity(G))
# degree_centralities.append(nx.degree_centrality(G).values()) # how is this possible?
#
# # # plot assortativity against p
# # plt.plot(probabilities, assortativities)
# # plt.xlabel("probability of edge connection")
# # plt.ylabel("assortativity")
# # plt.show()
# #
# # # plot transitivity against p
# # plt.plot(probabilities, transitivities)
# # plt.xlabel("probability of edge connection")
# # plt.ylabel("transitivity")
# # plt.show()
#
# # # plot degree_centrality against p
# # plt.plot(probabilities, degree_centralities)
# # plt.xlabel("probability of edge connection")
# # plt.ylabel("degree centrality")
# # plt.show()


# # 2.b Small World graph
# # repeat as in a but graph WRT rewiring probability
# G = nx.watts_strogatz_graph(300, 4, p)
