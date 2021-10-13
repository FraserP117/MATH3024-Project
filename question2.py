import networkx as nx
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

'''
Do we calculate 10, 50, 100 diferent graphs and make the necessary meausrements for each,
graphing the results? or do we simply chose a few values for p as exemplars?

We graph centralities aginst average degree instead?
'''

def average_degree(G):
    degrees = dict(nx.degree(G)).values()
    return sum(degrees)/len(degrees)

# GENERATE 100 GRAPHS VERSION
# 2.a - Erdos-Renyi Graph with 300 vertices and 0.5 probability edge conection
assortativities = []
probabilities = []
transitivities = []
degree_centralities = []
eigenvector_centralities = []
average_degrees = []

# compute all relevant metrics for 100 random graphs with harmoically
# decreasing edge connection probabilities
for i in range(1, 10):
    G = nx.erdos_renyi_graph(300, 1/i)
    probabilities.append(1/i)
    # print(f"rounded avg degree: {round(average_degree(G))}")
    average_degrees.append(round(average_degree(G)))
    # degrees = G.degree()
    # assortativities.append(nx.degree_assortativity_coefficient(G))
    # transitivities.append(nx.transitivity(G))
    degree_centralities.append(nx.degree_centrality(G).values()) # how is this possible?
    # eigenvector_centralities.append(nx.eigenvector_centrality(G).values())


# # plot assortativity against p
# plt.plot(probabilities, assortativities)
# plt.xlabel("probability of edge connection")
# plt.ylabel("assortativity")
# plt.show()
#
# # plot transitivity against p
# plt.plot(probabilities, transitivities)
# plt.xlabel("probability of edge connection")
# plt.ylabel("transitivity")
# plt.show()

# print(average_degrees)
# print(degree_centralities)

# plot degree_centrality against average degree
plt.plot(average_degrees, degree_centralities)
plt.xlabel("average degree")
plt.ylabel("degree centrality")
plt.show()

# # plot eigenvector_centralities against p
# plt.plot(probabilities, eigenvector_centralities)
# plt.xlabel("probability of edge connection")
# plt.ylabel("eigenvector centrality")
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
