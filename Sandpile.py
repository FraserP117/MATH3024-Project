import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import random
from itertools import chain


'''
Initialises the sandpile
'''
def create_sandpile(num_grains):
    sandpile = {}
    for node in range(0, num_grains):
        # sandpile[node] = random.randint(0, 3)
        sandpile[node] = 0

    return sandpile

'''
Return True iff the entire sandpile is saturated.
'''
def is_saturated(sandpile, threshold):
    for node, num_grains in sandpile.items():
        if sandpile[node] < threshold:
            return False
    return True

'''
Add num_grains_to_add grains of sand to the sandpile, at a random vertex in the lattice.
'''
def add_sand(num_grains_to_add, sandpile, lattice):
    nodes = list(lattice.nodes())
    node_to_add_sand = random.choice(nodes)
    sandpile[node_to_add_sand] = sandpile[node_to_add_sand] + num_grains_to_add
    print(f"{sandpile[node_to_add_sand]} grains on node: {node_to_add_sand}")
    return sandpile

'''
return True iff there exists at least 1 critical node in the lattice.
node is critical iff its number of grains > threshold and at least
one of its neighbors has less than threshold grains.
'''
def any_critical(sandpile, threshold, lattice):
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
        print(f"number critical: {num_critical}")
        return True

'''
Returns True iff node is critical; False otherwise.
'''
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
def simulate_sandpile(threshold, lattice, sandpile):
    num_avalanches = 0
    ctr = 0
    while any_critical(sandpile, threshold, lattice): # getting stuck here
        for node, num_grains in sandpile.items():
            if is_critical(node, sandpile, threshold):
                num_avalanches += 1
                sandpile[node] -= 4
                # print(f"critical node: {node}, num_grains: {num_grains}")
                neighbors = [n for n in lattice.neighbors(node)]
                # print(f"neighbors: {neighbors}\n\n")
                for neighbor in neighbors:
                    sandpile[neighbor] += 1

        if is_saturated(sandpile, threshold):
            print("SATURATED")
            # break
            return num_avalanches
        else:
            print("NOT SATUTRATED")

    # print(f"sandpile: {sandpile.values()}")
    return num_avalanches


if __name__ == '__main__':
    avalanche_counter = 0
    number_avalanches = []
    threshold = 4
    num_grains = 400
    lattice = nx.random_regular_graph(4, num_grains)
    sandpile = create_sandpile(num_grains)
    iter = 0
    iters = []

    # iterate adding grains now
    for i in range(800):
        print(iter)
        critical_sandpile = add_sand(1, sandpile, lattice)
        num_av = simulate_sandpile(4, lattice, sandpile)
        number_avalanches.append(num_av)
        avalanche_counter += num_av
        iters.append(iter)
        iter += 1

    # plt.plot(iters, number_avalanches, linewidth = 0.75, color = "green")
    (markers, stemlines, baseline) = plt.stem(iters, number_avalanches)
    plt.title("Avalanche Size Per Iteration")
    plt.xlabel("Iteration")
    plt.ylabel("Avalanche Size")
    plt.setp(baseline, visible = False)
    plt.setp(stemlines, color = "lime", linewidth = 0.5)
    plt.setp(markers, color = "green", markersize = 0.75)
    plt.show()

    print(f"number of total avalanches: {avalanche_counter}")

    # nx.draw(lattice, node_size = 2.0, width = 0.75, node_color = "lime", edge_color = "navy")
    # plt.show()
