def create_sandpile(num_grains):
    sandpile = {}
    # create the sandpile:
    # for node in range(0, 401):
    for node in range(0, num_grains):
        sandpile[node] = random.randint(0, 3)

    return sandpile

def add_sand(num_grains_to_add, sandpile):
    for node, num_grains in sandpile.items():
        if num_grains == 3:
            # num_grains += num_grains_to_add
            sandpile[node] = num_grains + num_grains_to_add
            print(f"{num_grains} grains on node: {node}")
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
        return True
    else:
        return False

'''
if is_critical(node, sandpile, threshold):
    subtract 4 grains from node
    neighbors = [n for n in lattice.neighbors(node)]
    for neighbor in neighbors:
        sandpile[neighbor] += 1
'''

num_avalanches = 0
def sandpile_sim(threshold, lattice, sandpile):
    print(f"sandpile: {sandpile.values()}")
    global num_avalanches
    while any_critical(sandpile, threshold):
        num_crit = 0
        for node, num_grains in sandpile.items():
            if is_critical(node, sandpile, threshold):
                num_avalanches += 1
                sandpile[node] -= 4
                print(f"\n\ncritical node: {node}, num_grains: {num_grains}")
                neighbors = [n for n in lattice.neighbors(node)]
                print(f"neighbors: {neighbors}\n\n")
                for neighbor in neighbors:
                    sandpile[neighbor] += 1
    print(f"sandpile: {sandpile.values()}")
    return num_avalanches
