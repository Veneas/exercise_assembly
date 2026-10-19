# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Old functions from previous exercise
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def degrees(graph):
    """
    Calculating the degree of every vertex in a undirected graph.
    Number of connection vertex has.
    :param graph:
    :return:
    """

    # Dict comprehension to map each vertex to the length of its neighbour
    return {vertex: len(neighbour  s) for vertex, neighbours in graph.items()}


def in_out_degree(graph):
    """
    Calculating the in-degree and out-degree for every vertex in a directed graph.
    Returns a tuple of two dicts: (in_degree, out_degree)
    :param graph:
    :return:
    """
    # Initialising counts for all vertices to 0
    in_degree = {vertex: 0 for vertex in graph}

    # Out-degree is simply the outgoing edges a vertex has
    out_degree = {vertex: len(neighbours) for vertex, neighbours in graph.items()}

    # To find in-degrees, checking every edge in the graph and counting how
    # many items each vertex is a "target" (v)
    for u, neighbours in graph.items():
        for v in neighbours:
            in_degree[v] += 1

    return in_degree, out_degree


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def dfs(graph, start):
    """
    A depth-first search starting from a given vertex.
    Return a list of vertices in the exact order they were first visited.
    :param graph:
    :param start:
    :return:
    """

    visited = set()     # Lookup to check if the vertex was seen
    order = []          # Keeping track of the exact visitation seq

    # Defining recursive helper func to do the walking
    def explore(vertex):
        visited.add(vertex)
        order.append(vertex)

        # Tie-break rule: sorting neighbours alphabetically
        for neighbour in sorted(graph[vertex]):
            # Only visiting if havent been there yet
            if neighbour not in visited:
                explore(neighbour)

    # Kicking off the recursion from the start vertex
    explore(start)

    return order


def connected_components(graph):
    """
    Finds all maximal connected components in an undirected graph.
    Returns a list of sets, where each set contains the vertices of one component.
    :param graph:
    :return:
    """

    visited_global = set()
    components = []

    # Checking every vertex in the graph
    for vertex in graph:
        # If havents visited this vettex in any previous DFS run
        if vertex not in visited_global:
            # start a new DFS to find everything connected to it
            reachable = dfs(graph, vertex)

            # Adding all newly found vertices to global visited list
            visited_global.update(reachable)

            # Storing this components as a set
            components.append(set(reachable))

    return components


def is_weakly_connected(graph):
    """
    Checks if a directed graph is weakly connected.
    It does this by "undirected shadow" of the graph and checking if that shadow has
    exactly 1 connected component.
    :param graph:
    :return:
    """
    if not graph:
        return True

    # Inizialising a new shadow graph with empty lists of all vertices
    shadow = {vertex: [] for vertex in graph}

    #Populating the shadow graph with undirected edges
    for u, neighbours in graph.items():
        for v in neighbours:
            shadow[u].append(v)
            shadow[v].append(u) # Adding reverse edge to get undirected


     # Checking the connected components of shadow graph
    components = connected_components(shadow)

    # its weakly connected if all vertices from exactly one large group
    return len(components) == 1

# =-=-=-=-=-=-=-=-=-=-=-= Task 4 =-=-=-=-=-=-=-=-=-=-=-=

def eulerian_type(graph, directed=False):
    """
    Determines if a graph has an Eulerian circuit, path, or none.
    :param graph:
    :param directed:
    :return:
    """

    if not directed:
        # Undirected graph logic
        # Checking connectivity
        comps = connected_components(graph)
        degs = degrees(graph)

        # Caring only about components that actually have edges.
        # Isolated vertices dont break Eulerian paths
        active_comps = [c for c in comps if any(degs[v] > 0 for v in c)]
        if len(active_comps) > 1:
            return "none"

        # Counting vertices with an odd degree
        odd_count = sum(1 for v, d in degs.items() if d % 2 != 0)

        # Applying eulerian rules
        if odd_count == 0:
            return "circuit"
        elif odd_count == 2:
            return "path"
        else:
            return "none"
    else:
        # Directed graph logic

        # Checking weak connectivity
        in_deg, out_deg = in_out_degree(graph)

        # Similarly to undirectred, ensuring all edges from a signle connected component
        # builing a shawdow to check active components
        shadow = {v: [] for v in graph}
        for u, neighs in graph.items():
            for v in neighs:
                shadow[u].append(v)
                shadow[v].append(u)

        comps = connected_components(shadow)
        active_comps = [c for c in comps if any(len(shadow[v]) > 0 for v in c)]
        if len(active_comps) > 1:
            return "none"

        # COunting degree imbalances
        start_nodes = 0
        end_nodes = 0

        for v in graph:
            diff = out_deg[v] - in_deg[v]

            if diff == 0:
                continue
            elif diff == 1:
                start_nodes += 1
            elif diff == -1:
                end_nodes += 1
            else:
                # If any vertex has an excess, an eulerian path is impossible
                return "none"
        # Applying eulerian rules
        if start_nodes == 0 and end_nodes == 0:
            return "circuit"
        elif start_nodes == 1 and end_nodes == 1:
            return "path"
        else:
            return "None"


    # =-=-=-=-=-=-=-=-=-=-=-= Task 2 =-=-=-=-=-=-=-=-=-=-=-=

def kmer_composition(text, k):
    """
    Returns every substring of length k from the text, in order.
    :param text:
    :param k:
    :return:
    """

    # Range ensures not going out of bounds
    return [text[i:i+k] for i in range(len(text)-k+1)]


def debruijn_graph(kmers):
    """
    Builds a directed de Bruijn graph from a list of kmers.
    Each k-mer represents an edge from its prefix to its suffix.
    :param kmers:
    :return:
    """

    graph = {}

    for kmer in kmers:
        prefix = kmer[:-1] # All letter except the last one
        suffix = kmer[1:] # All letter except the first one

        # Ensuring every prefix and suffix exists as key, even if empty
        if prefix not in graph:
            graph[prefix] = []
        if suffix not in graph:
            graph[suffix] = []

        # the k-mer bridges the prefix to the suffix
        graph[prefix].append(suffix)

    return graph


def spell_path(nodes):
    """
    Spells out the string represented by a walk through a de Bruijn graph.
    :param nodes:
    :return:
    """

    if not nodes:
        return ""

    # First full node (start)
    spelled = nodes[0]

    # For every subseq node, adding last character
    for node in nodes[1:]:
        spelled += node[-1]

    return spelled


    # =-=-=-=-=-=-=-=-=-=-=-= Task 3 =-=-=-=-=-=-=-=-=-=-=-=

def hierholzer(graph, start, verbose=False):
    """
    Finds an Eulerian path using Hierholzer's algorithm.
    :param graph:
    :param start:
    :param verbose:
    :return:
    """
    # Creating a copy of the graph to keep track of unused edges
    # without modifying the original!
    # Sorting the lists alphabetically (ascending) ensuring that pop(0)
    # always gives the alphabetically first available node.

    unused_edges = {u: sorted(v_list) for u, v_list in graph.items()}

    # Helper walk routine
    def walk_from(node):
        new_walk = [node]
        curr = node
        # Keeping walk until the current node has no unused outgoing edges
        while unused_edges.get(curr):
            next_node = unused_edges[curr].pop(0)
            new_walk.append(next_node)
            curr = next_node
        return new_walk

    # Initial walk from the start node
    walk = walk_from(start)

    if verbose:
        print(f"Walk 1 (from the start node): {walk}")

    # Looping until every edge is used
    # any() checks if there are still edges left in any of the lists
    while any(unused_edges.values()):
        # Finding the first node in the curr walk that still has unused edges
        restart_index = -1
        restart_node = ""

        for i, node in enumerate(walk):
            if unused_edges.get(node):
                restart_index = i
                restart_node = node
                break

        # Performing a new walk starting from this node.
        new_path = walk_from(restart_node)

        if verbose:
            print(f"Restart at: {restart_node} walk: {new_path}")

        # Splicing the new walk into the old walk
        # Replacing single instance of the restart node with the entire new_path
        walk = walk[:restart_index] + new_path + walk[restart_index + 1:]

        if verbose:
            print(f"Walk after splicing: {walk}")

    return walk

def eulerian_path(graph, verbose=False):
    """
    Determines the correct starting node, then runs Hierholzer's algorithm.
    :param graph:
    :param verbose:
    :return:
    """

    in_deg, out_deg = in_out_degree(graph)

    start_node = None

    # Checking for a node with out_deg - in_deg == 1 (start of path)
    for v in graph:
        if out_deg.get(v,0) - in_deg.get(v, 0) == 1:
            start_node = v
            break
    # If all nodes are balanced, its an eulerian circuit.
    # the rule says to start at the alphabetically first node.

    if start_node is None:
        start_node = min(graph.keys())

    return hierholzer(graph, start_node, verbose=verbose)


# # =-=-=-=-=-=-=-=-=-= Main =-=-=-=-=-=-=-=-=-=-=-=-=-=
if __name__ == "__main__":
    print("=-=-=-=-=-=-=-=-=-= Task 2 =-=-=-=-=-=-=-=-=-=-=-=-=-= ")
    test_string = "GATTACA"
    k = 4

    kmers = kmer_composition(test_string, k)
    print(f"{k}-mers of {test_string}: {kmers}")

    graph = debruijn_graph(kmers)
    print("\n De Bruijn Graph:")
    for node, edges in graph.items():
        print(f"{node}: {edges}")

    # Testing spell path
    nodes_walk = [kmers[0][:-1]] + [kmer[1:] for kmer in kmers]
    print(f"\nSpelled path: {spell_path(nodes_walk)}")

    print("=-=-=-=-=-=-=-=-=-= Task 3 =-=-=-=-=-=-=-=-=-=-=-=-=-= ")
    from collections import Counter

    toy_b_kmers = ["ACT", "AGC", "ATT", "CCG", "CGA", "CGT", "CTG", "GAC",
                   "GAG", "GAT", "GTC", "TCC", "TCG", "TGA", "TTC"]

    graph_b = debruijn_graph(toy_b_kmers)

    # running the pathfinding with verbose to compare with the hand trace
    path = eulerian_path(graph_b, verbose=True)

    print("\nFinal spelled genome:")
    genome = spell_path(path)
    print(genome)

    # Self-checks
    edges_count = sum(len(v) for v in graph_b.values())
    print("\n Self Checks")
    print(f"1. Path nodes == edges + 1? {len(path) == edges_count + 1}")

    path_edges = list(zip(path, path[1:]))
    graph_edges = [(u,v) for u,v_list in graph_b.items() for v in v_list]
    print(f"2. Every edge used exactly once? {Counter(path_edges) == Counter(graph_edges)}")

    spelled_kmers = kmer_composition(genome, 3)
    print(f"3. Spelled kmers match original k-mers ? {Counter(spelled_kmers) == Counter(toy_b_kmers)}")

