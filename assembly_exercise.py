# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Old functions from previous exercise
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

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

