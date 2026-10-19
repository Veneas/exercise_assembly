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
    return {vertex: len(neighbours) for vertex, neighbours in graph.items()}


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
    Depth-First Search rewritten to use a loop and an explicit stack
    to avoid recursion depth limits on massive genomes.
    """
    visited = set()
    order = []
    stack = [start]

    while stack:
        vertex = stack.pop()

        if vertex not in visited:
            visited.add(vertex)
            order.append(vertex)

            # We sort in REVERSE alphabetical order before adding to the stack,
            # so the alphabetically FIRST neighbor is at the top of the stack
            # and gets popped/visited next.
            # graph.get(vertex, []) prevents KeyErrors if a node has no outgoing edges
            for neighbour in sorted(graph.get(vertex, []), reverse=True):
                if neighbour not in visited:
                    stack.append(neighbour)

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

# =-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=

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


    # =-=-=-=-=-=-=-=-=-=-=-= Task 4 =-=-=-=-=-=-=-=-=-=-=-=


def contigs(graph):
    """
    Finds all maximal non-branching paths (contigs) in a de Bruijn graph.
    Returns a sorted list of the strings they spell.
    :param graph:
    """
    in_deg, out_deg = in_out_degree(graph)

    # Identifying all nodes that are NOT 1-in-1-out
    # These are branching nodes (starts and ends of contigs)
    non_1_in_1_out = set()
    for node in graph:
        if in_deg.get(node, 0) != 1 or out_deg.get(node, 0) != 1:
            non_1_in_1_out.add(node)

    # Keep track of unused edges to ensure every edge is used exactly once
    unused_edges = {u: list(v_list) for u, v_list in graph.items()}
    result_contigs = []

    # Extracting all paths that start from a branching node
    for start_node in non_1_in_1_out:
        # As long as this branching node has outgoing edges...
        while unused_edges.get(start_node):
            path = [start_node]
            curr = start_node

            # Take the first step along an unused edge
            next_node = unused_edges[curr].pop(0)
            path.append(next_node)
            curr = next_node

            # Keep walking forward as long as the node it lands on IS 1-in-1-out
            while curr not in non_1_in_1_out:
                next_node = unused_edges[curr].pop(0)
                path.append(next_node)
                curr = next_node

            # Hitting a branching node, so the contig ends here.
            result_contigs.append(spell_path(path))

    # Handling any isolated cycles left over
    # If the graph has isolated loops of only 1-in-1-out nodes, they
    # wouldn't be caught by step 2 because they have no branching nodes.
    for start_node in graph:
        while unused_edges.get(start_node):
            path = [start_node]
            curr = start_node

            # Walk the cycle until it returns to the start_node
            next_node = unused_edges[curr].pop(0)
            path.append(next_node)
            curr = next_node

            while curr != start_node:
                next_node = unused_edges[curr].pop(0)
                path.append(next_node)
                curr = next_node

            # The instructions require the cycle to start/end at the
            # alphabetically first node in the cycle.
            cycle_nodes = path[:-1]  # Remove the repeated end node for a moment
            first_node = min(cycle_nodes)
            idx = cycle_nodes.index(first_node)

            # Rotate the list so it starts and ends at the alphabetically first node
            rotated_path = cycle_nodes[idx:] + cycle_nodes[:idx] + [first_node]
            result_contigs.append(spell_path(rotated_path))

    # The instructions require returning the contigs sorted alphabetically
    return sorted(result_contigs)


    # =-=-=-=-=-=-=-=-=-=-=-= Task 5 =-=-=-=-=-=-=-=-=-=-=-=



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


    print("\n=-=-=-=-=-=-=-=-=-= Task 4 =-=-=-=-=-=-=-=-=-=-=-=-=-= ")
    toy_c_kmers = [
        "ATC", "CCT", "CTT", "GGT", "GTA", "GTG",
        "TAT", "TCT", "TGG", "TGT", "TTG"
    ]

    graph_c = debruijn_graph(toy_c_kmers)
    c_contigs = contigs(graph_c)

    print("Contigs found in Toy C:")
    for c in c_contigs:
        print(" -", c)

    print("\nMatches handwritten expectation?")
    expected = sorted(["CCT", "CTTG", "TGGT", "TGT", "GTG", "GTATCT"])
    print(c_contigs == expected)


    print("\n=-=-=-=-=-=-=-=-=-= Task 5 =-=-=-=-=-=-=-=-=-=-=-=-=-= ")
    from Bio import SeqIO
    import os

    # Reading the reads from the FASTA file
    try:
        reads = [str(record.seq) for record in SeqIO.parse("reads.fasta", "fasta")]
        print(f"Successfully loaded {len(reads)} reads.")
    except FileNotFoundError:
        print("Could not find 'reads.fasta'. Please ensure it is in the same directory.")
        exit()

    k_values = [9, 12, 13, 14, 21]

    # Storing the trimmed spelled strings to compare circularity later
    trimmed_genomes = {}
    contigs_dict = {}
    # Printing a header for our output table
    print(
        f"\n{'k':<4} | {'|V|':<6} | "
        f"{'|E|':<6} | {'Eulerian Type':<15} | "
        f"{'Contigs':<8} | {'Max Contig (Trimmed)':<20} | "
        f"{'Spelled String (Trimmed)'}")

    # Processing each k
    for k in k_values:
        # Collecting distinct k-mers using a set
        distinct_kmers_set = set()
        for read in reads:
            # kmer_composition from Task 2
            distinct_kmers_set.update(kmer_composition(read, k))

        # Passing sorted distinct k-mers to the graph builder
        sorted_kmers = sorted(list(distinct_kmers_set))
        graph = debruijn_graph(sorted_kmers)

        # Calculating |V| and |E|
        V = len(graph)
        E = sum(len(edges) for edges in graph.values())

        # Eulerian type
        e_type = eulerian_type(graph, directed=True)

        # Contigs
        raw_contigs = contigs(graph)
        num_contigs = len(raw_contigs)

        # Trim cycles for contigs
            # If a contig repeats its first k-1 letters at the end, it's a cycle -> Cut them off.
        trimmed_contigs = []
        for c in raw_contigs:
            if c[:k - 1] == c[-(k - 1):]:
                trimmed_contigs.append(c[:-(k - 1)])
            else:
                trimmed_contigs.append(c)

        contigs_dict[k] = trimmed_contigs
        longest_contig_len = max((len(c) for c in trimmed_contigs), default=0)

        # Finding Eulerian path if it exists
        spelled_len_str = "—"
        if e_type != "none":
            path = eulerian_path(graph)
            spelled = spell_path(path)

            # Trimming the cycle off the spelled genome
            if spelled[:k - 1] == spelled[-(k - 1):]:
                trimmed_spelled = spelled[:-(k - 1)]
            else:
                trimmed_spelled = spelled

            spelled_len_str = str(len(trimmed_spelled))
            trimmed_genomes[k] = trimmed_spelled

        # Printing the row for the table
        print(
            f"{k:<4} | {V:<6} | {E:<6} | "
            f"{e_type:<15} | {num_contigs:<8} | "
            f"{longest_contig_len:<20} | {spelled_len_str}")

    # Comparing the trimmed spelled strings (check for circularity)
    print("Checking if spelled strings describe the same circle...")
    # Getting all the generated genomes (should be for k=13, 14, 21)
    genomes_list = list(trimmed_genomes.values())

    if len(genomes_list) > 1:
        all_match = True
        base_genome = genomes_list[0]
        for other_genome in genomes_list[1:]:
            # Two strings describe the same circle if they have the same length
            # and one is a substring of the other concatenated with itself.
            if len(base_genome) != len(other_genome) or base_genome not in (other_genome + other_genome):
                all_match = False
                break
        print(f"Do all spelled strings describe the same circle? {all_match}")
    else:
        print("Not enough spelled strings to compare.")

    # Preparing sequences for BLAST
    print("\nSequences for BLAST")

    target_k = None
    for k in k_values:
        if len(contigs_dict[k]) == 1:
            target_k = k
            break

    if target_k is not None:
        print(f"\n(a) For k={target_k} (Smallest k with exactly 1 contig), here is the single trimmed contig to BLAST:")
        print(contigs_dict[target_k][0])

    print("\n(b) Trimmed spelled string of every smaller k at which you got one:")
    for k in k_values:
        if target_k and k < target_k and k in trimmed_genomes:
            print(f"\nSpelled string for k={k}:")
            print(trimmed_genomes[k])