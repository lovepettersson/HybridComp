import matplotlib.pyplot as plt
import networkx as nx
import json

def local_complementation_on_graphs(G, target):
    G_target = nx.ego_graph(G, target, 1, center=False)
    G.remove_edges_from(list(G_target.edges()))
    G_target_comp = nx.complement(G_target)
    G = nx.compose(G, G_target_comp)
    H = nx.Graph(G)
    return H


def X_measure(target, nb, graph):
    graph = local_complementation_on_graphs(graph, nb)
    graph = local_complementation_on_graphs(graph, target)
    graph.remove_node(target)
    graph = local_complementation_on_graphs(graph, nb)
    return graph


def gen_concat_ring_no_parity_code(num_nodes, numb_graphs):
    graph = nx.Graph()
    nodes = range(num_nodes)
    graph.add_nodes_from(nodes)
    iteration = [x for x in range(0, 30)]
    ring_edges = []
    ring_nb = []
    non_ring_edges = []
    for idx in range(numb_graphs):
        zero, one, two, three, four = iteration[idx * 5:(idx + 1) * (5)]
        edges = [(zero, one), (one, two), (two, three),
                 (three, four), (four, zero)]
        non_ring_edges.append(one)
        non_ring_edges.append(two)
        non_ring_edges.append(three)
        non_ring_edges.append(four)
        ring_edges.append(zero)
        ring_nb.append(one)
        graph.add_edges_from(edges)
    # graph.add_edges_from(edges)
    zero, one, two, three, four, five = ring_edges
    edges = [(zero, one), (one, two), (two, three), (three, four), (four, five), (five, zero)]
    graph.add_edges_from(edges)
    for idx in range(len(ring_nb)):
        graph = X_measure(ring_edges[idx], ring_nb[idx], graph)
    return graph



def gen_concat_ring(num_nodes, numb_graphs):
    graph = nx.Graph()
    nodes = range(num_nodes)
    graph.add_nodes_from(nodes)
    iteration = [x for x in range(48, 78)]
    ring_edges = []
    ring_nb = []
    non_ring_edges = []
    parity_qbts = [x for x in range(0, 48)]
    for idx in range(numb_graphs):
        zero, one, two, three, four = iteration[idx * 5:(idx + 1) * (5)]
        edges = [(zero, one), (one, two), (two, three),
                 (three, four), (four, zero)]
        non_ring_edges.append(one)
        non_ring_edges.append(two)
        non_ring_edges.append(three)
        non_ring_edges.append(four)
        ring_edges.append(zero)
        ring_nb.append(one)
        graph.add_edges_from(edges)
    # graph.add_edges_from(edges)
    zero, one, two, three, four, five = ring_edges
    edges = [(zero, one), (one, two), (two, three), (three, four), (four, five), (five, zero)]
    graph.add_edges_from(edges)
    for idx in range(len(ring_nb)):
        graph = X_measure(ring_edges[idx], ring_nb[idx], graph)

    for i in range(len(non_ring_edges)):
        p_edge_1 = parity_qbts[2 * i]
        p_edge_2 = parity_qbts[2 * i + 1]
        ring_node = non_ring_edges[i]
        edges = [(ring_node, p_edge_1), (p_edge_1, p_edge_2)]
        graph.add_edges_from(edges)
        graph = X_measure(ring_node, p_edge_1, graph)
    return graph



def gen_four_large_parity_code(qubits_per_graph, numb_graphs, concat_flag=True):
    graph = nx.Graph()
    iteration = [x for x in
                 range(4 * 6 * (qubits_per_graph - 1), 4 * 6 * (qubits_per_graph - 1) + 6 * qubits_per_graph)]
    ring_edges = []
    ring_nb = []
    non_ring_edges = []
    parity_qbts = [x for x in range(0, 4 * 6 * (qubits_per_graph - 1))]
    for idx in range(numb_graphs):
        zero, one, two, three, four = iteration[idx * qubits_per_graph:(idx + 1) * qubits_per_graph]
        edges = [(zero, one), (one, two), (two, three),
                 (three, four), (four, zero)]
        non_ring_edges.append(one)
        non_ring_edges.append(two)
        non_ring_edges.append(three)
        non_ring_edges.append(four)
        ring_edges.append(zero)
        ring_nb.append(one)
        graph.add_edges_from(edges)
    zero, one, two, three, four, five = ring_edges
    edges = [(zero, one), (one, two), (two, three), (three, four), (four, five), (five, zero)]
    graph.add_edges_from(edges)

    if concat_flag:
        for idx in range(len(ring_nb)):
            graph = X_measure(ring_edges[idx], ring_nb[idx], graph)
        for i in range(len(non_ring_edges)):
            p_edge_1 = parity_qbts[4 * i]
            p_edge_2 = parity_qbts[4 * i + 1]
            p_edge_3 = parity_qbts[4 * i + 2]
            p_edge_4 = parity_qbts[4 * i + 3]
            ring_node = non_ring_edges[i]
            edges = [(ring_node, p_edge_1), (p_edge_1, p_edge_2), (p_edge_1, p_edge_3), (p_edge_1, p_edge_4)]
            graph.add_edges_from(edges)
            graph = X_measure(ring_node, p_edge_1, graph)
    else:
        for i in range(len(non_ring_edges)):
            p_edge_1 = parity_qbts[4 * i]
            p_edge_2 = parity_qbts[4 * i + 1]
            p_edge_3 = parity_qbts[4 * i + 2]
            p_edge_4 = parity_qbts[4 * i + 3]
            ring_node = non_ring_edges[i]
            edges = [(ring_node, p_edge_1), (p_edge_1, p_edge_2), (p_edge_1, p_edge_3), (p_edge_1, p_edge_4)]
            graph.add_edges_from(edges)

    return graph


def gen_cascade_ring(num_nodes, numb_graphs):
    graph = nx.Graph()
    nodes = range(num_nodes)
    graph.add_nodes_from(nodes)
    iteration = [x for x in range(48, 78)]
    ring_edges = []
    ring_nb = []
    non_ring_edges = []
    # ring_edges = [24, 25, 29, 28, 27, 26]
    parity_qbts = [x for x in range(0, 48)]
    for idx in range(numb_graphs):
        zero, one, two, three, four = iteration[idx * 5:(idx + 1) * (5)]
        edges = [(zero, one), (one, two), (two, three),
                 (three, four), (four, zero)]
        non_ring_edges.append(one)
        non_ring_edges.append(two)
        non_ring_edges.append(three)
        non_ring_edges.append(four)
        ring_edges.append(zero)
        ring_nb.append(one)
        graph.add_edges_from(edges)
    # graph.add_edges_from(edges)
    zero, one, two, three, four, five = ring_edges
    edges = [(zero, one), (one, two), (two, three), (three, four), (four, five), (five, zero)]
    graph.add_edges_from(edges)
    # for idx in range(len(ring_nb)):
    #     graph = X_measure(ring_edges[idx], ring_nb[idx], graph)

    for i in range(len(non_ring_edges)):
        p_edge_1 = parity_qbts[2 * i]
        p_edge_2 = parity_qbts[2 * i + 1]
        ring_node = non_ring_edges[i]
        edges = [(ring_node, p_edge_1), (p_edge_1, p_edge_2)]
        graph.add_edges_from(edges)
        # graph = X_measure(ring_node, p_edge_1, graph)
    return graph



def optimized_four_code():
    edges = [(46, 47),
             (44, 45),
             (42, 43),
             (38, 39),
             (36, 37),
             (34, 35),
             (30, 31),
             (28, 29),
             (26, 27),
             (24, 25),
             (22, 23),
             (20, 21),
             (17, 41),
             (14, 15),
             (12, 13),
             (8, 9),
             (6, 7),
             (4, 5),
             (2, 3),
             (0, 33),
             (40, 47),
             (32, 38),
             (19, 22),
             (16, 17),
             (10, 30),
             (18, 22),
             (2, 7),
             (40, 43),
             (32, 34),
             (10, 26),
             (0, 40),
             (10, 25),
             (10, 33),
             (10, 40),
             (10, 19),
             (10, 18),
             (40, 41),
             (10, 17),
             (17, 32),
             (32, 33),
             (17, 40),
             (0, 32),
             (0, 23),
             (0, 22),
             (8, 10),
             (10, 14),
             (10, 11),
             (0, 2),
             (0, 1),
             (42, 44),
             (41, 44),
             (35, 36),
             (33, 36),
             (26, 29),
             (24, 29),
             (18, 20),
             (17, 20),
             (10, 13),
             (9, 13), (3, 5), (0, 5)]
    graph = nx.Graph()
    graph.add_edges_from(edges)
    return graph



def gen_encoded_five_code(qubits_per_graph, numb_graphs, concat_flag=True):
    graph = nx.Graph()
    iteration = [x for x in range(2 * 6 * (qubits_per_graph - 1), 2 * 6 * (qubits_per_graph - 1) + 6 * qubits_per_graph)]
    ring_edges = []
    ring_nb = []
    non_ring_edges = []
    parity_qbts = [x for x in range(0, 2 * 6 * (qubits_per_graph- 1))]
    for idx in range(numb_graphs):
        zero, one, two, three, four, five = iteration[idx * qubits_per_graph:(idx + 1) * (qubits_per_graph)]
        edges = [(zero, one), (one, two), (two, three), (three, four),
                 (five, four), (five, zero), (five, two)]
        non_ring_edges.append(one)
        non_ring_edges.append(two)
        non_ring_edges.append(three)
        non_ring_edges.append(four)
        non_ring_edges.append(five)
        ring_edges.append(zero)
        ring_nb.append(one)
        graph.add_edges_from(edges)

    zero, one, two, three, four, five = ring_edges
    edges = [(zero, one), (one, two), (two, three), (three, four), (four, five), (five, zero)]
    graph.add_edges_from(edges)
    if concat_flag:
        for idx in range(len(ring_nb)):
            graph = X_measure(ring_edges[idx], ring_nb[idx], graph)

        for i in range(len(non_ring_edges)):
            p_edge_1 = parity_qbts[2 * i]
            p_edge_2 = parity_qbts[2 * i + 1]
            ring_node = non_ring_edges[i]
            edges = [(ring_node, p_edge_1), (p_edge_1, p_edge_2)]
            graph.add_edges_from(edges)
            graph = X_measure(ring_node, p_edge_1, graph)
    else:
        for i in range(len(non_ring_edges)):
            p_edge_1 = parity_qbts[2 * i]
            p_edge_2 = parity_qbts[2 * i + 1]
            ring_node = non_ring_edges[i]
            edges = [(ring_node, p_edge_1), (p_edge_1, p_edge_2)]
            graph.add_edges_from(edges)
    return graph



def optimized_encoded_five_code():
    edges = [(58, 59),
(56, 57),
(52, 53),
(48, 49),
(44, 45),
(43, 47),
(40, 41),
(38, 39),
(36, 37),
(34, 35),
(32, 33),
(28, 29),
(26, 27),
(24, 25),
(22, 23),
(18, 19),
(12, 13),
(11, 17),
(6, 7),
(5, 9),
(2, 3),
(55, 58),
(50, 57),
(43, 46),
(30, 36),
(20, 26),
(15, 25),
(11, 16),
(1, 7),
(54, 58),
(5, 29),
(50, 53),
(50, 51),
(40, 43),
(42, 43),
(30, 32),
(30, 31),
(20, 22),
(20, 21),
(11, 12),
(10, 11),
(1, 2),
(0, 1),
(15, 28),
(5, 24),
(15, 48),
(15, 45),
(15, 38),
(15, 34),
(15, 55),
(15, 54),
(5, 15),
(15, 18),
(14, 15),
(8, 9),
(4, 5),
(42, 45),
(40, 49),
(22, 24),
(21, 28),
(11, 18),
(12, 15),
(1, 9),
(2, 5),
(52, 55),
(50, 58),
(32, 35),
(30, 39)]
    graph = nx.Graph()
    graph.add_edges_from(edges)
    return graph

def gen_encoded_six_code(qubits_per_graph, numb_graphs, concat_flag=True):
    graph = nx.Graph()
    iteration = [x for x in range(2 * 6 * (qubits_per_graph - 1), 2 * 6 * (qubits_per_graph - 1) + 6 * qubits_per_graph)]
    ring_edges = []
    ring_nb = []
    non_ring_edges = []
    parity_qbts = [x for x in range(0, 2 * 6 * (qubits_per_graph - 1))]

    for idx in range(numb_graphs):
        zero, one, two, three, four, five, six = iteration[idx * qubits_per_graph:(idx + 1) * (qubits_per_graph)]
        edges = [(zero, one), (one, two), (two, three), (three, four), (four, five), (five, six), (six, zero)]
        graph.add_edges_from(edges)
        non_ring_edges.append(one)
        non_ring_edges.append(two)
        non_ring_edges.append(three)
        non_ring_edges.append(four)
        non_ring_edges.append(five)
        non_ring_edges.append(six)
        ring_edges.append(zero)
        ring_nb.append(one)
        graph.add_edges_from(edges)
    zero, one, two, three, four, five = ring_edges
    edges = [(zero, one), (one, two), (two, three), (three, four), (four, five), (five, zero)]
    graph.add_edges_from(edges)
    if concat_flag:
        for idx in range(len(ring_nb)):
            graph = X_measure(ring_edges[idx], ring_nb[idx], graph)

        for i in range(len(non_ring_edges)):
            p_edge_1 = parity_qbts[2 * i]
            p_edge_2 = parity_qbts[2 * i + 1]
            ring_node = non_ring_edges[i]
            edges = [(ring_node, p_edge_1), (p_edge_1, p_edge_2)]
            graph.add_edges_from(edges)
            graph = X_measure(ring_node, p_edge_1, graph)
    else:
        for i in range(len(non_ring_edges)):
            p_edge_1 = parity_qbts[2 * i]
            p_edge_2 = parity_qbts[2 * i + 1]
            ring_node = non_ring_edges[i]
            edges = [(ring_node, p_edge_1), (p_edge_1, p_edge_2)]
            graph.add_edges_from(edges)
    return graph



def optimized_six_code():
    edges = [(68, 69), (66, 67), (64, 65), (62, 63),
             (60, 70), (56, 59), (54, 55), (52, 53), (50, 51), (46, 47),
             (44, 45), (42, 43), (40, 41), (39, 57), (36, 37), (34, 35), (32, 33), (30, 31), (28, 29), (26, 27), (24, 25), (20, 21), (18, 19), (16, 17), (15, 23), (13, 22), (10, 58), (8, 9), (6, 7), (4, 5), (2, 3), (0, 71), (60, 61), (49, 55), (38, 39), (14, 15), (12, 22), (0, 11), (0, 1), (62, 65), (43, 45), (30, 59), (26, 32), (25, 34), (2, 5), (49, 51), (48, 49), (14, 21), (62, 70), (61, 65), (43, 57), (38, 44), (27, 29), (35, 56), (25, 31), (14, 16), (3, 11), (1, 4), (46, 57), (37, 44), (36, 57), (53, 55), (66, 70), (36, 40), (39, 40), (6, 11), (68, 70), (11, 58), (9, 70), (22, 57), (9, 11), (7, 70), (6, 70), (70, 71), (0, 70), (34, 57), (31, 57), (22, 59), (13, 59), (1, 60), (10, 60), (11, 35), (11, 30), (57, 58), (11, 22), (11, 13), (11, 60), (10, 57), (0, 58), (0, 57), (22, 71), (22, 70), (22, 60), (18, 22), (30, 32), (26, 35), (26, 34), (22, 23), (19, 20), (15, 22), (25, 31), (25, 30), (25, 29), (55, 56), (54, 56), (49, 56), (56, 57), (56, 58), (13, 19), (13, 18), (13, 16), (49, 59)]
    graph = nx.Graph()
    graph.add_edges_from(edges)
    return graph



def gen_encoded_seven_code(qubits_per_graph, numb_graphs, concat_flag=True):
    graph = nx.Graph()
    iteration = [x for x in range(2 * 6 * (qubits_per_graph - 1), 2 * 6 * (qubits_per_graph - 1) + 6 * qubits_per_graph)]
    ring_edges = []
    ring_nb = []
    non_ring_edges = []
    parity_qbts = [x for x in range(0, 2 * 6 * (qubits_per_graph - 1))]

    for idx in range(numb_graphs):
        zero, one, two, three, four, five, six, seven = iteration[idx * qubits_per_graph:(idx + 1) * (qubits_per_graph)]
        edges = [(zero, one), (zero, two), (zero, three), (one, six),
                 (one, five), (five, three), (five, four),
                 (four, three), (four, two), (two ,seven), (seven, three),
                 (seven, six), (six, three)]
        graph.add_edges_from(edges)
        non_ring_edges.append(one)
        non_ring_edges.append(two)
        non_ring_edges.append(three)
        non_ring_edges.append(four)
        non_ring_edges.append(five)
        non_ring_edges.append(six)
        non_ring_edges.append(seven)
        ring_edges.append(zero)
        ring_nb.append(one)
        graph.add_edges_from(edges)
    zero, one, two, three, four, five = ring_edges
    edges = [(zero, one), (one, two), (two, three), (three, four), (four, five), (five, zero)]
    graph.add_edges_from(edges)
    if concat_flag:
        for idx in range(len(ring_nb)):
            graph = X_measure(ring_edges[idx], ring_nb[idx], graph)

        for i in range(len(non_ring_edges)):
            p_edge_1 = parity_qbts[2 * i]
            p_edge_2 = parity_qbts[2 * i + 1]
            ring_node = non_ring_edges[i]
            edges = [(ring_node, p_edge_1), (p_edge_1, p_edge_2)]
            graph.add_edges_from(edges)
            graph = X_measure(ring_node, p_edge_1, graph)
    else:
        for i in range(len(non_ring_edges)):
            p_edge_1 = parity_qbts[2 * i]
            p_edge_2 = parity_qbts[2 * i + 1]
            ring_node = non_ring_edges[i]
            edges = [(ring_node, p_edge_1), (p_edge_1, p_edge_2)]
            graph.add_edges_from(edges)
    return graph


def gen_encoded_eigth_code(qubits_per_graph, numb_graphs, concat_flag=True):
    graph = nx.Graph()
    iteration = [x for x in range(2 * 6 * (qubits_per_graph - 1), 2 * 6 * (qubits_per_graph - 1) + 6 * qubits_per_graph)]
    ring_edges = []
    ring_nb = []
    non_ring_edges = []
    parity_qbts = [x for x in range(0, 2 * 6 * (qubits_per_graph - 1))]

    for idx in range(numb_graphs):
        zero, one, two, three, four, five, six, seven, eigth = iteration[idx * qubits_per_graph:(idx + 1) * (qubits_per_graph)]
        edges = [(zero, three), (zero, two), (zero, four), (zero, five), (three, seven),
                 (three, two), (two, six),
                 (six, four), (six, eigth), (four, seven), (seven, one),
                 (one, five), (five, eigth)]
        graph.add_edges_from(edges)
        non_ring_edges.append(one)
        non_ring_edges.append(two)
        non_ring_edges.append(three)
        non_ring_edges.append(four)
        non_ring_edges.append(five)
        non_ring_edges.append(six)
        non_ring_edges.append(seven)
        non_ring_edges.append(eigth)
        ring_edges.append(zero)
        ring_nb.append(three)
        graph.add_edges_from(edges)
    zero, one, two, three, four, five = ring_edges
    edges = [(zero, one), (one, two), (two, three), (three, four), (four, five), (five, zero)]
    graph.add_edges_from(edges)
    if concat_flag:
        for idx in range(len(ring_nb)):
            graph = X_measure(ring_edges[idx], ring_nb[idx], graph)

        for i in range(len(non_ring_edges)):
            p_edge_1 = parity_qbts[2 * i]
            p_edge_2 = parity_qbts[2 * i + 1]
            ring_node = non_ring_edges[i]
            edges = [(ring_node, p_edge_1), (p_edge_1, p_edge_2)]
            graph.add_edges_from(edges)
            graph = X_measure(ring_node, p_edge_1, graph)
    else:
        for i in range(len(non_ring_edges)):
            p_edge_1 = parity_qbts[2 * i]
            p_edge_2 = parity_qbts[2 * i + 1]
            ring_node = non_ring_edges[i]
            edges = [(ring_node, p_edge_1), (p_edge_1, p_edge_2)]
            graph.add_edges_from(edges)
    return graph



def gen_encoded_nine_code(qubits_per_graph, numb_graphs, concat_flag=True):
    graph = nx.Graph()
    iteration = [x for x in
                 range(2 * 6 * (qubits_per_graph - 1), 2 * 6 * (qubits_per_graph - 1) + 6 * qubits_per_graph)]
    ring_edges = []
    ring_nb = []
    non_ring_edges = []
    parity_qbts = [x for x in range(0, 2 * 6 * (qubits_per_graph - 1))]
    for idx in range(numb_graphs):
        zero, one, two, three, four, five, six, seven, eigth, nine = iteration[idx * qubits_per_graph:(idx + 1) * (qubits_per_graph)]
        edges = [(zero, seven), (zero, three), (zero, nine), (three, nine), (three, eigth), (seven, four), (seven, one), (one, two),
                        (two, nine), (two, six), (six, five), (six, eigth), (eigth, four), (four, nine)]
        graph.add_edges_from(edges)
        non_ring_edges.append(one)
        non_ring_edges.append(two)
        non_ring_edges.append(three)
        non_ring_edges.append(four)
        non_ring_edges.append(five)
        non_ring_edges.append(six)
        non_ring_edges.append(seven)
        non_ring_edges.append(eigth)
        non_ring_edges.append(nine)
        ring_edges.append(zero)
        ring_nb.append(one)
        graph.add_edges_from(edges)
    zero, one, two, three, four, five = ring_edges
    edges = [(zero, one), (one, two), (two, three), (three, four), (four, five), (five, zero)]
    graph.add_edges_from(edges)
    if concat_flag:
        for idx in range(len(ring_nb)):
            graph = X_measure(ring_edges[idx], ring_nb[idx], graph)

        for i in range(len(non_ring_edges)):
            p_edge_1 = parity_qbts[2 * i]
            p_edge_2 = parity_qbts[2 * i + 1]
            ring_node = non_ring_edges[i]
            edges = [(ring_node, p_edge_1), (p_edge_1, p_edge_2)]
            graph.add_edges_from(edges)
            graph = X_measure(ring_node, p_edge_1, graph)
    else:
        for i in range(len(non_ring_edges)):
            p_edge_1 = parity_qbts[2 * i]
            p_edge_2 = parity_qbts[2 * i + 1]
            ring_node = non_ring_edges[i]
            edges = [(ring_node, p_edge_1), (p_edge_1, p_edge_2)]
            graph.add_edges_from(edges)
    return graph


def gen_encoded_ten_code(qubits_per_graph, numb_graphs, concat_flag=True):
    graph = nx.Graph()
    iteration = [x for x in range(2 * 6 * (qubits_per_graph - 1), 2 * 6 * (qubits_per_graph - 1) + 6 * qubits_per_graph)]
    ring_edges = []
    ring_nb = []
    non_ring_edges = []
    parity_qbts = [x for x in range(0, 2 * 6 * (qubits_per_graph - 1))]
    for idx in range(numb_graphs):
        zero, one, two, three, four, five, six, seven, eigth, nine, ten = iteration[idx * qubits_per_graph:(idx + 1) * (
            qubits_per_graph)]
        edges = [(zero, four), (zero, nine), (zero, eigth),
         (four, seven), (four, nine),
         (seven, five), (seven, one),
         (five, ten), (ten, six),
         (ten, eigth), (ten, nine),
         (one, six),
         (six, two), (two, eigth),
         (two, three), (three, eigth),
         (three, nine)]
        graph.add_edges_from(edges)
        non_ring_edges.append(one)
        non_ring_edges.append(two)
        non_ring_edges.append(three)
        non_ring_edges.append(four)
        non_ring_edges.append(five)
        non_ring_edges.append(six)
        non_ring_edges.append(seven)
        non_ring_edges.append(eigth)
        non_ring_edges.append(nine)
        non_ring_edges.append(ten)
        ring_edges.append(zero)
        ring_nb.append(four)
        graph.add_edges_from(edges)

    zero, one, two, three, four, five = ring_edges
    edges = [(zero, one), (one, two), (two, three), (three, four), (four, five), (five, zero)]
    graph.add_edges_from(edges)
    if concat_flag:
        for idx in range(len(ring_nb)):
            graph = X_measure(ring_edges[idx], ring_nb[idx], graph)

        for i in range(len(non_ring_edges)):
            p_edge_1 = parity_qbts[2 * i]
            p_edge_2 = parity_qbts[2 * i + 1]
            ring_node = non_ring_edges[i]
            edges = [(ring_node, p_edge_1), (p_edge_1, p_edge_2)]
            graph.add_edges_from(edges)
            graph = X_measure(ring_node, p_edge_1, graph)
    else:
        for i in range(len(non_ring_edges)):
            p_edge_1 = parity_qbts[2 * i]
            p_edge_2 = parity_qbts[2 * i + 1]
            ring_node = non_ring_edges[i]
            edges = [(ring_node, p_edge_1), (p_edge_1, p_edge_2)]
            graph.add_edges_from(edges)
    return graph





def general_parity_codes_concat(qubits_per_graph, numb_graphs=6, concat_flag=True):
    graph = nx.Graph()
    iteration = [x for x in range(2 * 6 * (qubits_per_graph - 1), 2 * 6 * (qubits_per_graph - 1) + 6 * qubits_per_graph)]
    ring_edges = []
    ring_nb = []
    non_ring_edges = []
    parity_qbts = [x for x in range(0, 2 * 6 * (qubits_per_graph - 1))]
    for idx in range(numb_graphs):
        qbts = iteration[idx * qubits_per_graph:(idx + 1) * (
            qubits_per_graph)]
        zero = qbts[0]
        edges = []
        for qbt in qbts[1:]:
            edges.append((zero, qbt))
            non_ring_edges.append(qbt)
        graph.add_edges_from(edges)
        ring_edges.append(zero)
        ring_nb.append(qbts[1])

    zero, one, two, three, four, five = ring_edges
    edges = [(zero, one), (one, two), (two, three), (three, four), (four, five), (five, zero)]
    graph.add_edges_from(edges)
    if concat_flag:
        for idx in range(len(ring_nb)):
            graph = X_measure(ring_edges[idx], ring_nb[idx], graph)

        for i in range(len(non_ring_edges)):
            p_edge_1 = parity_qbts[2 * i]
            p_edge_2 = parity_qbts[2 * i + 1]
            ring_node = non_ring_edges[i]
            edges = [(ring_node, p_edge_1), (p_edge_1, p_edge_2)]
            graph.add_edges_from(edges)
            graph = X_measure(ring_node, p_edge_1, graph)
    else:
        for i in range(len(non_ring_edges)):
            p_edge_1 = parity_qbts[2 * i]
            p_edge_2 = parity_qbts[2 * i + 1]
            ring_node = non_ring_edges[i]
            edges = [(ring_node, p_edge_1), (p_edge_1, p_edge_2)]
            graph.add_edges_from(edges)
    return graph


def optimized_graphs(code_size):
    path_graphs_dict = r"C:\Users\Admin\Desktop\HeraldedResourceStateGen\graphsdata\hemant_opt_graphs.json"
    f = open(path_graphs_dict)
    graphs_dict = json.load(f)
    key = str(code_size - 4)
    edges = graphs_dict[key]
    graph = nx.Graph()
    graph.add_edges_from(edges)
    return graph


def optimized_four_parity_code():
    path_graphs_dict = r"C:\Users\Admin\Downloads\193238_cz_data.json"
    f = open(path_graphs_dict)
    graphs_dict = json.load(f)
    edges = graphs_dict["1"]
    graph = nx.Graph()
    graph.add_edges_from(edges)
    return graph


def concatenated_cube():
    graph = nx.Graph()
    iteration = [x for x in range(8 * 8)]
    outer_cube = []
    outer_cube_nb = []
    for idx in range(8):
        zero, one, two, three, four, five, six, seven = iteration[idx * 8:(idx + 1) * (
            8)]
        edges = [(zero, one), (zero, six), (five, six), (six, seven), (seven, one), (seven, four),
                 (five, three), (five, four), (three, zero), (three, two), (two, one), (two, four)]
        outer_cube.append(zero)
        outer_cube_nb.append(one)
        graph.add_edges_from(edges)

    zero, one, two, three, four, five, six, seven = outer_cube
    edges = [(zero, one), (zero, six), (five, six), (six, seven), (seven, one), (seven, four),
             (five, three), (five, four), (three, zero), (three, two), (two, one), (two, four)]
    graph.add_edges_from(edges)
    for idx in range(len(outer_cube_nb)):
        graph = X_measure(outer_cube[idx], outer_cube_nb[idx], graph)

    return graph

if __name__ == '__main__':
    import json

    # path_graphs_dict = r"C:\Users\Admin\Desktop\HeraldedResourceStateGen\graphsdata\hemant_opt_graphs.json"
    # f = open(path_graphs_dict)
    # graphs_dict = json.load(f)

    # for key in graphs_dict.keys():
    #     edges = graphs_dict[key]
    #     print(key, len(edges))

    # graph = concatenated_cube()  # gen_concat_ring(24, 6)
    # # graph = local_complementation_on_graphs(graph, 8)
    # nx.draw(graph)
    # print(len(list(graph.edges)))
    # print(list(graph.edges))
    # plt.show()
    graph = general_parity_codes_concat(6, 6, concat_flag=True)
    nx.draw(graph, with_labels=True)
    plt.show()
