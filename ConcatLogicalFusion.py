import numpy as np
import copy
from itertools import product


def generate_layer_comb(numb_layers):
    combinations = product('23456789', repeat=numb_layers)
    # combinations = product('6789', repeat=numb_layers)
    return list(combinations)

def log_fusion_probs(p_s, p_f_x, p_f_z, p_f_y, p_l, etas, losses_sing,  dicts_in, f_type):
    log_fusion_succ_prob = 0

    for trajectories in dicts_in[f_type]:
        key, value = trajectories
        fusions, ffx, ffy, ffz, losses = key
        prob_f = (p_s) ** len(fusions) * (p_f_x ** len(ffx)) * (p_f_y ** len(ffy)) * (p_f_z ** len(ffz)) * (p_l ** len(losses))  #    + len(ffy) + len(ffz)) * (1 - eta ** (1/pfail)) ** len(losses)
        pp = 0
        for item in value:
            nx = item['x']
            ny = item['y']
            nz = item['z']
            nxf = item['xf']
            nyf = item['yf']
            nzf = item['zf']
            pp += (etas[0] ** nx) * (etas[1] ** nz) * (etas[2] ** ny) * (losses_sing[0] ** nxf) * (losses_sing[1] ** nzf) * (losses_sing[2] ** nyf)

        traj_prob = prob_f * (pp ** 2)
        log_fusion_succ_prob += traj_prob
    return log_fusion_succ_prob



def get_prob(list_of_probs, etas, losses):
    tot_prob = 0
    for traj in list_of_probs:
        x, xf, z, zf, y, yf = traj
        prob = (etas[0] ** x) * (losses[0] ** xf) * (etas[1] ** z) * (losses[1] ** zf) * (etas[2] ** y) * (losses[2] ** yf)
        tot_prob += prob
    return tot_prob




def concat_recursion(fusion_dict_list, pauli_dict_list, eta, numb_concat_layers, spin_photon_gate_flag=False):
    # p_s, p_f_x, p_f_z, p_f_y, p_l = (1 / 2) * (eta ** 2), (1 / 2) * (eta ** 2), \
    #                                 (1 / 2) * (eta ** 2), (1 / 2) * (eta ** 2), 1 - eta ** 2
    if spin_photon_gate_flag:
        p_s, p_f_x, p_f_z, p_f_y, p_l = eta ** 2, 0, 0, 0, 1 - eta ** 2
    else:
        p_s, p_f_x, p_f_z, p_f_y, p_l = (1 / 2) * (eta ** 2), (1 / 2) * (eta ** 2), \
                                          (1 / 2) * (eta ** 2), (1 / 2) * (eta ** 2), 1 - eta ** 2
    init_etas = [eta for _ in range(3)]
    init_losses = [1 - eta for _ in range(3)]
    bases = ["x", "z", "y"]

    fusion_dict = fusion_dict_list[0]
    pauli_dict = pauli_dict_list[0]
    p_s_next_layer = log_fusion_probs(p_s, p_f_x, p_f_z, p_f_y, p_l, init_etas, init_losses, fusion_dict, "fusion")
    p_f_x_next_layer = log_fusion_probs(p_s, p_f_x, p_f_z, p_f_y, p_l, init_etas, init_losses, fusion_dict, "xrec")
    p_f_z_next_layer = log_fusion_probs(p_s, p_f_x, p_f_z, p_f_y, p_l, init_etas, init_losses, fusion_dict, "zrec")
    p_f_y_next_layer = 0
    p_l_next_layer = 1 - p_s_next_layer - p_f_x_next_layer - p_f_y_next_layer - p_f_z_next_layer
    etas_next_layer = [get_prob(pauli_dict[b], init_etas, init_losses) for b in bases]
    losses_next_layer = [1 - etas_next_layer[i] for i in range(len(etas_next_layer))]

    etas = copy.deepcopy(etas_next_layer)
    losses = copy.deepcopy(losses_next_layer)
    p_s = copy.deepcopy(p_s_next_layer)
    p_f_x = copy.deepcopy(p_f_x_next_layer)
    p_f_z = copy.deepcopy(p_f_z_next_layer)
    p_f_y = 0
    p_l = copy.deepcopy(p_l_next_layer)
    for l_idx in range(numb_concat_layers - 1):
        fusion_dict = fusion_dict_list[l_idx + 1]
        pauli_dict = pauli_dict_list[l_idx + 1]
        p_s_next_layer = log_fusion_probs(p_s, p_f_x, p_f_z, p_f_y, p_l, etas, losses, fusion_dict, "fusion")
        p_f_x_next_layer = log_fusion_probs(p_s, p_f_x, p_f_z, p_f_y, p_l, etas, losses, fusion_dict, "xrec")
        p_f_z_next_layer = log_fusion_probs(p_s, p_f_x, p_f_z, p_f_y, p_l, etas, losses, fusion_dict, "zrec")
        p_f_y_next_layer = 0
        p_l_next_layer = 1 - p_s_next_layer - p_f_x_next_layer - p_f_y_next_layer - p_f_z_next_layer
        etas_next_layer = [get_prob(pauli_dict[b], etas, losses) for b in bases]
        losses_next_layer = [1 - etas_next_layer[i] for i in range(len(etas_next_layer))]

        etas = copy.deepcopy(etas_next_layer)
        losses = copy.deepcopy(losses_next_layer)
        p_s = copy.deepcopy(p_s_next_layer)
        p_f_x = copy.deepcopy(p_f_x_next_layer)
        p_f_z = copy.deepcopy(p_f_z_next_layer)
        p_l = copy.deepcopy(p_l_next_layer)

    erasure = 1 - (2 * p_s + p_f_x + p_f_z) / 2
    return erasure, p_s

def bisection_search(interval, func, depth=10, y0=None, y1=None):
    k0, k1 = interval[0], interval[1]
    mid = 0.5 * (k0 + k1)
    if depth == 0:
        return (k0 + k1)/2
    else:
        if y0 is None:
            y0 = func(k0)
            y1 = func(k1)
        ymid = func(mid)
        if y0 * ymid < 0:
            return bisection_search([k0, mid], func, depth=depth-1, y0=y0, y1=ymid)
        elif ymid * y1 < 0:
            return bisection_search([mid, k1], func, depth=depth-1, y0=ymid, y1=y1)
        else:
            raise ValueError('No root in interval')


def fusion_threshold_from_dict(fusion_dict_list, pauli_dict_list, numb_concat_layers, spin_photon_gate_flag=False, pthresh=0.88):
    """
    Get the FBQC threshold for the 6-ring approach to generating the Raussendorf lattice in http://arxiv.org/abs/2101.09310
    :param dict:
    :param pf:
    :param take_min:
    :param pthresh:
    :return:
    """

    def func_to_optimise(eta):
        erasure, _ = concat_recursion(fusion_dict_list, pauli_dict_list, eta, numb_concat_layers, spin_photon_gate_flag)  # prob_from_dict(eta, pf, dict, decoder_type=decoder_type)
        t = (1 - erasure) - pthresh
        return t

    try:
        threshold = bisection_search([0.55, 0.95], func_to_optimise)
    except ValueError:
        threshold = 1
    return threshold





def find_opt_layer_config(fusion_dicts_in, pauli_dicts_in, numb_layers, spin_photon_gate_flag=False, max_size=500):
    list_of_combs = generate_layer_comb(numb_layers)
    best_thres = 1
    best_combination = 0
    x_scatter = []
    y_scatter = []
    print(len(list_of_combs))
    # list_of_combs = [('6', '6', '6', '6', '6', '9'), ('6', '6', '6', '6', '6', '8')]
    for comb in list_of_combs:
        fusion_dicts = []
        pauli_dicts = []
        size = 1
        for key in comb:
            size *= (int(key) + 1)
            fusion_dicts.append(fusion_dicts_in[str(int(key) + 1)])  # Goes from 2-9, when it should go from 3-10, thus the "+1".
            pauli_dicts.append(pauli_dicts_in[str(int(key) - 2)])
        if size < max_size:
            threshold_val = fusion_threshold_from_dict(fusion_dicts, pauli_dicts, numb_layers, spin_photon_gate_flag)
            x_scatter.append(size)
            y_scatter.append(threshold_val)
            if threshold_val < best_thres:
                best_thres = threshold_val
                best_combination = comb
                print("New best threshold: ", best_thres, ", with layers: ", best_combination, numb_layers)

    plt.scatter(x_scatter, y_scatter)
    plt.show()
    return x_scatter, y_scatter


def threshold_fixed_input(fusion_dicts_in, pauli_dicts_in, numb_layers, comb, spin_photon_gate_flag=False):
    fusion_dicts = []
    pauli_dicts = []
    size = 1
    for key in comb:
        size *= (int(key) + 1)
        fusion_dicts.append(fusion_dicts_in[str(int(key) + 1)])
        pauli_dicts.append(pauli_dicts_in[str(int(key) - 2)])
    threshold_val = fusion_threshold_from_dict(fusion_dicts, pauli_dicts, numb_layers, spin_photon_gate_flag)

    return threshold_val


def parse_loss_size_data(losses, sizes):
    loss_tolerance_dict = {}
    final_losses = []
    final_sizes = []
    for idx in range(len(losses)):
        key = sizes[idx]
        loss = losses[idx]
        if key in loss_tolerance_dict.keys():
            loss_tolerance_dict[key].append(loss)
        else:
            loss_tolerance_dict[key] = [loss]
    size_list = [key for key in loss_tolerance_dict.keys()]
    size_list.sort()
    best_trans = 1
    for size in size_list:
        min_trans = min(loss_tolerance_dict[size])
        if min_trans < best_trans:
            final_losses.append(1 - min_trans)
            final_sizes.append(size * 6)
            best_trans = min_trans
    # final_sizes.append((7 ** 4) * 9 * 6)
    # final_losses.append(1 - 0.60529785)
    return final_losses, final_sizes


def psi_quantum_numbers():
    sizes_la = [24, 36, 48, 60, 72, 90, 108, 120, 144, 168]
    loss_la = [0.057, 0.068, 0.092, 0.105, 0.107, 0.115, 0.119, 0.125, 0.133, 0.14]

    sizes_sb = [36, 72, 168, 600, 4200, 60000, 780000, 9600000]
    loss_sb = [0.051, 0.077, 0.113, 0.167, 0.208, 0.236, 0.249, 0.256]

    sizes_ba = [12, 24, 72, 168]
    loss_gb = [0.026, 0.075, 0.139, 0.174]
    return sizes_la, loss_la, sizes_sb, loss_sb, sizes_ba, loss_gb


def rings_quantum_numbers():
    sizes = [4 * 6, (4 ** 2) * 6, (4 ** 3) * 6, (4 ** 4) * 6, (4 ** 5) * 6, (4 ** 6) * 6, (4 ** 7) * 6]
    losses = [1 - 0.9545150501672242, 1 - 0.8876254180602007, 1 - 0.8234113712374582, 1 - 0.7752508361204014, 1 - 0.7364548494983277, 1 - 0.7083612040133779, 1 - 0.6869565217391305]
    return sizes, losses

if __name__ == '__main__':
    import json
    import matplotlib.pyplot as plt

    # path_fusion_dict = r"C:\Users\Admin\Desktop\HeraldedResourceStateGen\graphsdata\save_fusion_dicts.json"
    # path_pauli_dict = r"C:\Users\Admin\Desktop\HeraldedResourceStateGen\graphsdata\pauli_tree.json"

    path_fusion_dict = r"C:\Users\Admin\Desktop\HeraldedResourceStateGen\graphsdata\PauliAndFusionMoreGraphs\save_fusion_dicts.json"
    path_pauli_dict = r"C:\Users\Admin\Desktop\HeraldedResourceStateGen\graphsdata\PauliAndFusionMoreGraphs\pauli_tree.json"


    # Loop through all combinations of three layers ?
    f = open(path_fusion_dict)
    fusion_dict = json.load(f)

    f = open(path_pauli_dict)
    pauli_dict = json.load(f)
    print(fusion_dict.keys())

    final_x_scatter = []
    final_y_scatter = []
    for numb_layers in range(1, 5):
        x_scatter, y_scatter = find_opt_layer_config(fusion_dict, pauli_dict, numb_layers, spin_photon_gate_flag=False, max_size=10**9)
        for idx in range(len(x_scatter)):
            final_x_scatter.append(x_scatter[idx])
            final_y_scatter.append(y_scatter[idx])
    plt.scatter(final_x_scatter, final_y_scatter, color="red")
    plt.show()

    parse_y, parse_x = parse_loss_size_data(final_y_scatter, final_x_scatter)
    plt.plot(parse_x, parse_y, "-o", color="blue", label="Us")

    sizes_la, loss_la, sizes_sb, loss_sb, sizes_ba, loss_gb = psi_quantum_numbers()
    plt.plot(sizes_la, loss_la, "-o", color="red", label="Psi LA")
    plt.plot(sizes_sb, loss_sb, "-o", color="purple", label="Psi SB")
    plt.plot(sizes_ba, loss_gb, "-o", color="orange", label="Psi GB")

    ring_size, ring_loss = rings_quantum_numbers()
    plt.plot(ring_size, ring_loss, "-o", color="pink", label="Rings")

    best_psi = [max(loss_sb) for _ in range(len(sizes_sb))]
    plt.plot(sizes_sb, best_psi, "--", color="black")
    plt.legend()
    plt.xscale("log")
    plt.ylabel("Loss threshold")
    plt.xlabel("# Qubits in resource state")
    plt.show()

    # N_layers = [1, 2, 3, 4, 5]
    # loss_thresholds = []
    # for N in N_layers:
    #     comb = ("5" for _ in range(N))
    #     threshold_val = threshold_fixed_input(fusion_dict, pauli_dict, N, comb, spin_photon_gate_flag=False)
    #     loss_thresholds.append(1 - threshold_val)
    # plt.plot(N_layers, loss_thresholds, "-o", color="black")
    # plt.show()
    '''
    fusion_dict_key = "4"
    pauli_dict_key = "0"
    numb_concat_layers = 2
    fusion_dicts = [fusion_dict["7"], fusion_dict["10"]]
    pauli_dicts = [pauli_dict["3"], pauli_dict["6"]]

    transmissions = np.linspace(0.5, 1, 25)
    fusion_success_list = []
    erasure_list = []

    for eta in transmissions:
        erasure, f_succ = concat_recursion(fusion_dicts, pauli_dicts, eta, numb_concat_layers)
        fusion_success_list.append(f_succ)
        erasure_list.append(erasure)
    # plt.plot(transmissions, fusion_success_list, color="black")
    plt.plot(transmissions, fusion_success_list, color="black")
    plt.plot(transmissions, transmissions, "--", color="red")
    plt.show()
    '''
    # threshold_val = fusion_threshold_from_dict(fusion_dicts, pauli_dicts, numb_concat_layers)
    # print("threshold: ", threshold_val)