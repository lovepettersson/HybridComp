import matplotlib.pyplot as plt
import numpy as np
from MiscPlotsPaper import *
from RepetionCodes import get_threshold_given_eps_in_optimize_func
from ErrorsFromHigherOrderPhotonComponents import effective_loss_rate, error_rate


def get_eta_thresholds(pair_gens):
    eta_id = [0.8, 0.9]
    etas_s = 0.99
    eps_spin = 0.01
    rep_codes_sizes = [3, 4, 5, 6, 7, 8]
    epses = [0.01 / 100, 0.05 / 100]
    list_thresholds = []
    for eta in eta_id:
        for eps in epses:
            thresholds = {3 : [], 4:[], 5:[], 6:[], 7:[], 8:[]}
            for p in pair_gens:
                

def mux_plot_optimize_p(t_bells, linestyles, p_th=1):
    colors = ["red", "blue"]
    fig, ax = plt.subplots(figsize=(9, 6))

    # The success probability of heralding a resource state given spin error which has to be multiplexed to average to one
    succ_lists = [[0.7866335111212663, 0.7261559912925004, 0.6703280704865696, 0.6187922808189714, 0.5712186370521632, 0.527302523689994],
                  [0.7866335111212663, 0.7261559912925004, 0.6703280704865696, 0.6187922808189714, 0.5712186370521632, 0.527302523689994]]

    # The tolerable transmission rates under 1% spin error and with eps = 0.05% and eps = 0.01% spin-photon error
    # etas_codes = [[0.9996984924623116, 0.9903517587939699, 0.987537688442211, 0.9874371859296482, 0.988140703517588,
    #                0.9892462311557789],
    #               [1, 0.9960804020100502, 0.9936683417085427, 0.9936683417085427, 0.9944723618090452, 0.9954773869346734]]
    eps_spin = 0.01
    rep_codes_sizes = [3, 4, 5, 6, 7, 8]
    sizes_code = [3 * 6, 4 * 6, 5 * 6, 6 * 6, 7 * 6, 8 * 6]
    epses = [0.01 / 100, 0.05 / 100]
    eta_spin_gates = [0.8, 0.9]
    handles = []
    dBS = np.linspace(0.0001, 0.005, 50)
    dBS = [0.0001, 0.0005]
    pair_gen_probs = np.linspace(0.03, 0.1, 10)
    pair_gen_probs = [0.03, 0.06, 0.09]
    # dBS = np.linspace(0.0001, 0.008, 400)
    for idx_eps, eps in enumerate(epses):
        succ_list = succ_lists[idx_eps]
        t_bell = t_bells[idx_eps]
        for ix_eta, eta_spin_gate_in in enumerate(eta_spin_gates):
            eta_spin_gate = eta_spin_gate_in * eta_spin_gate_in
            code_counts = {"0":0, "1":0, "2":0, "3":0, "4":0}
            multiplexing = []
            switch_loss_x_axis = []
            size_spat_mux = []
            for db in dBS:
                switch_loss = 1 - 10 ** (-db / 10)
                multiplex_one_instance = []
                spat_mux_one_instance = []
                best_idxs = []
                print("At DB: ", db)
                for ix_r, rep_code_size in enumerate(rep_codes_sizes):
                    number_of_photons = sizes_code[ix_r]
                    best_multiplex = 10 ** 12
                    best_spat_mux = 0
                    for pair_s in pair_gen_probs:
                        for k in range(2, 8):
                            temporal_size = k
                            for Ntemp in range(2 ** (k - 1), 2 ** k):
                                if Ntemp > 211:
                                    continue
                                else:
                                    temp_numb_input = Ntemp
                                    p_succ = ((1 - (1 - pair_s * eta_spin_gate) ** temp_numb_input) ** number_of_photons) # The total success probability of generating N photons with "temp_numb_input" temporal multiplexing attempts
                                    temp_max_delay = Ntemp - 1
                                    loss = delay_loss(temp_max_delay, t_bell)
                                    # Loss induced from delay lines and treating two photon clicks as photon loss
                                    # if idx_eps == 1 and ix_eta == 0:
                                    #     loss = delay_loss(temp_max_delay, t_bell) + 0.0022129873919999998
                                    # elif idx_eps == 1 and ix_eta == 1:
                                    #     loss = delay_loss(temp_max_delay, t_bell) + 0.001478206422
                                    # elif idx_eps == 0 and ix_eta == 1:
                                    #     loss = delay_loss(temp_max_delay, t_bell) + 0.0003695516055
                                    # else:
                                    #     loss = delay_loss(temp_max_delay, t_bell) + 0.0005532468479999999

                                    tot_prob = p_succ * (succ_list[ix_r]) # Multiply success probability with the error detection heralding probability

                                    # Find the binary size of a log tree for the spatial mux above
                                    # Want on average to get out 1 resource state, that is we need N*p = 1, giving N = 1 / p.
                                    # However, we use logarithmic trees, thus we find the smallest sized tree that satisfy this

                                    binary_sized_multi = -1
                                    numb_spin_inputs = 0
                                    for k_1 in range(1, 16):
                                        binary_size = 2 ** k_1
                                        mux_prob = binary_size * tot_prob
                                        if mux_prob >= 1:
                                            binary_sized_multi = k_1
                                            numb_spin_inputs = (2 ** k_1) / mux_prob
                                            break

                                    if binary_sized_multi > -1:
                                        multiplexing_loss = 1 - ((1 - switch_loss) ** (binary_sized_multi))  # Loss from the spatial multiplexing
                                        total_loss_with_multi = 1 - (1 - loss) * (1 - multiplexing_loss) * ((1 - switch_loss) ** (temporal_size + 1))  # total loss when doing spatial mux and temporal mux
                                        loss_from_detection = effective_loss_rate(pair_s, eta_spin_gate, 1 - total_loss_with_multi)
                                        total_loss_with_multi += loss_from_detection
                                        total_loss_with_photon_as_erase = (1 - p_succ * (1 - loss) * ((1 - switch_loss) ** (temporal_size + 1)))  # Total loss when when only doing temporal mux
                                        loss_from_detection = effective_loss_rate(pair_s, eta_spin_gate, 1 - total_loss_with_photon_as_erase)
                                        total_loss_with_photon_as_erase += loss_from_detection

                                        error_rate_spat = error_rate(pair_s, eta_spin_gate, 1 - total_loss_with_multi)
                                        error_rate_erase = error_rate(pair_s, eta_spin_gate, 1 - total_loss_with_photon_as_erase)

                                        spin_photon_error = (eps + error_rate_spat) * (2 / 3)
                                        threshold_loss_spat = get_threshold_given_eps_in_optimize_func(eps_spin, rep_code_size ,spin_err=spin_photon_error)

                                        spin_photon_error = (eps + error_rate_erase) * (2 / 3)
                                        threshold_loss_erase = get_threshold_given_eps_in_optimize_func(eps_spin, rep_code_size, spin_err=spin_photon_error)

                                        if total_loss_with_photon_as_erase < threshold_loss_erase:
                                            average_numb = 2 * number_of_photons  # Times two due to the spin encoding
                                            if average_numb < best_multiplex:
                                                best_multiplex = average_numb
                                        if total_loss_with_multi < threshold_loss_spat:
                                            average_numb = 2 * number_of_photons * numb_spin_inputs # (2 ** binary_sized_multi)  # Times two due to the spin encoding and the overhead from spatial mux
                                            if average_numb < best_multiplex:
                                                best_multiplex = average_numb
                                                best_spat_mux = numb_spin_inputs
                    if best_multiplex < 10 ** 6:
                        multiplex_one_instance.append(best_multiplex)
                        best_idxs.append(ix_r)
                        spat_mux_one_instance.append(best_spat_mux)

                if len(multiplex_one_instance) > 0:
                    multiplexing.append(min(multiplex_one_instance))
                    switch_loss_x_axis.append(db * 1000)
                    code_idx_idx = multiplex_one_instance.index(min(multiplex_one_instance))
                    code_size = best_idxs[code_idx_idx]
                    size_spat_mux.append(spat_mux_one_instance[code_idx_idx])
                    code_counts[str(code_size)] += 1


            if idx_eps == 0:
                line, = ax.plot(switch_loss_x_axis, multiplexing, color=colors[ix_eta],
                        linestyle=linestyles[idx_eps],
                        linewidth=2.0, label="$\eta_{sp} = $" + str(eta_spin_gate_in)
                        )
                handles.append(line)
            else:
                print(size_spat_mux)
                ax.plot(switch_loss_x_axis, multiplexing, color=colors[ix_eta], linestyle=linestyles[idx_eps],
                    linewidth=2.0
                    )
    #         print(code_counts)
    # ax.plot([x * 1000 for x in dBS], [-1 for _ in dBS])

    line1, = ax.plot(1, -2, "-", color="black", linewidth=2.5, markersize=9.5,
                     label="$\lambda_{lp} = $ 0.01 %")
    line2, = ax.plot(1, -2, linestyle=linestyles[1], color="black", linewidth=2.5, markersize=9.5,
                     label="$\lambda_{lp} = $ 0.05 %")
    first_legend = ax.legend(handles=[line1, line2], bbox_to_anchor=(0, 0.35),
                             loc="lower left",
                             fontsize=16, title_fontsize=14)
    ax.add_artist(first_legend)
    second_legend = ax.legend(handles=handles, bbox_to_anchor=(0, 0.64), loc="lower left", fontsize=16,
                              title_fontsize=14)
    ax.add_artist(second_legend)


    ax.set_xscale("log")
    ax.set_yscale("log")
    plt.tick_params(axis='both', which='major', labelsize=30)
    # fig.savefig('FinalOverheadPlotFixedSwtichLossScanDiffP.pdf', dpi=600)
    plt.show()




def mux_plot_optimize_p_with_fixed_eta_s(t_bells, linestyles, p_th=1):
    colors = ["red", "blue"]
    fig, ax = plt.subplots(figsize=(9, 6))

    # The success probability of heralding a resource state given spin error which has to be multiplexed to average to one
    succ_lists = [[0.7866335111212663, 0.7261559912925004, 0.6703280704865696, 0.6187922808189714, 0.5712186370521632, 0.527302523689994],
                  [0.7866335111212663, 0.7261559912925004, 0.6703280704865696, 0.6187922808189714, 0.5712186370521632, 0.527302523689994]]

    # The tolerable transmission rates under 1% spin error and with eps = 0.05% and eps = 0.01% spin-photon error
    # etas_codes = [[0.9996984924623116, 0.9903517587939699, 0.987537688442211, 0.9874371859296482, 0.988140703517588,
    #                0.9892462311557789],
    #               [1, 0.9960804020100502, 0.9936683417085427, 0.9936683417085427, 0.9944723618090452, 0.9954773869346734]]
    eps_spin = 0.01
    rep_codes_sizes = [3, 4, 5, 6, 7, 8]
    sizes_code = [3 * 6, 4 * 6, 5 * 6, 6 * 6, 7 * 6, 8 * 6]
    epses = [0.01 / 100, 0.05 / 100]
    eta_spin_gates = [0.8, 0.9]
    handles = []
    dBS = np.linspace(0.0001, 0.005, 50)
    dBS = [0.0001, 0.0005]
    pair_gen_probs = np.linspace(0.03, 0.1, 10)
    pair_gen_probs = [0.03, 0.06, 0.09]
    # dBS = np.linspace(0.0001, 0.008, 400)
    for idx_eps, eps in enumerate(epses):
        succ_list = succ_lists[idx_eps]
        t_bell = t_bells[idx_eps]
        for ix_eta, eta_spin_gate_in in enumerate(eta_spin_gates):
            eta_spin_gate = eta_spin_gate_in * eta_spin_gate_in
            code_counts = {"0":0, "1":0, "2":0, "3":0, "4":0}
            multiplexing = []
            switch_loss_x_axis = []
            size_spat_mux = []
            for db in dBS:
                switch_loss = 1 - 10 ** (-db / 10)
                multiplex_one_instance = []
                spat_mux_one_instance = []
                best_idxs = []
                print("At DB: ", db)
                for ix_r, rep_code_size in enumerate(rep_codes_sizes):
                    number_of_photons = sizes_code[ix_r]
                    best_multiplex = 10 ** 12
                    best_spat_mux = 0
                    for pair_s in pair_gen_probs:
                        for k in range(2, 8):
                            temporal_size = k
                            for Ntemp in range(2 ** (k - 1), 2 ** k):
                                if Ntemp > 211:
                                    continue
                                else:
                                    temp_numb_input = Ntemp
                                    p_succ = ((1 - (1 - pair_s * eta_spin_gate) ** temp_numb_input) ** number_of_photons) # The total success probability of generating N photons with "temp_numb_input" temporal multiplexing attempts
                                    temp_max_delay = Ntemp - 1
                                    loss = delay_loss(temp_max_delay, t_bell)
                                    # Loss induced from delay lines and treating two photon clicks as photon loss
                                    # if idx_eps == 1 and ix_eta == 0:
                                    #     loss = delay_loss(temp_max_delay, t_bell) + 0.0022129873919999998
                                    # elif idx_eps == 1 and ix_eta == 1:
                                    #     loss = delay_loss(temp_max_delay, t_bell) + 0.001478206422
                                    # elif idx_eps == 0 and ix_eta == 1:
                                    #     loss = delay_loss(temp_max_delay, t_bell) + 0.0003695516055
                                    # else:
                                    #     loss = delay_loss(temp_max_delay, t_bell) + 0.0005532468479999999

                                    tot_prob = p_succ * (succ_list[ix_r]) # Multiply success probability with the error detection heralding probability

                                    # Find the binary size of a log tree for the spatial mux above
                                    # Want on average to get out 1 resource state, that is we need N*p = 1, giving N = 1 / p.
                                    # However, we use logarithmic trees, thus we find the smallest sized tree that satisfy this

                                    binary_sized_multi = -1
                                    numb_spin_inputs = 0
                                    for k_1 in range(1, 16):
                                        binary_size = 2 ** k_1
                                        mux_prob = binary_size * tot_prob
                                        if mux_prob >= 1:
                                            binary_sized_multi = k_1
                                            numb_spin_inputs = (2 ** k_1) / mux_prob
                                            break

                                    if binary_sized_multi > -1:
                                        multiplexing_loss = 1 - ((1 - switch_loss) ** (binary_sized_multi))  # Loss from the spatial multiplexing
                                        total_loss_with_multi = 1 - (1 - loss) * (1 - multiplexing_loss) * ((1 - switch_loss) ** (temporal_size + 1))  # total loss when doing spatial mux and temporal mux
                                        loss_from_detection = effective_loss_rate(pair_s, eta_spin_gate, 1 - total_loss_with_multi)
                                        total_loss_with_multi += loss_from_detection
                                        total_loss_with_photon_as_erase = (1 - p_succ * (1 - loss) * ((1 - switch_loss) ** (temporal_size + 1)))  # Total loss when when only doing temporal mux
                                        loss_from_detection = effective_loss_rate(pair_s, eta_spin_gate, 1 - total_loss_with_photon_as_erase)
                                        total_loss_with_photon_as_erase += loss_from_detection

                                        error_rate_spat = error_rate(pair_s, eta_spin_gate, 1 - total_loss_with_multi)
                                        error_rate_erase = error_rate(pair_s, eta_spin_gate, 1 - total_loss_with_photon_as_erase)

                                        spin_photon_error = (eps + error_rate_spat) * (2 / 3)
                                        threshold_loss_spat = get_threshold_given_eps_in_optimize_func(eps_spin, rep_code_size ,spin_err=spin_photon_error)

                                        spin_photon_error = (eps + error_rate_erase) * (2 / 3)
                                        threshold_loss_erase = get_threshold_given_eps_in_optimize_func(eps_spin, rep_code_size, spin_err=spin_photon_error)

                                        if total_loss_with_photon_as_erase < threshold_loss_erase:
                                            average_numb = 2 * number_of_photons  # Times two due to the spin encoding
                                            if average_numb < best_multiplex:
                                                best_multiplex = average_numb
                                        if total_loss_with_multi < threshold_loss_spat:
                                            average_numb = 2 * number_of_photons * numb_spin_inputs # (2 ** binary_sized_multi)  # Times two due to the spin encoding and the overhead from spatial mux
                                            if average_numb < best_multiplex:
                                                best_multiplex = average_numb
                                                best_spat_mux = numb_spin_inputs
                    if best_multiplex < 10 ** 6:
                        multiplex_one_instance.append(best_multiplex)
                        best_idxs.append(ix_r)
                        spat_mux_one_instance.append(best_spat_mux)

                if len(multiplex_one_instance) > 0:
                    multiplexing.append(min(multiplex_one_instance))
                    switch_loss_x_axis.append(db * 1000)
                    code_idx_idx = multiplex_one_instance.index(min(multiplex_one_instance))
                    code_size = best_idxs[code_idx_idx]
                    size_spat_mux.append(spat_mux_one_instance[code_idx_idx])
                    code_counts[str(code_size)] += 1


            if idx_eps == 0:
                line, = ax.plot(switch_loss_x_axis, multiplexing, color=colors[ix_eta],
                        linestyle=linestyles[idx_eps],
                        linewidth=2.0, label="$\eta_{sp} = $" + str(eta_spin_gate_in)
                        )
                handles.append(line)
            else:
                print(size_spat_mux)
                ax.plot(switch_loss_x_axis, multiplexing, color=colors[ix_eta], linestyle=linestyles[idx_eps],
                    linewidth=2.0
                    )
    #         print(code_counts)
    # ax.plot([x * 1000 for x in dBS], [-1 for _ in dBS])

    line1, = ax.plot(1, -2, "-", color="black", linewidth=2.5, markersize=9.5,
                     label="$\lambda_{lp} = $ 0.01 %")
    line2, = ax.plot(1, -2, linestyle=linestyles[1], color="black", linewidth=2.5, markersize=9.5,
                     label="$\lambda_{lp} = $ 0.05 %")
    first_legend = ax.legend(handles=[line1, line2], bbox_to_anchor=(0, 0.35),
                             loc="lower left",
                             fontsize=16, title_fontsize=14)
    ax.add_artist(first_legend)
    second_legend = ax.legend(handles=handles, bbox_to_anchor=(0, 0.64), loc="lower left", fontsize=16,
                              title_fontsize=14)
    ax.add_artist(second_legend)


    ax.set_xscale("log")
    ax.set_yscale("log")
    plt.tick_params(axis='both', which='major', labelsize=30)
    # fig.savefig('FinalOverheadPlotFixedSwtichLossScanDiffP.pdf', dpi=600)
    plt.show()

if __name__ == '__main__':
    gridresolution = 100
    t_bells = 10 * 10 ** (-9)
    mux_flag = True
    linestyles = ["solid", "dashed", "dashdot"]

    ###################################################################
    ####  PLOT MULTIPLEXING FOR PARITY CODES WITH SPIN ERROR RATES ####
    ###################################################################


    # pair_source_prob_vs_multiplexing_parity_codes_spin_error_fixed_scan_switch_loss(t_bells, linestyles)
    # pair_source_prob_vs_multiplexing_parity_codes_spin_error_fixed_scan_switch_loss_test(t_bells, linestyles)

    t_bells = [9.5 * (10 ** (-9)), 9.5 * (10 ** (-9))]
    mux_plot_optimize_p(t_bells, linestyles)