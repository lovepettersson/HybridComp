import numpy as np
import matplotlib.pyplot as plt
from ErrorsFromHigherOrderPhotonComponents import error_rate, effective_loss_rate
from RepetionCodes import get_threshold_given_eps, get_threshold_given_eps_one_code


def delay_loss(N, t_bell):
    # Delay loss using telecom photons and having a attenuation length of L_att = 20 km
    time = N * t_bell / 2   # Time is number of tries times the time per gate divided by two. This is the average time following a uniform distribution.
    speed = 2 * (10 ** 8)
    L = time * speed
    transmission = np.exp(-L / 20000)
    loss = 1 - transmission
    return loss



# def pair_source_prob_vs_multiplexing_parity_codes_spin_error_fixed(t_bell, linestyles, switch_loss = 0.6884 / 100, gridresolution=100,mux_flag=True):
#     colors = ["red", "blue"]
#     fig, ax = plt.subplots(figsize=(9, 6))
#     pair_succ = np.linspace(0.005, 0.1, gridresolution)

#     etas_codes = [[0.991859296482412, 0.9883417085427135, 0.9872361809045226, 0.9872361809045226, 0.9876381909547739],
#                   [0.9926633165829145, 0.989145728643216, 0.9880402010050251, 0.9880402010050251, 0.9885427135678392]]  # Transmission thresholds for the spin error and two different photon link error rates
#     succ_lists = [[0.7866335111212663, 0.7261559912925004, 0.6703280704865696, 0.6187922808189714, 0.5712186370521632, 0.527302523689994],
#                   [0.7866335111212663, 0.7261559912925004, 0.6703280704865696, 0.6187922808189714, 0.5712186370521632, 0.527302523689994]]  # The success probability which has to be multiplexed to average to one

#     sizes_code = [3 * 6, 4 * 6, 5 * 6, 6 * 6, 7 * 6, 8 * 6]
#     epses = ["0.001 %", "0.01 %"]

#     eta_spin_gates = [0.8, 0.9]
#     handles = []
#     numb_tries = [x for x in range(1, 700)]
#     for idx_eps, eps in enumerate(epses):
#         etas_code = etas_codes[idx_eps]
#         succ_list = succ_lists[idx_eps]
#         for ix_eta, eta_spin_gate_in in enumerate(eta_spin_gates):
#             eta_spin_gate = eta_spin_gate_in * eta_spin_gate_in
#             code_counts = {"0":0, "1":0, "2":0, "3":0}
#             multiplexing = []
#             spat_mux_flags = []
#             pair_succ_plot = []
#             for pair_s in pair_succ:
#                 multiplex_one_instance = []
#                 spat_mux = []
#                 best_idxs = []
#                 for ix_r in range(len(etas_code)):
#                     threshold_loss = 1 - etas_code[ix_r]
#                     number_of_photons = sizes_code[ix_r]
#                     best_multiplex = 10 ** 12
#                     flag = True
#                     for N in numb_tries:
#                         # binary_sized_time_mux = N
#                         for k in range(1, 16):
#                             binary_size = 2 ** k
#                             if N - binary_size < 1:
#                                 binary_sized_time_mux = k
#                                 break
#                         temp_numb_input = 2 ** binary_sized_time_mux
#                         if mux_flag:
#                             p_succ = (1 - ((1 - pair_s * eta_spin_gate) ** (temp_numb_input))) * ((1 - switch_loss) ** (binary_sized_time_mux + 1))  # The total loss from time mux including non unity succ prob and switch loss
#                         else:
#                             p_succ = (1 - ((1 - pair_s * eta_spin_gate) ** (temp_numb_input)))
#                         # loss = delay_loss(N, t_bell)
#                         temp_max_delay = 2 ** binary_sized_time_mux - 1
#                         loss = delay_loss(temp_max_delay, t_bell)  # max number of delay
#                         multiplxing_size = 1 / (p_succ ** number_of_photons)  # The required multiplexing size to get one resource state out on average
#                         binary_sized_multi = multiplxing_size


                        # Find the binary size of a log tree for the spatial mux above
#                         for k in range(1, 10):
#                             binary_size = 2 ** k
#                             if multiplxing_size - binary_size < 1:
#                                 binary_sized_multi = k
#                                 break

#                         # Multiplexing overhead from throwing away the the error detection
#                         if succ_list[ix_r] != 1:
#                             spat_from_detect_size = 1 / succ_list[ix_r]
#                             for k in range(1, 10):
#                                 binary_size_detect = 2 ** k
#                                 if spat_from_detect_size - binary_size_detect < 1:
#                                     binary_sized_multi_detect = k
#                                     break
#                         else:
#                             binary_sized_multi_detect = 0
#                             spat_from_detect_size = 1 / succ_list[ix_r]


#                         multiplexing_loss = 1 - ((1 - switch_loss) ** binary_sized_multi) * ((1 - switch_loss) ** binary_sized_multi_detect)  # Loss from the spatial multiplexing
#                         total_loss_with_multi = 1 - (1 - loss) * (1 - multiplexing_loss)  # total loss when doing spatial mux
#                         total_loss_with_photon_as_erase = (1 - p_succ * (1 - loss) * ((1 - switch_loss) ** binary_sized_multi_detect))  # Total loss when when only doing temporal mux

#                         if total_loss_with_photon_as_erase < threshold_loss:
#                             average_numb = 2 * number_of_photons * spat_from_detect_size  # Times two due to the spin encoding
#                             if average_numb < best_multiplex:
#                                 best_multiplex = average_numb
#                                 flag = True
#                         if total_loss_with_multi < threshold_loss:
#                             average_numb = 2 * number_of_photons * multiplxing_size * spat_from_detect_size  # Times two due to the spin encoding
#                             if average_numb < best_multiplex:
#                                 best_multiplex = average_numb
#                                 flag = False
#                     if best_multiplex < 10 ** 12:
#                         multiplex_one_instance.append(best_multiplex)
#                         spat_mux.append(flag)
#                         best_idxs.append(ix_r)
#                 if len(multiplex_one_instance) > 0:
#                     multiplexing.append(min(multiplex_one_instance))
#                     pair_succ_plot.append(pair_s)
#                     code_idx_idx = multiplex_one_instance.index(min(multiplex_one_instance))
#                     code_size = best_idxs[code_idx_idx]
#                     code_counts[str(code_size)] += 1
#                     flag = spat_mux[multiplex_one_instance.index(min(multiplex_one_instance))]
#                    if flag:
#                         spat_mux_flags.append(-1)
#                     else:
#                         spat_mux_flags.append(1)

#             if idx_eps == 0:
#                line, = ax.plot(pair_succ_plot, multiplexing, color=colors[ix_eta],
#                         linestyle=linestyles[idx_eps],
#                         linewidth=2.0, label="$\eta_{sp} = $" + str(eta_spin_gate_in)
#                         )
#                 handles.append(line)
#             else:
#                 ax.plot(pair_succ_plot, multiplexing, color=colors[ix_eta], linestyle=linestyles[idx_eps],
#                     linewidth=2.0
#                     )
#             print(code_counts)
#     ax.plot(pair_succ, [-1 for x in pair_succ])

#     line1, = ax.plot(0, -2, "-", color="black", linewidth=2.5, markersize=9.5,
#                      label="$\lambda_{lp} = $ 0.001 %")
#     line2, = ax.plot(0, -2, linestyle=linestyles[1], color="black", linewidth=2.5, markersize=9.5,
#                      label="$\lambda_{lp} = $ 0.01 %")
#     # first_legend = ax.legend(handles=[line1, line2, line3], bbox_to_anchor=(0.68, 0.35),
#     #                           loc="lower left",
#     #                          fontsize=16, title_fontsize=16)
#     first_legend = ax.legend(handles=[line1, line2], bbox_to_anchor=(0, 0.35),
#                              loc="lower left",
#                              fontsize=16, title_fontsize=14)
#     ax.add_artist(first_legend)
#     # second_legend = ax.legend(handles=handles, bbox_to_anchor=(0.71, 0.64), loc="lower left", fontsize=16,
#     #           title_fontsize=16)
#     second_legend = ax.legend(handles=handles, bbox_to_anchor=(0, 0.64), loc="lower left", fontsize=16,
#                               title_fontsize=14)
#     ax.add_artist(second_legend)


#     ax.set_xscale("log")
#     # ax.set_ylim((70, 2 * 10**2))
#     ax.set_yscale("log")
#     # ax.set_ylabel("# of spins")
#     # ax.set_xlabel("Photon pair prob.")
#     plt.tick_params(axis='both', which='major', labelsize=30)
#     # ax.legend(fontsize=18)
#     # fig.savefig('FinalOverheadPlotFixed.pdf', dpi=600)
#     plt.show()



# def pair_source_prob_vs_multiplexing_parity_codes_spin_error_fixed(t_bell, linestyles, switch_loss = 0.6884 / 100, gridresolution=100,mux_flag=True, p_th=1):
#     colors = ["red", "blue"]
#     fig, ax = plt.subplots(figsize=(9, 6))
#     pair_succ = np.linspace(0.005, 0.1, gridresolution)
#
#     succ_lists = [[0.7866335111212663, 0.7261559912925004, 0.6703280704865696, 0.6187922808189714, 0.5712186370521632, 0.527302523689994],
#                  [0.7866335111212663, 0.7261559912925004, 0.6703280704865696, 0.6187922808189714, 0.5712186370521632, 0.527302523689994]]  # The success probability which has to be multiplexed to average to one
#
#     etas_codes = [[0.9984924623115577, 0.9887437185929648, 0.9857286432160803, 0.9854271356783919, 0.9860301507537689, 0.9870351758793969],
#                   [0.9989949748743718, 0.9894472361809045, 0.9866331658291457, 0.9863316582914573, 0.9870351758793969, 0.988140703517588]]
#
#     sizes_code = [3 * 6, 4 * 6, 5 * 6, 6 * 6, 7 * 6, 8 * 6]
#     epses = ["0.001 %", "0.005 %"]
#
#     eta_spin_gates = [0.8, 0.9]
#     handles = []
#     for idx_eps, eps in enumerate(epses):
#         etas_code = etas_codes[idx_eps]
#         succ_list = succ_lists[idx_eps]
#         for ix_eta, eta_spin_gate_in in enumerate(eta_spin_gates):
#             eta_spin_gate = eta_spin_gate_in * eta_spin_gate_in
#             code_counts = {"0":0, "1":0, "2":0, "3":0, "4":0}
#             multiplexing = []
#             spat_mux_flags = []
#             pair_succ_plot = []
#             for pair_s in pair_succ:
#                 multiplex_one_instance = []
#                 probs_list = []
#                 spat_mux = []
#                 best_idxs = []
#                 for ix_r in range(len(etas_code)):
#                     threshold_loss = 1 - etas_code[ix_r]
#                     number_of_photons = sizes_code[ix_r]
#                     best_multiplex = 10 ** 12
#                     best_prob = 0
#                     flag = True
#                     for k in range(1, 16):
#                         temp_numb_input = 2 ** k
#                         if mux_flag:
#                             p_succ = (1 - ((1 - pair_s * eta_spin_gate) ** (temp_numb_input))) * ((1 - switch_loss) ** (k + 1))  # The total loss from time mux including non unity succ prob and switch loss
#                         else:
#                             p_succ = (1 - ((1 - pair_s * eta_spin_gate) ** (temp_numb_input)))
#                         # loss = delay_loss(N, t_bell)
#                         temp_max_delay = 2 ** k - 1
#                         loss = delay_loss(temp_max_delay, t_bell)  # max number of delay
#                         multiplxing_size = 1 / (p_succ ** number_of_photons)  # The required multiplexing size to get one resource state out on average
#                         binary_sized_multi = multiplxing_size
#
#                         tot_prob = (succ_list[ix_r] * (p_succ ** number_of_photons))

                        # Find the binary size of a log tree for the spatial mux above
                        # Want on average to get out 1 resource state, that is we need N*p = 1, giving N = 1 / p.
                        # However, we use logarithmic trees, thus we find the smallest sized tree that satisfy this
#                         for k_1 in range(1, 12):
#                             binary_size = 2 ** k_1
#                             mux_prob = binary_size * tot_prob
#                             if mux_prob > p_th:
#                                 binary_sized_multi = k_1
#                                 break

                        # Multiplexing overhead from throwing away the the error detection
#                         if succ_list[ix_r] != 1:
#                             detect_prob = succ_list[ix_r]
#                             for k_2 in range(1, 10):
#                                 binary_size_detect = 2 ** k_2
#                                 detect_prob_mux = binary_size_detect * (1 - (1 - detect_prob))  # ** (binary_size_detect)
#                                 if detect_prob_mux > p_th:
#                                     binary_sized_multi_detect = k_2
#                                     break
#                             else:
#                                 binary_sized_multi_detect = 0


#                         multiplexing_loss = 1 - ((1 - switch_loss) ** (binary_sized_multi + binary_sized_multi_detect))  # Loss from the spatial multiplexing
#                         total_loss_with_multi = 1 - (1 - loss) * (1 - multiplexing_loss)  # total loss when doing spatial mux
#                         total_loss_with_photon_as_erase = (1 - p_succ * (1 - loss) * ((1 - switch_loss) ** binary_sized_multi_detect))  # Total loss when when only doing temporal mux
#                         prob_erase = (1 - (1 - p_succ ** number_of_photons) ** (2 ** binary_sized_multi)) * (1 - total_loss_with_photon_as_erase)
#                         prob_spat = (1 - ((1 - p_succ ** number_of_photons) ** (2 ** binary_sized_multi))) * (1 - total_loss_with_multi)
#                         if total_loss_with_photon_as_erase < threshold_loss: #  and prob_erase > p_th:
#                             average_numb = 2 * number_of_photons * (2 ** binary_sized_multi_detect)  # Times two due to the spin encoding
#                             if average_numb < best_multiplex:
#                                 best_multiplex = average_numb
#                                 best_prob = prob_erase
#                                 flag = True
#                         if total_loss_with_multi < threshold_loss: # and prob_spat > p_th:
#                             average_numb = 2 * number_of_photons * (2 ** binary_sized_multi) * (2 ** binary_sized_multi_detect)  # Times two due to the spin encoding
#                             if average_numb < best_multiplex:
#                                 best_multiplex = average_numb
#                                 best_prob = prob_spat
#                                 flag = False
#                     if best_multiplex < 10 ** 12:
#                         multiplex_one_instance.append(best_multiplex)
#                         spat_mux.append(flag)
#                         best_idxs.append(ix_r)
#
#                 if len(multiplex_one_instance) > 0:
#                     multiplexing.append(min(multiplex_one_instance))
#                     pair_succ_plot.append(pair_s)
#                     code_idx_idx = multiplex_one_instance.index(min(multiplex_one_instance))
#                     code_size = best_idxs[code_idx_idx]
#                     code_counts[str(code_size)] += 1
#                     flag = spat_mux[multiplex_one_instance.index(min(multiplex_one_instance))]
#                     if flag:
#                         spat_mux_flags.append(-1)
#                     else:
#                         spat_mux_flags.append(1)
#
#             if idx_eps == 0:
#                 line, = ax.plot(pair_succ_plot, multiplexing, color=colors[ix_eta],
#                         linestyle=linestyles[idx_eps],
#                         linewidth=2.0, label="$\eta_{sp} = $" + str(eta_spin_gate_in)
#                         )
#                 handles.append(line)
#             else:
#                 ax.plot(pair_succ_plot, multiplexing, color=colors[ix_eta], linestyle=linestyles[idx_eps],
#                     linewidth=2.0
#                     )
#             print(code_counts)
#     ax.plot(pair_succ, [-1 for x in pair_succ])
#
#     line1, = ax.plot(0, -2, "-", color="black", linewidth=2.5, markersize=9.5,
#                      label="$\lambda_{lp} = $ 0.001 %")
#     line2, = ax.plot(0, -2, linestyle=linestyles[1], color="black", linewidth=2.5, markersize=9.5,
#                      label="$\lambda_{lp} = $ 0.005 %")
#     # first_legend = ax.legend(handles=[line1, line2, line3], bbox_to_anchor=(0.68, 0.35),
#     #                           loc="lower left",
#     #                          fontsize=16, title_fontsize=16)
#     first_legend = ax.legend(handles=[line1, line2], bbox_to_anchor=(0, 0.35),
#                              loc="lower left",
#                              fontsize=16, title_fontsize=14)
#     ax.add_artist(first_legend)
#     # second_legend = ax.legend(handles=handles, bbox_to_anchor=(0.71, 0.64), loc="lower left", fontsize=16,
#     #           title_fontsize=16)
#     second_legend = ax.legend(handles=handles, bbox_to_anchor=(0, 0.64), loc="lower left", fontsize=16,
#                               title_fontsize=14)
#     ax.add_artist(second_legend)
#
#    ax.set_xscale("log")
#    # ax.set_ylim((70, 2 * 10**2))
#     ax.set_yscale("log")
#     # ax.set_ylabel("# of spins")
#     # ax.set_xlabel("Photon pair prob.")
#     plt.tick_params(axis='both', which='major', labelsize=30)
#     # ax.legend(fontsize=18)
#     # fig.savefig('FinalOverheadPlotFixed.pdf', dpi=600)
#     plt.show()









def pair_source_prob_vs_multiplexing_parity_codes_spin_error_fixed_scan_switch_loss(t_bell, linestyles,
                                                                                    p_th=1):
    colors = ["red", "blue"]
    fig, ax = plt.subplots(figsize=(9, 6))

    # The success probability of heralding a resource state given spin error which has to be multiplexed to average to one
    succ_lists = [[0.7866335111212663, 0.7261559912925004, 0.6703280704865696, 0.6187922808189714, 0.5712186370521632, 0.527302523689994],
                  [0.7866335111212663, 0.7261559912925004, 0.6703280704865696, 0.6187922808189714, 0.5712186370521632, 0.527302523689994]]

    # The tolerable transmission rates under 1% spin error and with eps = 0.05% and eps = 0.01% spin-photon error
    etas_codes = [[0.9996984924623116, 0.9903517587939699, 0.987537688442211, 0.9874371859296482, 0.988140703517588,
                   0.9892462311557789],
                  [1, 0.9960804020100502, 0.9936683417085427, 0.9936683417085427, 0.9944723618090452, 0.9954773869346734]]


    sizes_code = [3 * 6, 4 * 6, 5 * 6, 6 * 6, 7 * 6, 8 * 6]
    epses = ["0.01 %", "0.05 %"]
    eta_spin_gates = [0.8, 0.9]
    handles = []
    dBS = np.linspace(0.0001, 0.005, 1000)
    for idx_eps, eps in enumerate(epses):
        if idx_eps == 0:
            pair_s = 0.035
        else:
            pair_s = 0.07
        etas_code = etas_codes[idx_eps]
        succ_list = succ_lists[idx_eps]
        for ix_eta, eta_spin_gate_in in enumerate(eta_spin_gates):
            eta_spin_gate = eta_spin_gate_in * eta_spin_gate_in
            code_counts = {"0":0, "1":0, "2":0, "3":0, "4":0}
            multiplexing = []
            switch_loss_x_axis = []
            for db in dBS:
                switch_loss = 1 - 10 ** (-db / 10)
                multiplex_one_instance = []
                best_idxs = []
                for ix_r in range(len(etas_code)):
                    threshold_loss = 1 - etas_code[ix_r]
                    number_of_photons = sizes_code[ix_r]
                    best_multiplex = 10 ** 12
                    for k in range(1, 8):
                    # for k in range(1, 16):
                        temp_numb_input = 2 ** k
                        temporal_size = k
                        p_succ = ((1 - (1 - pair_s * eta_spin_gate) ** temp_numb_input) ** number_of_photons) # The total success probability of generating N photons with "temp_numb_input" temporal multiplexing attempts
                        temp_max_delay = 2 ** k - 1
                        # loss = delay_loss(temp_max_delay, t_bell)

                        # Loss induced from delay lines and treating two photon clicks as photon loss
                        if idx_eps == 1 and ix_eta == 0:
                            loss = delay_loss(temp_max_delay, t_bell) + 0.0022129873919999998
                        elif idx_eps == 1 and ix_eta == 1:
                            loss = delay_loss(temp_max_delay, t_bell) + 0.001478206422
                        elif idx_eps == 0 and ix_eta == 1:
                            loss = delay_loss(temp_max_delay, t_bell) + 0.0003695516055
                        else:
                            loss = delay_loss(temp_max_delay, t_bell) + 0.0005532468479999999

                        tot_prob = p_succ * (succ_list[ix_r]) # Multiply success probability with the error detection heralding probability

                        # Find the binary size of a log tree for the spatial mux above
                        # Want on average to get out 1 resource state, that is we need N*p = 1, giving N = 1 / p.
                        # However, we use logarithmic trees, thus we find the smallest sized tree that satisfy this
                        binary_sized_multi = -1
                        for k_1 in range(1, 16):
                            binary_size = 2 ** k_1
                            mux_prob = binary_size * tot_prob
                            if mux_prob >= 1:
                                binary_sized_multi = k_1
                                break

                        if binary_sized_multi > -1:
                            multiplexing_loss = 1 - ((1 - switch_loss) ** (binary_sized_multi))  # Loss from the spatial multiplexing
                            total_loss_with_multi = 1 - (1 - loss) * (1 - multiplexing_loss) * ((1 - switch_loss) ** (temporal_size + 1))  # total loss when doing spatial mux and temporal mux
                            # total_loss_with_multi = 1 - (1 - multiplexing_loss) * ((1 - switch_loss) ** (temporal_size))  # Sanity check where we only have spat mux and find max switch eff.
                            total_loss_with_photon_as_erase = (1 - p_succ * (1 - loss) * ((1 - switch_loss) ** (temporal_size + 1)))  # Total loss when when only doing temporal mux
                            if total_loss_with_photon_as_erase < threshold_loss:
                                average_numb = 2 * number_of_photons  # Times two due to the spin encoding
                                if average_numb < best_multiplex:
                                    best_multiplex = average_numb
                            if total_loss_with_multi < threshold_loss:
                                average_numb = 2 * number_of_photons * (2 ** binary_sized_multi)  # Times two due to the spin encoding and the overhead from spatial mux
                                if average_numb < best_multiplex:
                                    best_multiplex = average_numb
                    if best_multiplex < 10 ** 12:
                        multiplex_one_instance.append(best_multiplex)
                        best_idxs.append(ix_r)

                if len(multiplex_one_instance) > 0:
                    multiplexing.append(min(multiplex_one_instance))
                    switch_loss_x_axis.append(db * 1000)
                    code_idx_idx = multiplex_one_instance.index(min(multiplex_one_instance))
                    code_size = best_idxs[code_idx_idx]
                    code_counts[str(code_size)] += 1


            if idx_eps == 0:
                line, = ax.plot(switch_loss_x_axis, multiplexing, color=colors[ix_eta],
                        linestyle=linestyles[idx_eps],
                        linewidth=2.0, label="$\eta_{sp} = $" + str(eta_spin_gate_in)
                        )
                handles.append(line)
            else:
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



def pair_source_prob_vs_multiplexing_parity_codes_spin_error_fixed_scan_switch_loss_new(t_bells, linestyles,
                                                                                    p_th=1):
    colors = ["red", "blue"]
    fig, ax = plt.subplots(figsize=(9, 6))

    # The success probability of heralding a resource state given spin error which has to be multiplexed to average to one
    succ_lists = [[0.7866335111212663, 0.7261559912925004, 0.6703280704865696, 0.6187922808189714, 0.5712186370521632, 0.527302523689994],
                  [0.7866335111212663, 0.7261559912925004, 0.6703280704865696, 0.6187922808189714, 0.5712186370521632, 0.527302523689994]]

    # The tolerable transmission rates under 1% spin error and with eps = 0.05% and eps = 0.01% spin-photon error
    etas_codes = [[0.9996984924623116, 0.9903517587939699, 0.987537688442211, 0.9874371859296482, 0.988140703517588,
                   0.9892462311557789],
                  [1, 0.9960804020100502, 0.9936683417085427, 0.9936683417085427, 0.9944723618090452, 0.9954773869346734]]


    sizes_code = [3 * 6, 4 * 6, 5 * 6, 6 * 6, 7 * 6, 8 * 6]
    epses = ["0.01 %", "0.05 %"]
    eta_spin_gates = [0.8, 0.9]
    handles = []
    dBS = np.linspace(0.0001, 0.005, 400)
    # dBS = np.linspace(0.0001, 0.008, 400)
    for idx_eps, eps in enumerate(epses):
        if idx_eps == 0:
            pair_s = 0.035
        else:
            pair_s = 0.07
        etas_code = etas_codes[idx_eps]
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
                for ix_r in range(len(etas_code)):
                    threshold_loss = 1 - etas_code[ix_r]
                    number_of_photons = sizes_code[ix_r]
                    best_multiplex = 10 ** 12
                    best_spat_mux = 0
                    for k in range(2, 8):
                        temporal_size = k
                        for Ntemp in range(2 ** (k - 1), 2 ** k):
                            if Ntemp > 211:
                                continue
                            else:
                                temp_numb_input = Ntemp
                                p_succ = ((1 - (1 - pair_s * eta_spin_gate) ** temp_numb_input) ** number_of_photons) # The total success probability of generating N photons with "temp_numb_input" temporal multiplexing attempts
                                temp_max_delay = Ntemp - 1
                                # loss = delay_loss(temp_max_delay, t_bell)

                                # Loss induced from delay lines and treating two photon clicks as photon loss
                                if idx_eps == 1 and ix_eta == 0:
                                    loss = delay_loss(temp_max_delay, t_bell) + 0.0022129873919999998
                                elif idx_eps == 1 and ix_eta == 1:
                                    loss = delay_loss(temp_max_delay, t_bell) + 0.001478206422
                                elif idx_eps == 0 and ix_eta == 1:
                                    loss = delay_loss(temp_max_delay, t_bell) + 0.0003695516055
                                else:
                                    loss = delay_loss(temp_max_delay, t_bell) + 0.0005532468479999999

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
                                    # total_loss_with_multi = 1 - (1 - multiplexing_loss) * ((1 - switch_loss) ** (temporal_size))  # Sanity check where we only have spat mux and find max switch eff.
                                    total_loss_with_photon_as_erase = (1 - p_succ * (1 - loss) * ((1 - switch_loss) ** (temporal_size + 1)))  # Total loss when when only doing temporal mux
                                    if total_loss_with_photon_as_erase < threshold_loss:
                                        average_numb = 2 * number_of_photons  # Times two due to the spin encoding
                                        if average_numb < best_multiplex:
                                            best_multiplex = average_numb
                                    if total_loss_with_multi < threshold_loss:
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





def pair_source_prob_vs_multiplexing_parity_codes_spin_error_fixed_scan_switch_loss_test(t_bell, linestyles,
                                                                                    p_th=1):
    colors = ["red", "blue"]
    fig, ax = plt.subplots(figsize=(9, 6))

    # The success probability of heralding a resource state given spin error which has to be multiplexed to average to one
    succ_lists = [[0.7866335111212663, 0.7261559912925004, 0.6703280704865696, 0.6187922808189714, 0.5712186370521632, 0.527302523689994],
                  [0.7866335111212663, 0.7261559912925004, 0.6703280704865696, 0.6187922808189714, 0.5712186370521632, 0.527302523689994]]

    # The tolerable transmission rates under 1% spin error and with eps = 0.05% and eps = 0.01% spin-photon error
    etas_codes = [[0.9996984924623116, 0.9903517587939699, 0.987537688442211, 0.9874371859296482, 0.988140703517588,
                   0.9892462311557789],
                  [1, 0.9960804020100502, 0.9936683417085427, 0.9936683417085427, 0.9944723618090452, 0.9954773869346734]]


    sizes_code = [3 * 6, 4 * 6, 5 * 6, 6 * 6, 7 * 6, 8 * 6]
    epses = ["0.01 %", "0.05 %"]
    eta_spin_gates = [0.8, 0.9]
    handles = []
    dBS = np.linspace(0.0001, 0.005, 1000)
    for idx_eps, eps in enumerate(epses):
        if idx_eps == 0:
            pair_s = 0.035
        else:
            pair_s = 0.07
        etas_code = etas_codes[idx_eps]
        succ_list = succ_lists[idx_eps]
        for ix_eta, eta_spin_gate_in in enumerate(eta_spin_gates):
            eta_spin_gate = eta_spin_gate_in * eta_spin_gate_in
            code_counts = {"0":0, "1":0, "2":0, "3":0, "4":0}
            multiplexing = []
            switch_loss_x_axis = []
            spat_muxes = []
            temp_muxes = []
            code_sizes_used = []
            p_gen_used = []
            for db in dBS:
                switch_loss = 1 - 10 ** (-db / 10)
                multiplex_one_instance = []
                best_idxs = []
                spat_ks = []
                temp_ks = []
                code_size_succ = []
                photon_gen_probs = []
                for ix_r in range(len(etas_code)):
                    threshold_loss = 1 - etas_code[ix_r]
                    number_of_photons = sizes_code[ix_r]
                    best_multiplex = 10 ** 12
                    best_spat_k = 0
                    best_temp_k = 0
                    best_code_size = 0
                    best_photon_gen = 0
                    for k in range(1, 8):
                    # for k in range(1, 16):
                        temp_numb_input = 2 ** k
                        temporal_size = k
                        p_succ = ((1 - (1 - pair_s * eta_spin_gate) ** temp_numb_input) ** number_of_photons) # The total success probability of generating N photons with "temp_numb_input" temporal multiplexing attempts
                        temp_max_delay = 2 ** k - 1
                        # loss = delay_loss(temp_max_delay, t_bell)

                        # Loss induced from delay lines and treating two photon clicks as photon loss
                        if idx_eps == 1 and ix_eta == 0:
                            loss = delay_loss(temp_max_delay, t_bell) + 0.0022129873919999998
                        elif idx_eps == 1 and ix_eta == 1:
                            loss = delay_loss(temp_max_delay, t_bell) + 0.001478206422
                        elif idx_eps == 0 and ix_eta == 1:
                            loss = delay_loss(temp_max_delay, t_bell) + 0.0003695516055
                        else:
                            loss = delay_loss(temp_max_delay, t_bell) + 0.0005532468479999999

                        tot_prob = p_succ * (succ_list[ix_r]) # Multiply success probability with the error detection heralding probability

                        # Find the binary size of a log tree for the spatial mux above
                        # Want on average to get out 1 resource state, that is we need N*p = 1, giving N = 1 / p.
                        # However, we use logarithmic trees, thus we find the smallest sized tree that satisfy this
                        binary_sized_multi = -1
                        for k_1 in range(1, 16):
                            binary_size = 2 ** k_1
                            mux_prob = binary_size * tot_prob
                            if mux_prob >= 1:
                                binary_sized_multi = k_1
                                break

                        if binary_sized_multi > -1:
                            multiplexing_loss = 1 - ((1 - switch_loss) ** (binary_sized_multi))  # Loss from the spatial multiplexing
                            total_loss_with_multi = 1 - (1 - loss) * (1 - multiplexing_loss) * ((1 - switch_loss) ** (temporal_size + 1))  # total loss when doing spatial mux and temporal mux
                            # total_loss_with_multi = 1 - (1 - multiplexing_loss) * ((1 - switch_loss) ** (temporal_size))  # Sanity check where we only have spat mux and find max switch eff.
                            total_loss_with_photon_as_erase = (1 - p_succ * (1 - loss) * ((1 - switch_loss) ** (temporal_size + 1)))  # Total loss when when only doing temporal mux
                            if total_loss_with_photon_as_erase < threshold_loss:
                                average_numb = 2 * number_of_photons  # Times two due to the spin encoding
                                if average_numb < best_multiplex:
                                    best_multiplex = average_numb
                            if total_loss_with_multi < threshold_loss:
                                average_numb = 2 * number_of_photons * (2 ** binary_sized_multi)  # Times two due to the spin encoding and the overhead from spatial mux
                                if average_numb < best_multiplex:
                                    best_multiplex = average_numb
                                best_spat_k = binary_sized_multi
                                best_temp_k = k
                                best_code_size = sizes_code[ix_r] / 6
                                best_photon_gen = p_succ
                    if best_multiplex < 10 ** 12:
                        multiplex_one_instance.append(best_multiplex)
                        best_idxs.append(ix_r)
                        spat_ks.append(best_spat_k)
                        temp_ks.append(best_temp_k)
                        code_size_succ.append(best_code_size)
                        photon_gen_probs.append(best_photon_gen)

                if len(multiplex_one_instance) > 0:
                    multiplexing.append(min(multiplex_one_instance))
                    switch_loss_x_axis.append(db * 1000)
                    code_idx_idx = multiplex_one_instance.index(min(multiplex_one_instance))
                    spat_muxes.append(spat_ks[code_idx_idx])
                    temp_muxes.append(temp_ks[code_idx_idx])
                    code_size = best_idxs[code_idx_idx]
                    code_counts[str(code_size)] += 1
                    code_sizes_used.append(code_size)
                    p_gen_used.append(photon_gen_probs[code_idx_idx])


            if idx_eps == 0:
                line, = ax.plot(switch_loss_x_axis, multiplexing, color=colors[ix_eta],
                        linestyle=linestyles[idx_eps],
                        linewidth=2.0, label="$\eta_{sp} = $" + str(eta_spin_gate_in)
                        )
                handles.append(line)
                print("Low error", str(eta_spin_gate_in))
                print(spat_muxes)
                print(temp_muxes)
                print(code_sizes_used)
                print()
            else:
                ax.plot(switch_loss_x_axis, multiplexing, color=colors[ix_eta], linestyle=linestyles[idx_eps],
                    linewidth=2.0
                    )
                print(spat_muxes)
                print(temp_muxes)
                print(code_sizes_used)
                print(p_gen_used)
                print()
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



def mux_plot_optimize_p(t_bells, linestyles, p_th=1):
    colors = ["red", "blue"]
    fig, ax = plt.subplots(figsize=(9, 6))

    # The success probability of heralding a resource state given spin error which has to be multiplexed to average to one
    succ_lists = [[0.7866335111212663, 0.7261559912925004, 0.6703280704865696, 0.6187922808189714, 0.5712186370521632, 0.527302523689994],
                  [0.7866335111212663, 0.7261559912925004, 0.6703280704865696, 0.6187922808189714, 0.5712186370521632, 0.527302523689994]]

    # The tolerable transmission rates under 1% spin error and with eps = 0.05% and eps = 0.01% spin-photon error
    etas_codes = [[0.9996984924623116, 0.9903517587939699, 0.987537688442211, 0.9874371859296482, 0.988140703517588,
                   0.9892462311557789],
                  [1, 0.9960804020100502, 0.9936683417085427, 0.9936683417085427, 0.9944723618090452, 0.9954773869346734]]


    sizes_code = [3 * 6, 4 * 6, 5 * 6, 6 * 6, 7 * 6, 8 * 6]
    epses = ["0.01 %", "0.05 %"]
    eta_spin_gates = [0.8, 0.9]
    handles = []
    dBS = np.linspace(0.0001, 0.005, 400)
    # dBS = np.linspace(0.0001, 0.008, 400)
    for idx_eps, eps in enumerate(epses):
        if idx_eps == 0:
            pair_s = 0.035
        else:
            pair_s = 0.07
        etas_code = etas_codes[idx_eps]
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
                for ix_r in range(len(etas_code)):
                    threshold_loss = 1 - etas_code[ix_r]
                    number_of_photons = sizes_code[ix_r]
                    best_multiplex = 10 ** 12
                    best_spat_mux = 0
                    for k in range(2, 8):
                        temporal_size = k
                        for Ntemp in range(2 ** (k - 1), 2 ** k):
                            if Ntemp > 211:
                                continue
                            else:
                                temp_numb_input = Ntemp
                                p_succ = ((1 - (1 - pair_s * eta_spin_gate) ** temp_numb_input) ** number_of_photons) # The total success probability of generating N photons with "temp_numb_input" temporal multiplexing attempts
                                temp_max_delay = Ntemp - 1
                                # loss = delay_loss(temp_max_delay, t_bell)

                                # Loss induced from delay lines and treating two photon clicks as photon loss
                                if idx_eps == 1 and ix_eta == 0:
                                    loss = delay_loss(temp_max_delay, t_bell) + 0.0022129873919999998
                                elif idx_eps == 1 and ix_eta == 1:
                                    loss = delay_loss(temp_max_delay, t_bell) + 0.001478206422
                                elif idx_eps == 0 and ix_eta == 1:
                                    loss = delay_loss(temp_max_delay, t_bell) + 0.0003695516055
                                else:
                                    loss = delay_loss(temp_max_delay, t_bell) + 0.0005532468479999999

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
                                    # total_loss_with_multi = 1 - (1 - multiplexing_loss) * ((1 - switch_loss) ** (temporal_size))  # Sanity check where we only have spat mux and find max switch eff.
                                    total_loss_with_photon_as_erase = (1 - p_succ * (1 - loss) * ((1 - switch_loss) ** (temporal_size + 1)))  # Total loss when when only doing temporal mux
                                    if total_loss_with_photon_as_erase < threshold_loss:
                                        average_numb = 2 * number_of_photons  # Times two due to the spin encoding
                                        if average_numb < best_multiplex:
                                            best_multiplex = average_numb
                                    if total_loss_with_multi < threshold_loss:
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
    pair_source_prob_vs_multiplexing_parity_codes_spin_error_fixed_scan_switch_loss_new(t_bells, linestyles)