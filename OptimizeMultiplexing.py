import json
from MiscPlotsPaper import *
from RepetionCodes import get_threshold_given_eps_in_optimize_func
from ErrorsFromHigherOrderPhotonComponents import effective_loss_rate, error_rate

                

def mux_plot_optimize_p(t_bells, linestyles, p_th=1):
    colors = ["red", "blue"]
    fig, ax = plt.subplots(figsize=(9, 6))

    # The success probability of heralding a resource state given spin error which has to be multiplexed to average to one
    succ_lists = [[0.7866335111212663, 0.7261559912925004, 0.6703280704865696, 0.6187922808189714, 0.5712186370521632, 0.527302523689994],
                  [0.7866335111212663, 0.7261559912925004, 0.6703280704865696, 0.6187922808189714, 0.5712186370521632, 0.527302523689994]]
    eps_spin = 0.01
    rep_codes_sizes = [3, 4, 5, 6, 7, 8]
    sizes_code = [3 * 6, 4 * 6, 5 * 6, 6 * 6, 7 * 6, 8 * 6]
    epses = [0.01 / 100, 0.05 / 100]
    eta_spin_gates = [0.8, 0.9]
    handles = []
    pair_gen_probs = [0.001 * x for x in range(50, 91)]
    save_dict = {0:[], 1:[], 2:[], 3:[], 4:[], 5:[], 6:[], 7:[], 8:[]}

    for idx_eps, eps in enumerate(epses):
        succ_list = succ_lists[idx_eps]
        t_bell = t_bells[idx_eps]
        for ix_eta, eta_spin_gate_in in enumerate(eta_spin_gates):
            if idx_eps == 0 and ix_eta == 0:
                dBS = np.linspace(0.0001, 0.003, 40)
            elif idx_eps == 0 and ix_eta == 1:
                dBS = np.linspace(0.0001, 0.0035, 40)
            elif idx_eps == 1 and ix_eta == 0:
                dBS = np.linspace(0.0001, 0.0012, 30)
            else:
                dBS = np.linspace(0.0001, 0.0014, 30)
            eta_spin_gate = eta_spin_gate_in * eta_spin_gate_in
            code_counts = {"0":0, "1":0, "2":0, "3":0, "4":0}
            multiplexing = []
            switch_loss_x_axis = []
            size_spat_mux = []
            chosen_p = []
            size_temp_mux = []
            spin_photon_errors_final = []
            mux_losses = []
            photon_succes = []
            thresholds_erase = []
            flag_stop_switch_loss = False
            for db in dBS:
                if flag_stop_switch_loss:
                    continue
                else:
                    switch_loss = 1 - 10 ** (-db / 10)
                    multiplex_one_instance = []
                    spat_mux_one_instance = []
                    best_idxs = []
                    p_one_instance = []
                    temp_one_instance = []
                    spin_photon_err_one_instance = []
                    mux_losses_one_instance = []
                    photon_succes_one_instance = []
                    th_one_instance = []
                    print("At DB: ", db)
                    for ix_r, rep_code_size in enumerate(rep_codes_sizes):
                        number_of_photons = sizes_code[ix_r]
                        best_multiplex = 10 ** 12
                        best_spat_mux = 0
                        best_p = 0
                        best_temp = 0
                        best_sp = 0
                        best_mux_loss = 0
                        best_photon_succ = 0
                        best_th = 0
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

                                        tot_prob = p_succ * (succ_list[ix_r]) # Multiply success probability with the error detection heralding probability

                                        # Find the binary size of a log tree for the spatial mux above
                                        # Want on average to get out 1 resource state, that is we need N*p = 1, giving N = 1 / p.

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

                                            p_succ_as_erase = (1 - (1 - pair_s * eta_spin_gate) ** temp_numb_input)
                                            total_loss_with_photon_as_erase = (1 - p_succ_as_erase * (1 - loss) * ((1 - switch_loss) ** (temporal_size + 1 + 1))) # An additional to del with the discarding of spin errors
                                            loss_from_detection = effective_loss_rate(pair_s, eta_spin_gate, 1 - total_loss_with_photon_as_erase)
                                            total_loss_with_photon_as_erase += loss_from_detection

                                            error_rate_spat = error_rate(pair_s, eta_spin_gate, 1 - total_loss_with_multi)
                                            error_rate_erase = error_rate(pair_s, eta_spin_gate, 1 - total_loss_with_photon_as_erase)

                                            spin_photon_error_in = 2 * eps * (2 / 3) * (1 - eps * (2 / 3)) # (error_rate_spat) * (2 / 3)
                                            spin_photon_error = spin_photon_error_in + error_rate_spat * (2 / 3) # 2 * spin_photon_error_in * (1 - spin_photon_error_in) + eps * (2 / 3)
                                            threshold_eta_spat = get_threshold_given_eps_in_optimize_func(eps_spin, rep_code_size ,spin_err=spin_photon_error)
                                            threshold_loss_spat = 1 - threshold_eta_spat

                                            # spin_photon_error_in = (error_rate_erase) * (2 / 3)
                                            # spin_photon_error_erase = 2 * spin_photon_error_in * (1 - spin_photon_error_in) + eps * (2 / 3)
                                            spin_photon_error_in = 2 * eps * (2 / 3) * (1 - eps * (2 / 3))  # The spin photon gate error is multiplied by two since we do the gate twice
                                            spin_photon_error_erase = spin_photon_error_in + error_rate_erase * (2 / 3)  # Error from two photon components are just added once
                                            threshold_eta_erase = get_threshold_given_eps_in_optimize_func(eps_spin, rep_code_size, spin_err=spin_photon_error_erase)
                                            threshold_loss_erase = 1 - threshold_eta_erase

                                            if total_loss_with_photon_as_erase < threshold_loss_erase and threshold_loss_erase != 2:
                                                average_numb = 2 * number_of_photons / succ_list[ix_r]  # Times two due to the spin encoding and divided by succ rate to deal with heralding of spin errors
                                                if average_numb < best_multiplex:
                                                    best_multiplex = average_numb
                                                    best_p = pair_s
                                                    best_temp = Ntemp
                                                    best_sp = spin_photon_error_erase
                                                    best_spat_mux = 0
                                                    best_mux_loss = total_loss_with_photon_as_erase
                                                    best_photon_succ = p_succ
                                                    best_th = threshold_loss_erase

                                            if total_loss_with_multi < threshold_loss_spat and threshold_loss_spat != 2:
                                                average_numb = 2 * number_of_photons * numb_spin_inputs  # Times two due to the spin encoding and the overhead from spatial mux
                                                if average_numb < best_multiplex:
                                                    best_multiplex = average_numb
                                                    best_spat_mux = numb_spin_inputs
                                                    best_p = pair_s
                                                    best_temp = Ntemp
                                                    best_sp = spin_photon_error
                                                    best_mux_loss = total_loss_with_multi
                                                    best_photon_succ = p_succ
                                                    best_th = threshold_loss_spat
                        if best_multiplex < 10 ** 6:
                            multiplex_one_instance.append(best_multiplex)
                            best_idxs.append(ix_r)
                            spat_mux_one_instance.append(best_spat_mux)
                            p_one_instance.append(best_p)
                            temp_one_instance.append(best_temp)
                            spin_photon_err_one_instance.append(best_sp)
                            mux_losses_one_instance.append(best_mux_loss)
                            photon_succes_one_instance.append(best_photon_succ)
                            th_one_instance.append(best_th)

                    if len(multiplex_one_instance) > 0:
                        multiplexing.append(min(multiplex_one_instance))
                        switch_loss_x_axis.append(db * 1000)
                        code_idx_idx = multiplex_one_instance.index(min(multiplex_one_instance))
                        code_size = best_idxs[code_idx_idx]
                        size_spat_mux.append(spat_mux_one_instance[code_idx_idx])
                        code_counts[str(code_size)] += 1
                        chosen_p.append(p_one_instance[code_idx_idx])
                        size_temp_mux.append(temp_one_instance[code_idx_idx])
                        spin_photon_errors_final.append(spin_photon_err_one_instance[code_idx_idx])
                        mux_losses.append(mux_losses_one_instance[code_idx_idx])
                        photon_succes.append(photon_succes_one_instance[code_idx_idx])
                        thresholds_erase.append(th_one_instance[code_idx_idx])
                    else:
                        flag_stop_switch_loss = True
            print("Spin photon gate eff. :", eta_spin_gate_in)
            print("Spin photon gate error: ", eps)
            print("size of spat mux: ", size_spat_mux)
            print("pair gens: ", chosen_p)
            print("size of temp mux: ", size_temp_mux)
            print("spin phootn errors: ", spin_photon_errors_final)
            print("code counts: ", code_counts)
            print("Number of required spins: ", multiplexing)
            print("Loss induced by mux network: ", mux_losses)
            print("probability of generating all photons: ", photon_succes)
            print("Photon loss thresholds: ", thresholds_erase)
            print()
            save_dict[0].append(size_spat_mux)
            save_dict[1].append(chosen_p)
            save_dict[2].append(size_temp_mux)
            save_dict[3].append(spin_photon_errors_final)
            save_dict[4].append(multiplexing)
            save_dict[5].append(mux_losses)
            save_dict[6].append(photon_succes)
            save_dict[7].append(thresholds_erase)
            save_dict[8].append(switch_loss_x_axis)
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
    # with open('data_0_8.json', 'w') as fp:
    #     json.dump(save_dict, fp)
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
    linestyles = ["solid", "dashed", "dashdot"]

    ###################################################################
    ####  PLOT MULTIPLEXING FOR PARITY CODES WITH SPIN ERROR RATES ####
    ###################################################################

    t_bells = [9.5 * (10 ** (-9)), 9.5 * (10 ** (-9))]
    mux_plot_optimize_p(t_bells, linestyles)


