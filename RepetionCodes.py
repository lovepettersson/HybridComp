import numpy as np
import matplotlib.pyplot as plt
import math


def error_detection_parity_code(eps):
    # We detect an error if either of the qubits suffer one,
    # and we get an error if both suffer an error because this is not detected.
    eps_mod = 2 * eps / 3
    detect = 2 * eps_mod * (1 - eps_mod)
    error = eps_mod ** 2
    return detect, error


def binom_coeff(n, k):
    return math.factorial(n) / (math.factorial(k) * math.factorial(n-k))


def erasure_rate_rep_code(eta, n_attempt):
    # Erasure rates for the XX and ZZ parity in a logical fusion for a parity code
    p_fail = (1 / 2) * (eta ** 2)
    p_succ = (1 / 2) * (eta ** 2)
    XX = (p_fail + p_succ) ** n_attempt
    ZZ = 1 - (1 - p_succ) ** n_attempt
    erase = ((1 - XX) + (1 - ZZ)) / 2
    return erase



def logical_ZZ_error(eps, n_attempt, eta):
    # The logical ZZ error rate after using majority voting.
    p_succ = (1 / 2) * (eta ** 2)
    p_no_succ = (1 / 2) * (eta ** 2) + (1 - eta ** 2)
    eps_mod = 2 * (1 - eps) * eps
    log_error = 0
    tot_succ = 0
    for f_att in range(1, n_attempt + 1):
        if f_att == 1:
            tot_succ += binom_coeff(n_attempt, 1) * p_succ * (p_no_succ ** (n_attempt - 1))
            log_error += binom_coeff(n_attempt, 1) * p_succ * (p_no_succ ** (n_attempt - 1)) * eps_mod
        elif f_att == 2:
            tot_succ += binom_coeff(n_attempt, 2) * (p_succ ** 2) * (p_no_succ ** (n_attempt - 2)) * (1 - (2 * eps_mod * (1 - eps_mod)))  # If we get an error here, we treat it as an erasure as we have no confidence in the outcome
            # log_error += binom_coeff(n_attempt, 2) * (p_succ ** 2) * ((1 - p_succ) ** (n_attempt - 2)) * (2 * eps_mod * (1 - eps_mod))
        else:
            inner_error = 0
            prob = binom_coeff(n_attempt, f_att) * (p_succ ** f_att) * (p_no_succ ** (n_attempt - f_att))
            tot_succ += prob
            for k in range(int((f_att + 1) / 2), f_att + 1):  # int((4 + 1) / 2) = 2, and thus four success is killed by two errors as it should
                err = binom_coeff(f_att, k) * (eps_mod ** k) * ((1 - eps_mod) ** (f_att - k))
                inner_error += prob * err
            log_error += inner_error
    return log_error / tot_succ / 2  # Divide by two because we have randomized bias


def ZZ_error_erase(eps, n_attempt, eta):
    # The logical ZZ erasure rate from one error events in two succesful fusions
    p_succ = (1 / 2) * (eta ** 2)
    eps_mod = 2 * (1 - eps) * eps
    log_erase = binom_coeff(n_attempt, 2) * (p_succ ** 2) * ((1 - p_succ) ** (n_attempt - 2)) * (2 * eps_mod * (1 - eps_mod))
    return log_erase / 2  # Divide by two because we have randomized bias


def logical_XX_error(n_attempt, eps):
    # The logical XX error rate
    eps_mod = 2 * eps * (1 - eps)
    log_err = 0
    for k in range(1, n_attempt + 1):
        if k % 2 != 0:
            log_err += binom_coeff(n_attempt, k) * (eps_mod ** k) * ((1 - eps_mod) ** (n_attempt - k))
    return log_err / 2


def no_spin_code_error_rate(n_attempt, eps, eta):
    # The average error rate of the two parities when not using a spin code to mitigate spin errors
    XX_succ = (eta ** 2) ** n_attempt
    XX_error = logical_XX_error(n_attempt, eps) 
    ZZ_error = logical_ZZ_error(eps, n_attempt, eta)  # Already normalized by success probability
    return ZZ_error + XX_error  # Both errors are already divided by two, so it is an average


def total_spin_code_error_rate(n, eps, eta, spin_err=0):
    # The average logical error rate of the two parities when using a spin code
    XX_succ = eta ** (2 * n)
    detect, error = error_detection_parity_code(eps)
    eps_in = 2 * (eps / 3) * (1 - eps / 3)  # Either of the spin code qubits can have an error, and X error leak hence the eps/ 3 factor
    eps_photonic = eps_in * (1 - spin_err) + spin_err * (1 - eps_in)
    XX_tot_err = error * (1 - spin_err) + spin_err * (1 - error)  # the X errors do not affect the XX parities as this is just combination of physical XX outcomes
    tot_eps = (eps_photonic) * (1 - error) + error * (1 - eps_photonic)
    XX_err = logical_XX_error(n, XX_tot_err)
    ZZ_err = logical_ZZ_error(tot_eps, n, eta)  # Already normalized by success probability
    return XX_err + ZZ_err  # Both errors are already divided by two, so it is an average



def erasure_meas_error_threshold_fit():
    # Phenom. error and erasure threshold for six qubit ring found through Gaussian elim. and MWPM
    error = [0.01046551724137931, 0.008863269777109873, 0.007630952468689918, 0.006450396110946164, 0.0055472183812063924, 0.0046879143162049, 0.0036724137931034495, 0.002724367348559638, 0.001490979180657437, 6.841130395236939e-19]
    plot_erasure = [0, 0.01649657687835838, 0.029183754561433895, 0.038965517241379304, 0.04863020743628507, 0.058391048236745284, 0.06650477842165327, 0.07857774501399406, 0.08907416331555247, 0.11889655172413793]
    z = np.polyfit(plot_erasure, error, 3)
    p = np.poly1d(z)
    return p


def optimize_error_plot_detect(depolarizing_noise_vals, eta, fitted_thres):
    # Optimizing the parity code size when using a spin code and throwing away error detection events.
    ns = [3, 4, 5, 6, 7, 8]
    eta_thres = [0.9826262626262626, 0.9777777777777777, 0.9781818181818182, 0.9797979797979798, 0.9818181818181818,
                  0.9834343434343434]
    error_thres_list = []
    n_succ = []
    for idx, n_attempt in enumerate(ns):
        for eps in depolarizing_noise_vals:
            if eta_thres[idx] < eta:
                # erase = erasure_rate_rep_code(eta, n_attempt)
                erase = erasure_rate_rep_code(eta, n_attempt) + ZZ_error_erase(2 * (eps / 3) * (1 - eps / 3), n_attempt, eta)  # Treating an error event when we get two successful fusions as an erasure
                log_err = total_spin_code_error_rate(n_attempt, eps, eta)
                if log_err > fitted_thres(erase):
                    error_thres_list.append(eps)
                    n_succ.append(n_attempt)
                    break
    idx_max = error_thres_list.index(max(error_thres_list))
    return max(error_thres_list), n_succ[idx_max]


def optimize_error_plot_erase(depolarizing_noise_vals, eta, fitted_thres):
    # Optimizing the parity code size when using a spin code and treating error detection events as erasures.
    ns = [3, 4, 5, 6, 7, 8]
    eta_thres = [0.9826262626262626, 0.9777777777777777, 0.9781818181818182, 0.9797979797979798, 0.9818181818181818,
                  0.9834343434343434]
    error_thres_list = []
    for idx, n_attempt in enumerate(ns):
        for eps in depolarizing_noise_vals:
            if eta_thres[idx] < eta:
                detect, error = error_detection_parity_code(eps)
                total_eta = (1 - detect) * eta
                erase = erasure_rate_rep_code(total_eta, n_attempt) + ZZ_error_erase(2 * (eps / 3) * (1 - eps / 3), n_attempt, total_eta)
                log_err = total_spin_code_error_rate(n_attempt, eps, total_eta)
                if log_err > fitted_thres(erase):
                    error_thres_list.append(eps)
                    break
    return max(error_thres_list)



def optimize_error_plot_bare(depolarizing_noise_vals, eta, fitted_thres):
    # Optimizing the parity code size when not using a spin code.
    ns = [3, 4, 5, 6, 7, 8]
    eta_thres = [0.9826262626262626, 0.9777777777777777, 0.9781818181818182, 0.9797979797979798, 0.9818181818181818,
                  0.9834343434343434]
    error_thres_list = []
    for idx, n_attempt in enumerate(ns):
        for eps in depolarizing_noise_vals:
            if eta_thres[idx] < eta:
                erase = erasure_rate_rep_code(eta, n_attempt) + ZZ_error_erase(2 * eps / 3, n_attempt, eta)
                log_err = no_spin_code_error_rate(n_attempt, 2 * eps / 3, eta)
                if log_err > fitted_thres(erase):
                    error_thres_list.append(eps)
                    break
    return max(error_thres_list)


def optimize_error_plot_code_spin_err(depolarizing_noise_vals, eta, fitted_thres):
    # Optimizing the parity code size when using a spin code for photonic link errors.
    ns = [3, 4, 5, 6, 7, 8]
    eta_thres = [0.9826262626262626, 0.9777777777777777, 0.9781818181818182, 0.9797979797979798, 0.9818181818181818,
                  0.9834343434343434]
    error_thres_list = []
    for idx, n_attempt in enumerate(ns):
        for eps in depolarizing_noise_vals:
            if eta_thres[idx] < eta:
                tot_err = 2 * (2 * eps / 3) * (1 - 2 * eps / 3)  # Two out of the three Paulis leaed to errors, and we have two spin code qubits
                erase = erasure_rate_rep_code(eta, n_attempt) + ZZ_error_erase(tot_err, n_attempt, eta)
                log_err = no_spin_code_error_rate(n_attempt, tot_err, eta)
                if log_err > fitted_thres(erase):
                    error_thres_list.append(eps)
                    break
    return max(error_thres_list)







def five_qbt_parity():
    # Five qubit code example
    fitted_thres = erasure_meas_error_threshold_fit()  # x is erasure
    eta_values = np.linspace(0.9781818181818182, 1, 10)
    depolarizing_noise_vals = np.linspace(0, 0.05, 300)
    log_error_thres = [0]
    n_attempt = 5
    for eta in eta_values[1:]:
        # erase = erasure_rate_rep_code(eta, n_attempt)
        for eps in depolarizing_noise_vals:
            erase = erasure_rate_rep_code(eta, n_attempt) + ZZ_error_erase(2 * (eps / 3) * (1 - eps / 3), n_attempt, eta)
            log_err = total_spin_code_error_rate(n_attempt, eps, eta)
            if log_err > fitted_thres(erase):
                log_error_thres.append(eps)
                break

    depolarizing_noise_vals = np.linspace(0, 0.02, 300)
    log_error_thres_erase = [0]
    for eta in eta_values[1:]:
        for eps in depolarizing_noise_vals:
            detect, _ = error_detection_parity_code(eps)
            total_eta = (1 - detect) * eta
            erase = erasure_rate_rep_code(total_eta, n_attempt) + ZZ_error_erase(2 * (eps / 3) * (1 - eps / 3), n_attempt, total_eta)
            log_err = total_spin_code_error_rate(n_attempt, eps, total_eta)
            if log_err > fitted_thres(erase):
                log_error_thres_erase.append(eps)
                break



    depolarizing_noise_vals = np.linspace(0, 0.004, 300)
    no_code_threhold = [0]
    for eta in eta_values[1:]:
        # erase = erasure_rate_rep_code(eta, n_attempt)
        for eps in depolarizing_noise_vals:
            erase = erasure_rate_rep_code(eta, n_attempt) + ZZ_error_erase(2 * eps / 3, n_attempt, eta)
            log_err = no_spin_code_error_rate(n_attempt, 2 * eps / 3, eta)
            if log_err > fitted_thres(erase):
                no_code_threhold.append(eps)
                break
    return eta_values, log_error_thres, log_error_thres_erase, no_code_threhold


def optimized_parity():
    # Optimizing the parity code size
    fitted_thres = erasure_meas_error_threshold_fit() # x is erasure

    eta_values = np.linspace(0.9781818181818182, 1, 10)
    loss_values = [1 - 0.9777777777777777]
    for eta in eta_values[1:]:
        loss_values.append(1 - eta)

    depolarizing_noise_vals = np.linspace(0, 0.05, 300)
    optimized_err_thres = [0]
    for eta in eta_values[1:]:
        err,_ = optimize_error_plot_detect(depolarizing_noise_vals, eta, fitted_thres)
        optimized_err_thres.append(err)

    depolarizing_noise_vals = np.linspace(0, 0.02, 300)
    optimized_thres_erase = [0]
    for eta in eta_values[1:]:
        err = optimize_error_plot_erase(depolarizing_noise_vals, eta, fitted_thres)
        optimized_thres_erase.append(err)



    depolarizing_noise_vals = np.linspace(0, 0.004, 300)
    no_code_threhold_optimized = [0]
    for eta in eta_values[1:]:
        err = optimize_error_plot_bare(depolarizing_noise_vals, eta, fitted_thres)
        no_code_threhold_optimized.append(err)
    return loss_values, optimized_err_thres, optimized_thres_erase, no_code_threhold_optimized


def optimized_parity_and_five_qbt_spin_photon_error():
    # Optimizing both five qubit code and over all parity codes for photon link errors
    fitted_thres = erasure_meas_error_threshold_fit() # x is erasure

    eta_values = np.linspace(0.9781818181818182, 1, 10)
    loss_values_no_opt = [1 - x for x in eta_values]
    loss_values_opt = [1 - 0.9777777777777777]
    for eta in eta_values[1:]:
        loss_values_opt.append(1 - eta)




    depolarizing_noise_vals = np.linspace(0, 0.004, 300)
    no_code_threhold_optimized = [0]
    for eta in eta_values[1:]:
        err = optimize_error_plot_bare(depolarizing_noise_vals, eta, fitted_thres)
        no_code_threhold_optimized.append(err)

    optimized_err_thres = [0]
    for eta in eta_values[1:]:
        err = optimize_error_plot_code_spin_err(depolarizing_noise_vals, eta, fitted_thres)
        optimized_err_thres.append(err)

    depolarizing_noise_vals = np.linspace(0, 0.004, 300)
    no_code_threhold = [0]
    n_attempt = 5
    for eta in eta_values[1:]:
        # erase = erasure_rate_rep_code(eta, n_attempt)
        for eps in depolarizing_noise_vals:
            erase = erasure_rate_rep_code(eta, n_attempt) + ZZ_error_erase(2 * eps / 3, n_attempt, eta)
            log_err = no_spin_code_error_rate(n_attempt, 2 * eps / 3, eta)
            if log_err > fitted_thres(erase):
                no_code_threhold.append(eps)
                break
    code_err_thres = [0]
    for eta in eta_values[1:]:
        # erase = erasure_rate_rep_code(eta, n_attempt)
        for eps in depolarizing_noise_vals:
            tot_err = 2 * (2 * eps / 3) * (1 - 2 * eps / 3)
            erase = erasure_rate_rep_code(eta, n_attempt) + ZZ_error_erase(tot_err, n_attempt, eta)
            log_err = no_spin_code_error_rate(n_attempt, tot_err, eta)
            if log_err > fitted_thres(erase):
                code_err_thres.append(eps)
                break



    return loss_values_opt ,loss_values_no_opt, optimized_err_thres, no_code_threhold_optimized, no_code_threhold, code_err_thres



def plot_thresholds_five_qbt_parity():
    # Plotting five_qbt example
    colors = plt.cm.gist_heat(np.linspace(0.9, 0, 3))
    eta_values, log_error_thres, log_error_thres_erase, no_code_threhold = five_qbt_parity()


    fig, ax = plt.subplots(figsize=(9, 6))
    loss_values = [1 - x for x in eta_values]
    # ax.plot(optimized_err_thres, loss_values, "-o", color=colors[3], linewidth=2.5,
    #         markersize=9.5, label="Optimized")
    ax.plot(log_error_thres, loss_values, "-o", color=colors[2], linewidth=2.5,
                                    markersize=9.5, label="Detect")

    ax.plot(log_error_thres_erase, loss_values, "-o", color=colors[1], linewidth=2.5,
                                    markersize=9.5, label="Erase")
    ax.plot(no_code_threhold, loss_values, "-o", color=colors[0], linewidth=2.5,
            markersize=9.5, label="Bare")

    ax.legend(fontsize=18)
    ax.tick_params(axis='both', which='major', labelsize=20)
    plt.show()


def plot_thresholds_full_set_parity():
    # Plotting both five qubit code and optimized codes
    colors = plt.cm.gist_heat(np.linspace(0.9, 0, 3))
    eta_values, log_error_thres, log_error_thres_erase, no_code_threhold = five_qbt_parity()
    loss_values_optimized, optimized_err_thres, optimized_thres_erase, no_code_threhold_optimized = optimized_parity()

    fig, ax = plt.subplots(figsize=(9, 6))
    loss_values = [1 - x for x in eta_values]
    handles = []


    ax.plot(optimized_err_thres, loss_values_optimized, "--D", color=colors[2], linewidth=2.5,
            markersize=9.5)

    ax.plot(optimized_thres_erase, loss_values_optimized, "--D", color=colors[1], linewidth=2.5,
            markersize=9.5)

    ax.plot(no_code_threhold_optimized, loss_values_optimized, "--D", color=colors[0], linewidth=2.5,
            markersize=9.5)





    line, = ax.plot(log_error_thres, loss_values, "-o", color=colors[2], linewidth=2.5,
                                    markersize=9.5, label="Detect")
    handles.append(line)
    line, = ax.plot(log_error_thres_erase, loss_values, "-o", color=colors[1], linewidth=2.5,
                                    markersize=9.5, label="Erase")
    handles.append(line)
    line, = ax.plot(no_code_threhold, loss_values, "-o", color=colors[0], linewidth=2.5,
            markersize=9.5, label="Bare")

    handles.append(line)
    line1, = ax.plot(0, -2, "-o", color="black", linewidth=2.5, markersize=9.5,
                     label="Five qbt code")
    line2, = ax.plot(0, -2, "--D", color="black", linewidth=2.5, markersize=9.5,
                     label="Optimized")
    first_legend = ax.legend(handles=[line1, line2], title="Linestyles", bbox_to_anchor=(0.68, 0.39),
                             loc="lower left",
                             fontsize=14, title_fontsize=14)  # 0.6, 0.25, 0.57, 0.2, slow:0.57, 0.09

    ax.add_artist(first_legend)

    ax.legend(handles=handles, title="Color", bbox_to_anchor=(0.78, 0.64), loc="lower left", fontsize=14,
              title_fontsize=14)

    ax.set_ylim((-0.0005, 0.024))
    ax.tick_params(axis='both', which='major', labelsize=20)
    # fig.savefig("rep_codes_spin_error.pdf", dpi=600)
    plt.show()





def plot_resource_succ_prob():
    # Plotting the success probability of generating a resource state when throwing away detection clicks
    colors = plt.cm.gist_heat(np.linspace(0.9, 0, 3))
    fitted_thres = erasure_meas_error_threshold_fit()

    eta_values = np.linspace(0.9781818181818182, 1, 10)
    loss_values_opt = [1 - 0.9777777777777777]
    for eta in eta_values[1:]:
        loss_values_opt.append(1 - eta)

    error_plot_opt = [0]
    depolarizing_noise_vals = np.linspace(0, 0.05, 300)
    plot_succ_prob_opt = [1]
    for eta in eta_values[1:]:
        err, n_size = optimize_error_plot_detect(depolarizing_noise_vals, eta, fitted_thres)
        detect, _ = error_detection_parity_code(err)
        plot_succ_prob_opt.append((1 - detect) ** (6 * n_size))
        error_plot_opt.append(err)

    error_plot = [0]
    eta_values = np.linspace(0.9781818181818182, 1, 10)
    depolarizing_noise_vals = np.linspace(0, 0.05, 300)
    plot_succ_prob = [1]
    n_attempt = 5
    for eta in eta_values[1:]:
        # erase = erasure_rate_rep_code(eta, n_attempt)
        for eps in depolarizing_noise_vals:
            erase = erasure_rate_rep_code(eta, n_attempt) + ZZ_error_erase(2 * (eps / 3) * (1 - eps / 3), n_attempt,
                                                                           eta)
            log_err = total_spin_code_error_rate(n_attempt, eps, eta)
            if log_err > fitted_thres(erase):
                detect, _ = error_detection_parity_code(eps)
                plot_succ_prob.append((1 - detect) ** (6 * n_attempt))
                error_plot.append(eps)
                break


    fig, ax = plt.subplots(figsize=(9, 6))
    handles = []



    print(error_plot_opt)
    line, = ax.plot(error_plot_opt, plot_succ_prob_opt, "--D", color=colors[2], linewidth=2.5,
            markersize=9.5, label="Optimized")
    handles.append(line)
    line, = ax.plot(error_plot, plot_succ_prob, "-o", color=colors[1], linewidth=2.5,
            markersize=9.5, label="Five qbt code")
    handles.append(line)



    ax.legend(handles=handles, title="Color", bbox_to_anchor=(0.68, 0.64), loc="lower left", fontsize=14,
              title_fontsize=14)

    ax.set_yscale("log")
    ax.tick_params(axis='both', which='major', labelsize=20)
    # fig.savefig("rep_codes_resource_state_succ_prob.pdf", dpi=600)
    plt.show()



def plot_spin_photon_gate_errors():
    # Plotting photon link error thresholds
    colors = plt.cm.gist_heat(np.linspace(0.9, 0, 3))
    loss_values_opt, loss_values_no_opt, optimized_err_thres, no_code_threhold_optimized, no_code_threhold, code_err_thres = optimized_parity_and_five_qbt_spin_photon_error()

    fig, ax = plt.subplots(figsize=(9, 6))
    handles = []

    ax.plot(optimized_err_thres, loss_values_opt, "--D", color=colors[2], linewidth=2.5,
            markersize=9.5)

    ax.plot(no_code_threhold_optimized, loss_values_opt, "--D", color=colors[0], linewidth=2.5,
            markersize=9.5)


    line, = ax.plot(code_err_thres, loss_values_no_opt, "-o", color=colors[2], linewidth=2.5,
                    markersize=9.5, label="Spin code")
    handles.append(line)
    line, = ax.plot(no_code_threhold, loss_values_no_opt, "-o", color=colors[0], linewidth=2.5,
                    markersize=9.5, label="Bare")

    handles.append(line)
    line1, = ax.plot(0, -2, "-o", color="black", linewidth=2.5, markersize=9.5,
                     label="Five qbt code")
    line2, = ax.plot(0, -2, "--D", color="black", linewidth=2.5, markersize=9.5,
                     label="Optimized")
    first_legend = ax.legend(handles=[line1, line2], title="Linestyles", bbox_to_anchor=(0.68, 0.39),
                             loc="lower left",
                             fontsize=14, title_fontsize=14)

    ax.add_artist(first_legend)

    ax.legend(handles=handles, title="Color", bbox_to_anchor=(0.72, 0.64), loc="lower left", fontsize=14,
              title_fontsize=14)

    ax.set_ylim((-0.0005, 0.024))
    ax.tick_params(axis='both', which='major', labelsize=20)
    # fig.savefig("rep_codes_spin_photon_error.pdf", dpi=600)
    plt.show()




def get_threshold_given_eps(eps, spin_err=0):
    # Get loss threshold given a spin error rate eps for all parity code sizes
    fitted_thres = erasure_meas_error_threshold_fit()
    ns = [3, 4, 5, 6, 7, 8]
    log_error_thres = []
    eta_values = np.linspace(0.98, 1, 200)
    succ_sizes = []
    for idx, n_attempt in enumerate(ns):
        for eta in eta_values[1:]:
            # erase = erasure_rate_rep_code(eta, n_attempt)
            err_erase_one = 2 * (eps / 3) * (1 - eps / 3)
            err_erase_in = err_erase_one * (1 - spin_err) + (1 - err_erase_one) * spin_err
            erase = erasure_rate_rep_code(eta, n_attempt) + ZZ_error_erase(err_erase_in, n_attempt, eta)
            log_err = total_spin_code_error_rate(n_attempt, eps, eta, spin_err=spin_err)  # logical_error_rate(eps / 3, n_attempt, eta)
            if log_err < fitted_thres(erase):
                log_error_thres.append(eta)
                succ_sizes.append(n_attempt)
                break
    print(succ_sizes)
    return log_error_thres



def get_threshold_given_eps_in_optimize_func(eps, n_attempt, spin_err=0):
    # Get loss threshold given a spin error rate eps for all parity code sizes
    fitted_thres = erasure_meas_error_threshold_fit()
    ns = [3, 4, 5, 6, 7, 8]
    log_error_thres = []
    eta_values = np.linspace(0.982, 1, 100)
    succ_sizes = []
    flag = False
    for eta in eta_values[1:]:
        # erase = erasure_rate_rep_code(eta, n_attempt)
        err_erase_one = 2 * (eps / 3) * (1 - eps / 3)
        err_erase_in = err_erase_one * (1 - spin_err) + (1 - err_erase_one) * spin_err
        erase = erasure_rate_rep_code(eta, n_attempt) + ZZ_error_erase(err_erase_in, n_attempt, eta)
        log_err = total_spin_code_error_rate(n_attempt, eps, eta, spin_err=spin_err)  # logical_error_rate(eps / 3, n_attempt, eta)
        if log_err < fitted_thres(erase):
            return eta
    if not flag:
        return -1

def get_threshold_given_eps_one_code(eps, n_attempt, spin_err=0):
    # Get loss threshold given a spin error rate eps for a specific parity code size
    fitted_thres = erasure_meas_error_threshold_fit()
    log_error_thres = []
    eta_values = np.linspace(0.98, 1, 100)
    for eta in eta_values[1:]:
        err_erase_one = 2 * (eps / 3) * (1 - eps / 3)
        err_erase_in = err_erase_one * (1 - spin_err) + (1 - err_erase_one) * spin_err
        erase = erasure_rate_rep_code(eta, n_attempt) + ZZ_error_erase(err_erase_in, n_attempt, eta)
        log_err = total_spin_code_error_rate(n_attempt, eps, eta, spin_err=spin_err)  # logical_error_rate(eps / 3, n_attempt, eta)
        if log_err < fitted_thres(erase):
            log_error_thres.append(eta)
            break
    return eta


def get_resource_state_overhead(eps):
    # Get the resource state overhead given a spin error rate eps for all parity code size
    detect, err = error_detection_parity_code(eps)
    resource_state_sizes = [3 * 6, 4 * 6, 5 * 6, 6 * 6, 7 * 6, 8 * 6]
    succ_list = []
    for size in resource_state_sizes:
        succ_list.append((1 - detect) ** size)
    return succ_list


if __name__ == '__main__':

    # plot_thresholds_five_qbt_parity()
    # plot_thresholds_full_set_parity()
    # plot_resource_succ_prob()
    # plot_spin_photon_gate_errors()

    ##########################################
    ######## LOSS THRESHOLD GIVEN EPS ########
    ##########################################

    eps = 0.01
    # spin_err_in = (10 ** (-5)) * 2 / 3
    spin_err_in = 10 * (10 ** (-4)) * 2 / 3
    spin_err = 2 * spin_err_in * (1 - spin_err_in)
    print(get_threshold_given_eps(eps, spin_err=spin_err))
    print(get_resource_state_overhead(eps))
    print(get_threshold_given_eps_in_optimize_func(eps, 5 ,spin_err=spin_err))

