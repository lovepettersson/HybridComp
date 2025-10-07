import numpy as np
import matplotlib.pyplot as plt
from RepetionCodes import *

def error_rate(p, eta_id, eta_s):
    # The effective error rate from two photon events
    error = (p ** 2) * 4 * eta_id * (1 - eta_id) * eta_s * (1 - eta_s)
    return error

def effective_loss_rate(p, eta_id, eta_s):
    # The effective loss rate from detecting more than two photons
    error = (p ** 2) * 2 * eta_id * (1 - eta_id) * eta_s * eta_s
    return error



def p_prob_threshold(p_vals, eta_thres, n_attempt, eta_id_in=0.75, eta_s=0.985, spin_code_flag=True):
    fitted_thres = erasure_meas_error_threshold_fit()

    if spin_code_flag:
        eta_id = eta_id_in ** 2
    else:
        eta_id = eta_id_in
    error_thres_list = []
    for p in p_vals:
        eps = error_rate(p, eta_id, eta_s)
        loss = effective_loss_rate(p, eta_id, eta_s)
        eta = (1 - loss) * eta_s
        if eta_thres < eta:
            erase = erasure_rate_rep_code(eta, n_attempt)
            if spin_code_flag:
                log_err = no_spin_code_error_rate(n_attempt, 2 * eps * (1 - eps), eta)
            else:
                log_err = no_spin_code_error_rate(n_attempt, eps, eta)
            if log_err > fitted_thres(erase):
                print("p: ", p, ", log err: ", log_err, ", erasure: ", erase, eta, ", phyiscal error: ", eps)
                error_thres_list.append(p)
                break
    return max(error_thres_list)


if __name__ == '__main__':
    eta_s = 0.99
    eta_id = 0.8 ** 2
    p_values = np.linspace(0.03, 0.04, 20)
    for p in p_values:
        error = error_rate(p, eta_id, eta_s)
        loss = effective_loss_rate(p, eta_id, eta_s)
        print("p=", p, ", error=", error, ", loss=", loss)
    p = 0.07
    error = error_rate(p, eta_id, eta_s)
    loss = effective_loss_rate(p, eta_id, eta_s)
    print("p=", p, ", error=", error, ", loss=", loss)
    p = 0.035
    error = error_rate(p, eta_id, eta_s)
    loss = effective_loss_rate(p, eta_id, eta_s)
    print("p=", p, ", error=", error, ", loss=", loss)
    ###############################################################
    ########## PLOTTING PAIR GEN PROBABILITY THRESHOLDS ###########
    ###############################################################
    '''
    eta_id = 0.8
    fig, ax = plt.subplots(figsize=(9, 6))
    ns = [4, 5, 6, 7, 8]
    eta_thres = [0.9777777777777777, 0.9781818181818182, 0.9797979797979798, 0.9818181818181818,
                 0.9834343434343434]
    p_values = np.linspace(10 ** (-2), 2 * 10 ** (-1), 500)
    threshold_p_values = []
    for idx in range(len(ns)):
        eta_th = eta_thres[idx]
        n_attempt = ns[idx]
        eta_s = eta_th + 0.007  # This is the loss rate for the signal photons, adding some buffer to be above the threshold
        threshold = p_prob_threshold(p_values, eta_th, n_attempt, eta_s=eta_s)
        print("Threshold p: ", threshold, "n_attempt: ", n_attempt)
        threshold_p_values.append(threshold)
    ax.plot(ns, threshold_p_values, "-o", linewidth=2.5, markersize=9.5, color="black", label="Spin code")

    threshold_p_values = []
    for idx in range(len(ns)):
        eta_th = eta_thres[idx]
        n_attempt = ns[idx]
        eta_s = eta_th + 0.007  # This is the loss rate for the signal photons, adding some buffer to be above the threshold
        threshold = p_prob_threshold(p_values, eta_th, n_attempt, eta_s=eta_s, spin_code_flag=False)
        print("Threshold p: ", threshold, "n_attempt: ", n_attempt)
        threshold_p_values.append(threshold)
    ax.plot(ns, threshold_p_values, "-o", linewidth=2.5, markersize=9.5, color="red", label="Bare")
    ax.tick_params(axis='both', which='major', labelsize=18)
    ax.legend(fontsize=18, bbox_to_anchor=(0.65, 0.64), loc="lower left")
    # fig.savefig("p_value_thresholds.pdf", dpi=600)
    plt.show()
    '''
