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


