import numpy as np
import math



def delay_loss(N, t_bin):
    # Delay loss using telecom photons and having a attenuation length of L_att = 20 km
    time = N * t_bin / 2   # Time is number of tries times the time per gate divided by two. This is the average time following a uniform distribution.
    speed = 2 * (10 ** 8)
    L = time * speed
    transmission = np.exp(-L / 20000)
    loss = 1 - transmission
    return loss




def estimate_opt_switch_boosting_parallel(p_init=0.026):
    # Amount of requried mux to generate a 4 qubit GHZ for boosting starting from 12 single photons that are
    # used to create two 3-qubit GHZ, that are then fused together into the final product
    best_temp_k = 0
    best_final_k = 0
    best_GHZ_k = 0
    best_final_N = 0
    best_tot_k = 10 ** 4
    for ghz_k in range(1, 14):
        N_ghz = 2 ** ghz_k
        GHZ_p = 1 - (1 - 1 / 32) ** N_ghz
        for k in range(1, 14):
            N = 2 ** k
            p = 1 - (1 - p_init) ** N
            prob_six_photons = p ** 12
            four_qbt_entangled_state = (GHZ_p ** 2)  # Need two GHZ states
            six_ent_state = 1 / 2  # Need one succesful fusions

            for k_final in range(1, 14):
                for N_final in range(2 ** (k_final - 1), 2 ** (k_final)):
                    tot_prob = four_qbt_entangled_state * prob_six_photons * six_ent_state
                    # tot_prob = prob_six_photons * six_ent_state
                    if tot_prob * N_final > 1:
                        tot_k = k + k_final + ghz_k
                        if tot_k < best_tot_k:
                            best_tot_k = tot_k
                            best_temp_k = k
                            best_final_k = k_final
                            best_GHZ_k = ghz_k
                            best_final_N = N_final


    print("Total number of required switches: ", best_tot_k, ", out of which ", best_temp_k,
          " are for photon generation and ", best_final_k, " are for fusing together 4 GHZ to the boosted entangled resource ", ", and ", best_GHZ_k, " are for GHZ generation")

    print(best_final_N, 2 ** best_final_k)

def estimate_opt_switch_ring_type_1(p_init=0.026):
    best_temp_k = 0
    best_ring_k = 0
    best_GHZ_k = 0
    best_ring_N = 0
    best_tot_k = 10 ** 4
    for ghz_k in range(1, 14):
        N_ghz = 2 ** ghz_k
        GHZ_p = 1 - (1 - 1 / 32) ** N_ghz
        for k in range(1, 14):
            N = 2 ** k
            p = 1 - (1 - p_init) ** N
            prob_six_photons = p ** 18  # I need 18 photons for three GHZ
            six_GHZ = (GHZ_p) ** 3  # Need six GHZ states
            ring_prob = (1 / 2) ** 3  # Need three type 1 fusions to create one ring from three GHZ

            for k_ring in range(1, 14):
                for N_ring in range(2 ** (k_ring - 1), 2 ** (k_ring)):
                    # N_ring = 2 ** k_ring
                    tot_prob = six_GHZ * prob_six_photons * ring_prob
                    if tot_prob * N_ring > 1:
                        tot_k = k + k_ring + ghz_k
                        if tot_k < best_tot_k:
                            best_tot_k = tot_k
                            best_temp_k = k
                            best_ring_k = k_ring
                            best_GHZ_k = ghz_k
                            best_ring_N = N_ring


    print("Total number of required switches: ", best_tot_k, ", out of which ", best_temp_k,
          " are for photon generation and ", best_ring_k, " are for ring generation ", ", and ", best_GHZ_k, " are for GHZ generation")
    print(best_ring_N, 2 ** best_ring_k)


def optimize_numb_of_photon_sources_with_switch_loss_type_1(numb_switches=17, tot_switch_boost=15, ratio_ring=62/64, ratio_boosting=63/64, db=0.0003):
    numb_entangled_states_boost = 6
    numb_photons_ring = 18
    numb_photons_boost = 12
    max_loss = 0.008  
    t_bin = 1 / (10 ** 9)
    switch_loss = 1 - 10 ** (-db / 10)
    best_overhead = 10 ** 12
    best_numb_temp = 0
    best_switch_loss = 0
    best_delay_loss = 0
    for numb_switch in range(1, numb_switches):
        N = (2 ** (numb_switch) - 1)
        loss_delay_line = delay_loss(N, t_bin)
        loss = 1 - (1 - switch_loss) ** (numb_switches + 1) * (1 - loss_delay_line)
        if loss < max_loss:
            overhead = numb_photons_ring * 2 ** (numb_switches - numb_switch)
            if overhead < best_overhead:
                best_overhead = overhead
                best_numb_temp = numb_switch
                best_switch_loss = 1 - (1 - switch_loss) ** (numb_switches + 1)
                best_delay_loss = loss_delay_line
    print("Number of photon sources overhead from spatial in generating six-ring: ", best_overhead, ", with ", best_numb_temp, " number of temp. mux")
    print("With adjusted ration for required number of inputs in spatial mux: ", best_overhead * ratio_ring)
    print("The amount of loss from switches is: ", best_switch_loss, ", and the delay loss is: ", best_delay_loss)
    overhead_ring = best_overhead * ratio_ring
    best_overhead = 10 ** 12
    best_switch_loss = 0
    best_delay_loss = 0
    for numb_switch in range(1, tot_switch_boost):
        N = (2 ** (numb_switch) - 1)
        loss_delay_line = delay_loss(N, t_bin)
        loss = 1 - (1 - switch_loss) ** (tot_switch_boost + 1) * (1 - loss_delay_line)
        if loss < max_loss:
            overhead = numb_photons_boost * (2 ** (tot_switch_boost - numb_switch))
            if overhead < best_overhead:
                best_overhead = overhead
                best_switch_loss = 1 - (1 - switch_loss) ** (tot_switch_boost + 1)
                best_delay_loss = loss_delay_line
    print("Number of photon sources overhead from spatial in boosted fusion entangled states: ",
            best_overhead)
    print("With adjusted ration for required number of inputs in spatial mux: ", best_overhead * ratio_boosting * numb_entangled_states_boost)
    print("The amount of loss from switches is: ", best_switch_loss, ", and the delay loss is: ", best_delay_loss)
    print("Total overhead: ", best_overhead * ratio_boosting * numb_entangled_states_boost + overhead_ring)

if __name__ == '__main__':
    p = 0.05
    db = 0.0015
    ratio_ring = 62 / 64
    ratio_boosted = 63 / 64
    numb_switches_ring = 17
    numb_switches_boost = 15
    estimate_opt_switch_ring_type_1(p)
    estimate_opt_switch_boosting_parallel(p)
    optimize_numb_of_photon_sources_with_switch_loss_type_1(numb_switches=numb_switches_ring, tot_switch_boost=numb_switches_boost,
                                                            ratio_ring=ratio_ring, ratio_boosting=ratio_boosted, db=db)

    print(2 ** 13 - 1, 2 ** 10 - 1)








