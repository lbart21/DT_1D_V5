import numpy as np

def generate_spatially_varying_time_step_scaling_for_finite_rate_chemistry(profile, n_cells):
    """
    profile = ["constant", time_step_scaling]
    profile = ["linear", time_step_scaling_left, time_step_scaling_right, x_pos_list, L]
    profile = ["tanh", time_step_scaling_left, time_step_scaling_right, x_pos_list, L, beta, x_mid]
    profile = ["custom", time_step_scaling_list]
    """
    if profile[0] == "constant":
        [_, time_step_scaling, x_pos_list, L] = profile
        time_step_scaling_list = [time_step_scaling for _ in range(n_cells)]

    elif profile[0] == "linear":
        [_, time_step_scaling_left, time_step_scaling_right, x_pos_list, L] = profile
        time_step_scaling_list = (time_step_scaling_left + \
                                  (time_step_scaling_right - time_step_scaling_left) \
                                                * np.array(x_pos_list) / L).tolist()

    elif profile[0] == "tanh":
        [_, time_step_scaling_left, time_step_scaling_right, beta, x_mid, x_pos_list, L] = profile
        time_step_scaling_list = (time_step_scaling_left + \
                                  0.5 * (time_step_scaling_right - time_step_scaling_left) \
                                * (1.0 + np.tanh(beta * (np.array(x_pos_list) - x_mid) / L)))\
                                    .tolist()

    elif profile[0] == "custom":
        [_, time_step_scaling_list, x_pos_list, L] = profile

    return time_step_scaling_list