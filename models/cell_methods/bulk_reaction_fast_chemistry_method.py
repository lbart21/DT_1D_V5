from numpy import array as array

def bulk_reaction_fast_chemsitry(massf_current, x, massf_goal):
    for ind, val in enumerate(massf_goal):
        if massf_current[ind] + x * (massf_goal[ind] - massf_current[ind]) < 0.0:
            x = min(x, massf_current[ind] / (massf_current[ind] - massf_goal[ind]))
    print("x", x)
    massf_new = (array(massf_current) + x * (array(massf_goal) - array(massf_current))).tolist()

    return massf_new