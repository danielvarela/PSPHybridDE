import numpy as np

IP_ADDRESS = "10.8.0.6"


def convert_range(OldValue, old_range, new_range):
    OldMin = old_range[0]
    OldMax = old_range[1]
    NewMin = new_range[0]
    NewMax = new_range[1]
    NewValue = (((OldValue - OldMin) * (NewMax - NewMin)) / (OldMax - OldMin)) + NewMin
    return NewValue


def convert_translation(OldValue, trans_max_magnitude):
    OldMin = -1
    OldMax = 1
    NewMin = -1 * trans_max_magnitude
    NewMax = trans_max_magnitude
    NewValue = (((OldValue - OldMin) * (NewMax - NewMin)) / (OldMax - OldMin)) + NewMin
    return NewValue


def convert_gamma(OldValue):
    return ((OldValue - (-1)) * (180)) / (1 - (-1))


def random_individual(max_value=180):
    ind = []
    for i in range(6):
        if i < 1:
            # ind.append(0)
            ind.append(np.random.random_sample() * max_value)
        else:
            # ind.append(0)
            ind.append(np.random.random_sample() * 1)
    return ind
