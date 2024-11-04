import numpy as np

def gaussian_kernel(x):
    # print(f"x^2 = {-x**2}")
    return np.exp(-x**2 / 2) / np.sqrt(2*np.pi)