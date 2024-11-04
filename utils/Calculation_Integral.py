import numpy as np

from src.Functions import gaussian_kernel

def cal_analytical_int(process, x, eta, gamma, n, h, T):
    
    analytical_int = np.zeros(len(x))
    for i in range(len(x)):
        sigma_square_part = (eta * process[1:n]**gamma)**2
        kernel_part = (1 / h) * gaussian_kernel((process[1:n] - x[i]) / h)
        analytical_int[i] += np.sum(sigma_square_part * kernel_part) * (T / n)
    
    return analytical_int

def cal_denorm_int(process, x, n, h, T):
    
    denorm_int = np.zeros(len(x))
    for i in range(len(x)):
        kernel_part = (1 / h) * gaussian_kernel((process - x[i]) / h)
        denorm_int[i] += np.sum(kernel_part) * (T / n)

    return denorm_int