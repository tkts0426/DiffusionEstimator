import numpy as np
import pandas as pd
import random

def _Cham_sigma(eta, gamma, pre_value):
    return eta * pre_value**gamma

def _Cham_mu(beta, alpha, pre_value):
    return beta * (alpha - pre_value)

def Cham(T, N, r_0, alpha, beta, gamma, eta):
    
    r_process = [r_0]
    timestamp = [0]
    dt = T/N
    for i in range(1, N+1):
        # random_term = random.normalvariate(0, np.sqrt(dt)) # normal distribution which mean is 0 and variance is dt
        random_term = random.gauss(0,1) * np.sqrt(dt)
        pre_r = r_process[i-1]

        r = pre_r + _Cham_sigma(eta, gamma, pre_r) * random_term + _Cham_mu(beta, alpha, pre_r) * dt
        if pd.isna(r) or pd.isna(_Cham_mu(beta, alpha, pre_r)) or pd.isna(_Cham_sigma(eta, gamma, pre_r)):
            print("nan is included in the process")
            print(f"  mu = {_Cham_mu(beta, alpha, pre_r)}")
            print(f"  sigma = {_Cham_sigma(eta, gamma, pre_r)}")
            print(r)
            break
        r_process.append(r)
        timestamp.append(timestamp[i-1] + dt)
    
    r_process = np.array(r_process)
    time_stamp = np.array(timestamp)

    non_noise_sigma = (eta * r_process**gamma)**2
    
    return r_process, time_stamp, non_noise_sigma

def _Sinh_sigma(x):
    return np.sqrt(1 + x**2)

def Sinh(T, N, x_0):
    
    sin_process = [x_0]
    timestamp = [0]
    dt = T / N
    for i in range(1, N+1):
        random_term = random.gauss(0, 1) * np.sqrt(dt)
        pre_x = sin_process[i-1]
        
        x = pre_x + _Sinh_sigma(pre_x) * random_term + (1/2) * pre_x * dt
        
        sin_process.append(x)
        timestamp.append(timestamp[i-1] + dt)

    sin_process = np.array(sin_process)
    timestamp = np.array(timestamp)

    non_noise_sigma = 1 + sin_process**2

    return sin_process, timestamp, non_noise_sigma
