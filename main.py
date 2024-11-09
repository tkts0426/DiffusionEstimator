# Library
import numpy as np
import time
import warnings
warnings.simplefilter('ignore')

# My code:
from utils.Euler_Maruyama_app import Cham 
from utils.FE_nonSpot import fft_cal_fourier_est
from utils.Calculation_Integral import cal_analytical_int
from utils.FZ_Estimator import cal_FZ_est
from utils.RV_Estimator import cal_RV_est
from utils.FE_spot import cal_FE_day_est, cal_FE_whole_est
from utils.Plot_Functions import hist_plot, cum_dist_plot, approximation_plot, approximation_plot_with_var

# fix random seed 
np.random.seed(326)

def main():

    print("********************************************************************************")
    print("*             Diffusion estimation on N estimation points over ND days         *")
    print("*                by means of Realized Volatility type methods                  *")
    print("*                           and Fourier estimator                              *")
    print("*                                                                              *")
    print("*                                                                              *")
    print("********************************************************************************")
    n_traj = 1 # the number of paths
    ND = 250 # the number of days
    T = ND  # the total period
    n_obs_intraday = 100 # the number of trading times per one day
    N = n_obs_intraday * ND # the number of trading times
    D = 1 # serves for possible scattered sampling
    x = np.arange(0.04, 0.1, 0.0001) # candidates of inital value of price process

    print("... calculation")
    """
    Estimation of the diffusion coefficient using all intraday prices
    """
    L = int(N/2) # D = 1 => L/n = 1/2
    # Parameter of Chan: ---------------------------------------------
    alpha = 0.079 
    beta=0.093
    gamma=1.474
    eta=0.794
    r_0=0.065
    # -----------------------------------------------------------------
    
    # Generating process
    r_processes = np.zeros((n_traj, N+1))
    timestamps = np.zeros((n_traj, N+1))
    obs_processes = np.zeros((n_traj, N+1))
    non_noise_sigmas = np.zeros((n_traj, N+1))
    fourier_estimators = np.zeros((n_traj, len(x)))
    
    analytical_ints = np.zeros((n_traj, len(x)))
    denorm_ints = np.zeros((n_traj, len(x)))
    var_of_var = np.zeros((n_traj, len(x)))
    
    FZ_estimators = np.zeros((n_traj, len(x)))
    RV_estimators = np.zeros((n_traj, len(x)))
    FE_day_estimators = np.zeros((n_traj, len(x)))
    FE_whole_estimators = np.zeros((n_traj, len(x)))
    
    total_obs_prices = np.array([])
    for k in range(n_traj):
        r_process, timestamp, non_noise_sigma = Cham(T, N, r_0, alpha, beta, gamma, eta)

        r_processes[k] = r_process
        timestamps[k] = timestamp
        non_noise_sigmas[k] = non_noise_sigma

        # sigma_eps = 2 * np.std(np.diff(r_process)) # epsiron
        sigma_eps = 0 # microstructure noise effect
        noise = sigma_eps * np.random.randn(len(r_process))
        obs_process = r_process + noise

        obs_processes[k] = obs_process
        total_obs_prices = np.append(total_obs_prices, obs_process)
    
    n = len(r_process)

    # plot histgram of total obserbed process
    sorted_total_obs_prices = np.sort(total_obs_prices)
    hist_plot(total_obs_prices)
    cum_dist_plot(sorted_total_obs_prices, "obserbed prices")
    
    """
    Calculate the fourier estimator
    """
    h_1 = 3
    S = np.std(obs_processes, axis=1)
    # h_2 = h_1 * S * N**(-1/5) # h = n^{-1/5}
    h_2 = h_1 * S * N**(-1/3) # h = n^{-1/5}
    h_FZs = h_1 * S * N**(-1/4)


    for k in range(n_traj):

        h = h_2[k]
        h_FZ = h_FZs[k]
        
        time_str = time.time()
        fourier_estimator, denorm_int = fft_cal_fourier_est(obs_processes[k], x, n, h, L, T) # fft
        time_end = time.time()
        fourier_time = time_end - time_str
        analytical_int = cal_analytical_int(obs_processes[k], x, eta, gamma, n, h, T)
        time_str = time.time()
        FZ_estimator = cal_FZ_est(obs_processes[k], x, ND, n,  T, n_obs_intraday, h_FZ)
        time_end = time.time()
        FZ_time = time_end - time_str
        time_str = time.time()
        RV_estimator = cal_RV_est(obs_processes[k], x, ND, n , T, n_obs_intraday, h)
        time_end = time.time()
        RV_time = time_end - time_str
        time_str = time.time()
        FE_day_estimator = cal_FE_day_est(obs_processes[k], x, timestamp, ND, T, n, h_2[k], n_obs_intraday)
        time_end = time.time()
        FE_day_time = time_end - time_str
        time_str = time.time()
        FE_whole_estimator = cal_FE_whole_est(obs_processes[k], x, timestamp, ND, T, n, h_2[k], n_obs_intraday, N)
        time_end = time.time()
        FE_whole_time = time_end - time_str

        fourier_estimators[k] = fourier_estimator
        denorm_ints[k] = denorm_int
        analytical_ints[k] = analytical_int
        # var_of_var[k] = fourier_estimator**2
        var_of_var[k] = fourier_estimator**2 * denorm_int

        FZ_estimators[k] = FZ_estimator
        RV_estimators[k] = RV_estimator
        FE_day_estimators[k] = FE_day_estimator
        FE_whole_estimators[k] = FE_whole_estimator
        
        
    if n_traj == 1:
        fourier_est = fourier_estimators[0]
        var = var_of_var[0]

        FZ_est = FZ_estimators[0]
        RV_est = RV_estimators[0]
        FE_day_est = FE_day_estimators[0]
        FE_whole_est = FE_whole_estimators[0]
    else:
        fourier_est = np.mean(fourier_estimators, axis=0)
        var = (np.sum(var_of_var, axis=0) - n_traj*(fourier_est**2)) / (n_traj - 1) # Specimen Dispersion

        FZ_est = np.mean(FZ_estimators, axis=0)
        RV_est = np.mean(RV_estimators, axis=0)
        FE_day_est = np.mean(FE_day_estimators, axis=0)
        FE_whole_est = np.mean(FE_whole_estimators, axis=0)  

    fourier_mse_error = np.mean(np.abs(fourier_est - (eta*x**gamma)**2)**2)
    FZ_mse_error = np.mean(np.abs(FZ_est - (eta*x**gamma)**2)**2)
    RV_mse_error = np.mean(np.abs(RV_est - (eta*x**gamma)**2)**2)
    FE_day_mse_error = np.mean(np.abs(FE_day_est - (eta*x**gamma)**2)**2)
    FE_whole_mse_error = np.mean(np.abs(FE_whole_est - (eta*x**gamma)**2)**2)
    
    print(" The Error:")
    print(" ")
    print(f" - The Mean Squared Error of Fourier estimator : {fourier_mse_error}")
    print(f" - The calculation time of Fourier estimator : {fourier_time}")
    print(" ")
    print(f" - The Mean Squared Error of FZ estimator : {FZ_mse_error}")
    print(f" - The calculation time of FZ estimator : {FZ_time}")
    print(" ")
    print(f" - The Mean Squared Error of RV estimator : {RV_mse_error}")
    print(f" - The calculation time of RV estimator : {RV_time}")
    print(" ")
    print(f" - The Mean Squared Error of FE[day by day] : {FE_day_mse_error}")
    print(f" - The calculation time of FE[day by day] : {FE_day_time}")
    print(" ")
    print(f" - The Mean Squared Error of FE[all in all]: {FE_whole_mse_error}")
    print(f" - The calculation timeof FE[all in all]: {FE_whole_time}")

    
    # approximation_plot(x, [fourier_est, FZ_est, RV_est], ["fourier", "Florens-Zmirou", "Realied volatility"], (eta * x**gamma)**2)
    approximation_plot(
        x, 
        [fourier_est, FZ_est, RV_est, FE_day_est, FE_whole_est], 
        ["fourier", "Florens-Zmirou", "Realied volatility", "FE day by day", "FE all in all"], 
        (eta * x**gamma)**2,
        label="original_"
    )
    # # approximation_plot(x, [fourier_est, RV_est], ["fourier", "Realized vol"], (eta * x**gamma)**2)
    # # approximation_plot_with_var(x, [fourier_est], ["fourier"], (eta * x**gamma)**2, var, n_traj, label="asympVar_")
    # # approximation_plot(x, [FE_whole_est], ["FE_whole"], (eta * x**gamma)**2)
    # approximation_plot(x, [fourier_est], ["fourier estimator"],  (eta * x**gamma)**2, label="original_only_") # direct culculation
    
    
if __name__ == "__main__":
    main()

