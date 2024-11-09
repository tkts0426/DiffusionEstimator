import matplotlib.pyplot as plt
import numpy as np

def hist_plot(process, label=""):
    
    fig, ax = plt.subplots()
    
    ax.hist(process, bins=200, alpha=0.5)
    ax.set_title("observed data")

    fig.savefig(f"figure/{label}Distribution_of_obs_prices.png")

def cum_dist_plot(sorted_process, x_label=None, label=""):
    
    fig, ax = plt.subplots()
    
    n_samples = len(sorted_process)
    prob = [k / n_samples for k in range(n_samples)]
    ax.plot(sorted_process, prob)
    ax.axhline(0.5, ls = "-.", color = "gray")
    ax.axhline(0.95, ls = "-.", color = "gray")

    ax.set_title("Cumulative Distribution")
    ax.set_xlabel(x_label)
    ax.set_ylabel('probability')

    fig.savefig(f"figure/{label}Cumlative_distribution_of_obs_prices.png")

def approximation_plot(x, est_list, est_names, real, label=""):
    
    fig, ax = plt.subplots()
    
    for est, est_name in zip(est_list, est_names):
        ax.plot(x, est, label=est_name, linestyle="-") # plots of estimator
    
    ax.plot(x, real, color = "blue", label="sigma^2", linestyle="--")# plots of real function

    ax.set_title(f"Compare between estimaters and sigma^2 function")
    ax.set_xlabel("x")
    ax.set_ylabel("values of sigma^2 function")
    plt.legend()
    
    fig.savefig(f"figure/{label}plots of estimators and sigma^2 function.png", dpi=300)

def approximation_plot_with_var(x, est_list, est_names, real, var, n_traj, label=""):
    
    fig, ax = plt.subplots()
    
    for est, est_name in zip(est_list, est_names):
        ax.plot(x, est, label=est_name, linestyle="-") # plots of estimator
    
    ax.plot(x, real, color = "blue", label="sigma^2", linestyle="--")# plots of real function

    ax.plot(x, est[0] + 1.96 * np.sqrt(var / n_traj), color = "gray", linestyle="-", label="upper CI") # upper confidential interval
    ax.plot(x, est[0] - 1.96 * np.sqrt(var / n_traj), color = "gray", linestyle="-", label="lower CI") # lower confidential interval
    ax.set_title(f"Compare between estimaters and sigma^2 function")
    ax.set_xlabel("x")
    ax.set_ylabel("values of sigma^2 function")
    plt.legend()
    
    fig.savefig(f"figure/{label}plots of estimators and sigma^2 function.png", dpi=300)