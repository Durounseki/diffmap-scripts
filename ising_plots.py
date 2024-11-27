import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

def plot_epsilonCriteria(SGE, epsilons, iopt1, iopt2, w=400, h=350, ms=10, label=False):
    fig, ax = plt.subplots(figsize=(w / 100, h / 100))

    ax.plot(epsilons, SGE, color='black', linewidth=3)
    ax.scatter(epsilons[iopt1], SGE[iopt1], marker='^', color='C1', s=ms)
    ax.scatter(epsilons[iopt2], SGE[iopt2], marker='D', color='C2', s=ms)

    ax.set_yscale('log')
    ax.set_ylim(1e-3, 10)
    ax.set_xlim(0, 1)

    if not label:
        ax.set_xticklabels([])
        ax.set_yticklabels([])

    return fig

def plot_epsilonCriteria(SGE, epsilons, iopt, w=400, h=350, ms=10, label=False):  # Overloaded function
    fig, ax = plt.subplots(figsize=(w / 100, h / 100))

    ax.plot(epsilons, SGE, color='black', linewidth=3)
    ax.scatter(epsilons[iopt], SGE[iopt], marker='^', color='C1', s=ms)

    ax.set_yscale('log')
    ax.set_ylim(1e-3, 10)
    ax.set_xlim(0, 1)

    if not label:
        ax.set_xticklabels([])
        ax.set_yticklabels([])

    return fig

def plot_projection(X, mode_epsilon_low, mode_epsilon_high, w=400, h=350, ms=10, label=False):
    fig, ax = plt.subplots(figsize=(w / 100, h / 100))

    ax.scatter(mode_epsilon_low, np.mean(X, axis=0), color='C1', marker='^', s=ms)
    ax.scatter(mode_epsilon_high, np.mean(X, axis=0), color='C2', marker='D', s=ms)

    # Linear regression using sklearn
    reg = LinearRegression().fit(mode_epsilon_high.reshape(-1, 1), np.mean(X, axis=0))
    line_x = np.array([np.min(mode_epsilon_high), np.max(mode_epsilon_high)])
    line_y = reg.predict(line_x.reshape(-1, 1))
    ax.plot(line_x, line_y, color='black', linewidth=3)

    if not label:
        ax.set_xticklabels([])
        ax.set_yticklabels([])

    return fig

def plot_projection2(X, mode_epsilon_low, mode_epsilon_high, w=400, h=350, ms=10, label=False):
    fig, ax = plt.subplots(figsize=(w / 100, h / 100))

    ax.scatter(mode_epsilon_low, X, color='C1', marker='^', s=ms)
    ax.scatter(mode_epsilon_high, X, color='C2', marker='D', s=ms)

    # Linear regression for mode_epsilon_low
    reg_low = LinearRegression().fit(mode_epsilon_low.reshape(-1, 1), np.mean(X, axis=0))
    line_x_low = np.array([-0.01, 0.01])
    line_y_low = reg_low.predict(line_x_low.reshape(-1, 1))
    ax.plot(line_x_low, line_y_low, color='black', linewidth=3)

    # Linear regression for mode_epsilon_high
    reg_high = LinearRegression().fit(mode_epsilon_high.reshape(-1, 1), np.mean(X, axis=0))
    line_x_high = np.array([-0.01, 0.01])
    line_y_high = reg_high.predict(line_x_high.reshape(-1, 1))
    ax.plot(line_x_high, line_y_high, color='black', linestyle='--', dashes=(10, 30), linewidth=3)

    if not label:
        ax.set_xticklabels([])
        ax.set_yticklabels([])

    return fig

def plot_projection(X, mode, w=400, h=350, ms=10, label=False):  # Overloaded function
    fig, ax = plt.subplots(figsize=(w / 100, h / 100))

    ax.scatter(mode, np.mean(X, axis=0), color='C1', marker='^', s=ms)

    # Linear regression using sklearn
    reg = LinearRegression().fit(mode.reshape(-1, 1), np.mean(X, axis=0))
    line_x = np.array([np.min(mode), np.max(mode)])
    line_y = reg.predict(line_x.reshape(-1, 1))
    ax.plot(line_x, line_y, color='black', linewidth=3)

    if not label:
        ax.set_xticklabels([])
        ax.set_yticklabels([])

    return fig

def plot_projection(mode_epsilon_low_2, mode_epsilon_low_3, mode_epsilon_hi_2, mode_epsilon_hi_3, order_param, w=400, h=350, ms=10, label=False):  # Overloaded function
    fig, ax = plt.subplots(figsize=(w / 100, h / 100))

    # Scatter plot with colormap
    cm = plt.cm.get_cmap('roma')
    sc1 = ax.scatter(mode_epsilon_low_2, mode_epsilon_low_3, c=order_param, cmap=cm, marker='^', s=ms, alpha=0.25)
    sc2 = ax.scatter(mode_epsilon_hi_2, mode_epsilon_hi_3, c=order_param, cmap=cm, marker='D', s=ms, alpha=0.25)

    # Add colorbar
    plt.colorbar(sc1)

    if not label:
        ax.set_xticklabels([])
        ax.set_yticklabels([])

    return fig

def plot_projection3(mode_2, mode_3, order_param, w=400, h=350, ms=10, label=False):
    fig, ax = plt.subplots(figsize=(w / 100, h / 100))

    # Scatter plot with colormap
    cm = plt.cm.get_cmap('roma')
    sc = ax.scatter(mode_2, mode_3, c=order_param, cmap=cm, marker='^', s=ms, alpha=0.75)

    # Add colorbar
    plt.colorbar(sc)

    if not label:
        ax.set_xticklabels([])
        ax.set_yticklabels([])

    return fig

def plot_series(X, mode_epsilon_low, mode_epsilon_high, w=400, h=350, ms=10, label=False):
    fig, ax = plt.subplots(figsize=(w / 100, h / 100))

    js = np.arange(len(X))
    ax.plot(js, np.mean(X, axis=0), color='black', linewidth=2)

    # Linear regression using sklearn
    reg = LinearRegression().fit(mode_epsilon_high.reshape(-1, 1), np.mean(X, axis=0))

    ax.scatter(js[::20], reg.predict(mode_epsilon_high[::20].reshape(-1, 1)), color='C2', marker='D', s=ms)
    ax.scatter(js[::20], mode_epsilon_low[::20] * np.sqrt(len(js)), color='C1', marker='^', s=ms)

    if not label:
        ax.set_xticklabels([])
        ax.set_yticklabels([])

    return fig

def plot_series2(X, mode_epsilon_low, mode_epsilon_high, w=400, h=350, ms=10, label=False):
    fig, ax = plt.subplots(figsize=(w / 100, h / 100))

    js = np.arange(len(X))
    ax.plot(js, np.mean(X, axis=0), color='black', linewidth=2)

    # Linear regression for mode_epsilon_low
    reg_low = LinearRegression().fit(mode_epsilon_low.reshape(-1, 1), X)
    ax.scatter(js[::2], reg_low.predict(mode_epsilon_low[::2].reshape(-1, 1)), color='C1', marker='^', s=ms)

    # Linear regression for mode_epsilon_high
    reg_high = LinearRegression().fit(mode_epsilon_high.reshape(-1, 1), X)
    ax.scatter(js[::2], reg_high.predict(mode_epsilon_high[::2].reshape(-1, 1)), color='C2', marker='D', s=ms)

    if not label:
        ax.set_xticklabels([])
        ax.set_yticklabels([])

    return fig

def plot_series2(X, mode, w=400, h=350, ms=10, label=False):  # Overloaded function
    fig, ax = plt.subplots(figsize=(w / 100, h / 100))

    js = np.arange(len(X))
    ax.set_ylim(-1.05, 1.05)
    ax.plot(js, np.mean(X, axis=0), color='black', linewidth=1)

    # Linear regression using sklearn
    reg = LinearRegression().fit(mode.reshape(-1, 1), X)
    ax.scatter(js[::2], reg.predict(mode[::2].reshape(-1, 1)), color='C1', marker='^', s=ms)

    if not label:
        ax.set_xticklabels([])
        ax.set_yticklabels([])

    return fig

def plot_seriesNorm(X, mode_epsilon_low, mode_epsilon_high, w=400, h=350, ms=10, label=False):
    fig, ax = plt.subplots(figsize=(w / 100, h / 100))

    js = np.arange(len(X))
    ax.plot(js, np.mean(X, axis=0), color='black', linewidth=2)

    # Normalize mode_epsilon_low
    x1min = np.min(mode_epsilon_low)
    x1max = np.max(mode_epsilon_low)
    x1 = (2 * (mode_epsilon_low) - (x1max + x1min)) / (x1max - x1min)

    # Normalize mode_epsilon_high
    x2min = np.min(mode_epsilon_high)
    x2max = np.max(mode_epsilon_high)
    x2 = (2 * (mode_epsilon_high) - (x2max + x2min)) / (x2max - x2min)

    ax.scatter(js[::10], x1[::10], color='C1', marker='^', s=ms)
    ax.scatter(js[::10], x2[::10], color='C2', marker='D', s=ms)

    if not label:
        ax.set_xticklabels([])
        ax.set_yticklabels([])

    return fig

def plot_series(X, mode, w=400, h=350, ms=10, label=False):  # Overloaded function
    fig, ax = plt.subplots(figsize=(w / 100, h / 100))

    js = np.arange(len(X))
    ax.set_ylim(-1.05, 1.05)

    # Linear regression using sklearn
    reg = LinearRegression().fit(mode.reshape(-1, 1), np.mean(X, axis=0))

    ax.plot(js, np.mean(X, axis=0), color='black', linewidth=2)
    ax.scatter(js, reg.predict(mode.reshape(-1, 1)), color='C1', marker='^', s=ms)

    if not label:
        ax.set_xticklabels([])
        ax.set_yticklabels([])

    return fig

def plot_spectrum(evals_epsilon_low, evals_epsilon_high, w=400, h=350, ms=10, label=False):
    fig, ax = plt.subplots(figsize=(w / 100, h / 100))

    ax.scatter(evals_epsilon_low[:10], marker='^', color='C1', s=ms)
    ax.scatter(evals_epsilon_high[:10], marker='D', color='C2', s=ms)

    if not label:
        ax.set_xticklabels([])
        ax.set_yticklabels([])

    return fig

def plot_spectrum(evals, w=400, h=350, ms=10, label=False):  # Overloaded function
    fig, ax = plt.subplots(figsize=(w / 100, h / 100))

    ax.scatter(evals[:10], marker='^', color='C1', s=ms)

    if not label:
        ax.set_xticklabels([])
        ax.set_yticklabels([])

    return fig