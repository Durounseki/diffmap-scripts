import numpy as np
import matplotlib.pyplot as plt

# ... (Previous functions) ...

def plot_transition_fixedL(Ns, L, etas, r_0, v_0, w=800, h=350, ms=10, label=False):
    fig, ax = plt.subplots(figsize=(w / 100, h / 100))
    markers = ['o', 's', 'D', 'H', '^', 'x', '*', 'p']

    jumps = []
    for i, N in enumerate(Ns):
        filename = f"./vicsekOrderN{N}L{L:.1f}r{r_0:.2f}v{v_0:.2f}.csv"
        order = np.loadtxt(filename, delimiter=",")
        threshold = (np.max(order) - np.min(order)) / 4
        jump_index = findSingleJump(order, threshold)
        jumps.append(etas[jump_index])
        ax.scatter(etas, order, marker=markers[i % len(markers)], s=ms, label=f"N={N}")
        ax.axvline(jumps[-1], linestyle='--', linewidth=2, color='black')

    ax.legend()

    if not label:
        ax.set_xticklabels([])
        ax.set_yticklabels([])

    return fig

def plot_scaling_fixedL(Ns, L, etas, r_0, v_0, w=800, h=350, ms=10, label=False):
    fig, ax = plt.subplots(figsize=(w / 100, h / 100))

    jumps = []
    for N in Ns:
        filename = f"./vicsekOrderN{N}L{L:.1f}r{r_0:.2f}v{v_0:.2f}.csv"
        order = np.loadtxt(filename, delimiter=",")
        threshold = (np.max(order) - np.min(order)) / 4
        jump_index = findSingleJump(order, threshold)
        jumps.append(etas[jump_index])

    ax.scatter(Ns, jumps, s=ms)

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(50, 2000)
    ax.set_ylim(np.pi / 4, 2 * np.pi + 1)
    ax.set_xticks([100, 400, 1600])
    ax.set_yticks([np.pi / 4, np.pi / 2, np.pi, 2 * np.pi], ["π/4", "π/2", "π", "2π"])

    if not label:
        ax.set_xticklabels([])
        ax.set_yticklabels([])

    return fig

def findSingleJump(X, threshold):
    for i in range(len(X) - 1):
        if abs(X[i + 1] - X[i]) > threshold:
            return i