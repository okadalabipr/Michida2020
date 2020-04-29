import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm


def model(N, Km, K1, K2, NFkB):

    return NFkB**N/(NFkB**N+K1*(Km**N+NFkB**N))*K2


def my_rc_params():
    plt.rcParams['font.size'] = 32
    plt.rcParams['axes.linewidth'] = 2
    plt.rcParams['xtick.major.width'] = 2
    plt.rcParams['ytick.major.width'] = 2
    plt.rcParams['lines.linewidth'] = 2
    plt.rcParams['lines.markersize'] = 12


if __name__ == '__main__':
    os.makedirs('./Fig3F', exist_ok=True)

    NFkB = np.linspace(0, 100, 1000)

    # Super Enhancer
    ## N
    my_rc_params()
    plt.figure(figsize=(8, 6))
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    
    Km = 53.83992218
    K1 = 3.379233326
    K2 = 7.973022134

    strength = np.linspace(0.5, 2, 30)

    for i, j in enumerate(strength):
        N = 4.223431997 * j

        plt.plot(
            NFkB, model(N, Km, K1, K2, NFkB), color=cm.viridis(i / len(strength))
        )

    N = 4.223431997
    plt.plot(NFkB, model(N, Km, K1, K2, NFkB), 'k', lw=3)

    plt.xlabel('Nuclear NF-κB level (A.U.)')
    plt.xticks([0, 25, 50, 75, 100])
    plt.ylabel('Fold change in\nmRNA level - 1')
    plt.title('')
    plt.savefig('./Fig3F/SE_N.pdf', bbox_inches='tight')
    plt.close()

    ## Km
    my_rc_params()
    plt.figure(figsize=(8, 6))
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)

    N = 4.223431997
    K1 = 3.379233326
    K2 = 7.973022134

    strength = np.linspace(0.5, 2, 30)

    for i, j in enumerate(strength):
        Km = 53.83992218 * j

        plt.plot(
            NFkB, model(N, Km, K1, K2, NFkB), color=cm.viridis(i / len(strength))
        )

    Km = 53.83992218
    plt.plot(NFkB, model(N, Km, K1, K2, NFkB), 'k', lw=3)

    plt.xlabel('Nuclear NF-κB level (A.U.)')
    plt.xticks([0, 25, 50, 75, 100])
    plt.ylabel('Fold change in\nmRNA level - 1')
    plt.title('')
    plt.savefig('./Fig3F/SE_Km.pdf', bbox_inches='tight')
    plt.close()

    # Typical Enhancer
    ## N
    my_rc_params()
    plt.figure(figsize=(8, 6))
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)

    Km = 242.4179786
    K1 = 2.16228223
    K2 = 10.40479864

    strength = np.linspace(0.5, 2, 30)

    for i, j in enumerate(strength):
        N = 1.300365511 * j

        plt.plot(
            NFkB, model(N, Km, K1, K2, NFkB), color=cm.viridis(i / len(strength))
        )

    N = 1.300365511
    plt.plot(NFkB, model(N, Km, K1, K2, NFkB), 'k', lw=3)

    plt.xlabel('Nuclear NF-κB level (A.U.)')
    plt.xticks([0, 25, 50, 75, 100])
    plt.ylabel('Fold change in\nmRNA level - 1')
    plt.title('')
    plt.savefig('./Fig3F/TE_N.pdf', bbox_inches='tight')
    plt.close()

    ## Km
    my_rc_params()
    plt.figure(figsize=(8, 6))
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)

    N = 1.300365511
    K1 = 2.16228223
    K2 = 10.40479864

    strength = np.linspace(0.5, 2, 30)

    for i, j in enumerate(strength):
        Km = 242.4179786 * j

        plt.plot(
            NFkB, model(N, Km, K1, K2, NFkB), color=cm.viridis(i / len(strength))
        )

    Km = 242.4179786
    plt.plot(NFkB, model(N, Km, K1, K2, NFkB), 'k', lw=3)

    plt.xlabel('Nuclear NF-κB level (A.U.)')
    plt.xticks([0, 25, 50, 75, 100])
    plt.ylabel('Fold change in\nmRNA level - 1')
    plt.title('')
    plt.savefig('./Fig3F/TE_Km.pdf', bbox_inches='tight')
    plt.close()