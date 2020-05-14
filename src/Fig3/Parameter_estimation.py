"""Related to Fig.3E"""
import os
import sys
import numpy as np
from scipy.optimize import differential_evolution
import matplotlib.pyplot as plt
import seaborn as sns


dose = np.array([0, 23.52, 30.858, 44, 100])
NFkB = np.linspace(0, 100, 100001)


def objective_se(x):
    """Super Enhancer"""
    average_se = np.array(
        [0, 0.100886654, 0.179533988, 0.65220359, 1.722985887]
    )
    model = NFkB**x[0]/(NFkB**x[0]+x[1]*(x[2]**x[0]+NFkB**x[0]))*x[3]
    rss = 0
    for i, j in enumerate(dose):
        rss += (model[int(j * 1000)] - average_se[i])**2

    return rss


def objective_te(x):
    """Typical Enhancer"""
    average_te = np.array(
        [0, 0.211417158, 0.307286311, 0.448678085, 1.040479864]
    )
    model = NFkB**x[0]/(NFkB**x[0]+x[1]*(x[2]**x[0]+NFkB**x[0]))*x[3]
    rss = 0
    for i, j in enumerate(dose):
        rss += (model[int(j * 1000)] - average_te[i])**2

    return rss


def my_rc_params():
    plt.rcParams['font.size'] = 24
    plt.rcParams['axes.linewidth'] = 2


if __name__ == '__main__':
    if not (os.path.isfile('params_se.npy') and 
            os.path.isfile('params_te.npy')):
        bounds = [(1, 12), (1, 100), (1, 1000), (1, 100)]
        n = 50
        params_se = np.empty((n, 4))
        params_te = np.empty((n, 4))
        for i in range(n):
            sys.stdout.write('\r{:d}/{:d}'.format(i + 1, n))
            result_se = differential_evolution(
                objective_se, bounds, strategy='best1bin', polish=True
            )
            result_te = differential_evolution(
                objective_te, bounds, strategy='best1bin', polish=True
            )
            params_se[i, :] = result_se.x
            params_te[i, :] = result_te.x
        np.save('params_se.npy', params_se)
        np.save('params_te.npy', params_te)
    else:
        params_se = np.load('params_se.npy')
        params_te = np.load('params_te.npy')

    os.makedirs('./Fig3E', exist_ok=True)

    """N
    """
    my_rc_params()
    fig = plt.figure(figsize=(3, 6))

    ax = sns.boxplot(
        data=np.hstack(
            (params_se[:, 0][:, np.newaxis], params_te[:, 0][:, np.newaxis])
        ), orient='v', linewidth=1, palette='Set2'
    )
    sns.despine()
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['SE', 'TE'])
    ax.set_ylim(0.5, 4.5)
    ax.set_yticks([1, 2, 3, 4])
    ax.set_ylabel('N')
    plt.savefig('./Fig3E/range_N.pdf', bbox_inches='tight')
    plt.close()

    """Km
    """
    my_rc_params()
    fig = plt.figure(figsize=(3, 6))

    ax = sns.boxplot(
        data=np.hstack(
            (params_se[:, 2][:, np.newaxis], params_te[:, 2][:, np.newaxis])
        ), orient='v', linewidth=1, palette='Set2'
    )
    sns.despine()
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['SE', 'TE'])
    ax.set_yticks([50, 100, 150, 200, 250, 300])
    ax.set_ylabel('Km')
    plt.savefig('./Fig3E/range_Km.pdf', bbox_inches='tight')
    plt.close()

    """K1
    """
    my_rc_params()
    fig = plt.figure(figsize=(3, 6))

    ax = sns.boxplot(
        data=np.hstack(
            (params_se[:, 1][:, np.newaxis], params_te[:, 1][:, np.newaxis])
        ), orient='v', linewidth=1, palette='Set2'
    )
    sns.despine()
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['SE', 'TE'])
    ax.set_yticks([0, 10, 20, 30, 40, 50])
    ax.set_ylabel('K1')
    plt.savefig('./Fig3E/range_K1.pdf', bbox_inches='tight')
    plt.close()

    """K2
    """
    my_rc_params()
    fig = plt.figure(figsize=(3, 6))

    ax = sns.boxplot(
        data=np.hstack(
            (params_se[:, 3][:, np.newaxis], params_te[:, 3][:, np.newaxis])
        ), orient='v', linewidth=1, palette='Set2'
    )
    sns.despine()
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['SE', 'TE'])
    ax.set_ylabel('K2')
    plt.savefig('./Fig3E/range_K2.pdf', bbox_inches='tight')
    plt.close()