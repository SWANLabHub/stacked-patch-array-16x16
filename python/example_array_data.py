"""Load the dataset, print basic checks and plot two pattern cuts.

Usage: python example_array_data.py [--no-plot] [--save DIR]
The printed values are the same as those of matlab/example_array_data.m.
Plots need matplotlib; with --save DIR they are written as PNG files instead of shown.
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from array_data import array_gain, copol_pattern, load  # noqa: E402


def cut_index(D, phc):
    """Directions of the cut phi = phc (theta >= 0) and phi = phc + 180 (returned as theta < 0)."""
    pos = np.flatnonzero((np.abs(D.phi_deg - phc) < 1e-6) | (D.theta_deg == 0))
    neg = np.flatnonzero((np.abs(D.phi_deg - (phc + 180) % 360) < 1e-6) & (D.theta_deg > 0))
    t = np.concatenate([-D.theta_deg[neg], D.theta_deg[pos]])
    i = np.concatenate([neg, pos])
    o = np.argsort(t)
    return i[o], t[o]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--no-plot', action='store_true')
    ap.add_argument('--save', metavar='DIR')
    args = ap.parse_args()

    D = load()
    q0 = int(np.argmin(np.abs(D.f_norm - 1)))       # centre frequency
    nf, N = D.f_norm.size, D.ix.size
    print('%d ports, %d frequencies (f/f0 = [%s]), %d pattern directions'
          % (N, nf, ' '.join('%.6g' % v for v in D.f_norm), D.theta_deg.size))

    # --- S-matrix checks ------------------------------------------------------------
    for q in range(nf):
        Sq = D.S[q]
        print("f/f0 = %.4f: reciprocity ||S - S.'||/||S|| = %.1e, largest singular value = %.4f"
              % (D.f_norm[q], np.linalg.norm(Sq - Sq.T) / np.linalg.norm(Sq),
                 np.linalg.svd(Sq, compute_uv=False)[0]))

    S0 = D.S[q0]
    dix = np.abs(D.ix[:, None] - D.ix[None, :])
    diy = np.abs(D.iy[:, None] - D.iy[None, :])
    for (ox, oy), name in zip([(1, 0), (0, 1), (1, 1), (2, 0), (0, 2)],
                              ['E-plane neighbour (1, 0)', 'H-plane neighbour (0, 1)', 'diagonal (1, 1)',
                               'second E-plane (2, 0)', 'second H-plane (0, 2)']):
        sel = (dix == ox) & (diy == oy)
        print('  %-26s median %6.1f dB' % (name, np.median(20 * np.log10(np.abs(S0[sel])))))

    a = np.ones(N) / np.sqrt(N)
    gact = (S0 @ a) / a
    print('active reflection, uniform broadside excitation: median %.1f dB, max %.1f dB'
          % (np.median(20 * np.log10(np.abs(gact))), np.max(20 * np.log10(np.abs(gact)))))

    # --- array patterns ---------------------------------------------------------------
    G, _, _ = array_gain(D, q0, np.ones(N))
    print('uniform excitation, realized gain at broadside: %.2f dBi'
          % (10 * np.log10(G[D.theta_deg == 0][0])))

    th0, ph0 = 30, 0                                   # steering direction
    u0 = np.sin(np.deg2rad(th0)) * np.cos(np.deg2rad(ph0))
    v0 = np.sin(np.deg2rad(th0)) * np.sin(np.deg2rad(ph0))
    a = np.exp(-2j * np.pi * D.f_norm[q0] * (D.x * u0 + D.y * v0))   # time convention exp(+j*omega*t)
    G, Gco, _ = array_gain(D, q0, a)
    k = int(np.argmax(G))
    print('steered to theta = %d, phi = %d deg: peak %.2f dBi at theta = %d, phi = %d deg (co-polar %.2f dBi)'
          % (th0, ph0, 10 * np.log10(G[k]), D.theta_deg[k], D.phi_deg[k], 10 * np.log10(Gco[k])))

    if args.no_plot:
        return
    try:
        import matplotlib
        if args.save:
            matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print('matplotlib not found: no plots')
        return

    kc = (8 - 1) * 16 + 8 - 1                          # element ix = 8, iy = 8 (zero-based index)
    Eco, _ = copol_pattern(D.Gtheta[:, kc, :], D.Gphi[:, kc, :], D.phi_deg)
    fig1, axs = plt.subplots(1, 2, figsize=(10, 4))
    for ax, pc in zip(axs, [0, 90]):
        i, t = cut_index(D, pc)
        for q in range(nf):
            ax.plot(t, 10 * np.log10(np.abs(Eco[q, i]) ** 2), label='f/f0 = %.4f' % D.f_norm[q])
        ax.set(xlabel='theta (deg)', ylabel='Co-polar realized gain (dBi)', title='phi = %d deg' % pc)
        ax.grid(True)
        ax.legend(loc='lower center')
    fig1.suptitle('Embedded element pattern, element (8, 8)')
    fig2, ax = plt.subplots(figsize=(6, 4))
    i, t = cut_index(D, 0)
    ax.plot(t, 10 * np.log10(Gco[i]))
    ax.set(xlabel='theta (deg)', ylabel='Co-polar realized gain (dBi)',
           title='Uniform excitation steered to theta = %d deg, phi = 0 (f/f0 = 1)' % th0)
    ax.grid(True)
    if args.save:
        fig1.savefig(os.path.join(args.save, 'element_pattern_cuts.png'), dpi=120, bbox_inches='tight')
        fig2.savefig(os.path.join(args.save, 'steered_array_pattern.png'), dpi=120, bbox_inches='tight')
    else:
        plt.show()


if __name__ == '__main__':
    main()
