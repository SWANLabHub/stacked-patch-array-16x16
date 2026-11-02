"""Read the S-matrix and embedded element patterns of the 16 x 16 stacked-patch array.

The data files are MATLAB version 7.3 MAT files, which are HDF5 files. MATLAB stores complex
arrays as a compound type with the fields 'real' and 'imag' and writes the axes in reverse order;
this module undoes both, so the arrays have the same shape and index order as in MATLAB.
Element and port indices are zero-based here: element k of the README is index k - 1.

Fields of the object that load() returns (q: frequency, m and n: ports, d: direction):
    f_norm     (nf,)          frequency divided by the centre frequency f0
    S          (nf, N, N)     complex S-matrix, S[q, m, n] = b_m / a_n (port n excited, port m received)
    ix, iy     (N,)           element column (along x) and row (along y), 1 to 16
    x, y       (N,)           element position relative to the array centre, in wavelengths at f0
    sigma_max  (nf,)          largest singular value of S at each frequency
    theta_deg  (Nd,)          pattern direction: theta from +z (broadside)
    phi_deg    (Nd,)          pattern direction: phi from +x toward +y
    Gtheta     (nf, N, Nd)    embedded element pattern, theta component (complex64)
    Gphi       (nf, N, Nd)    embedded element pattern, phi component (complex64)
    element    dict           element geometry and materials (lengths in wavelengths at f0)

|Gtheta|^2 + |Gphi|^2 is the realized gain of one port for unit incident power, with all other
ports terminated in 50 ohm. All patterns use the array centre as phase reference and the time
convention exp(+j*omega*t). See README.md for the conventions and the accuracy.

Requires numpy and h5py.
"""
import json
import os
from types import SimpleNamespace

import h5py
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))


def _read(f, name):
    """Read one MATLAB v7.3 variable in MATLAB axis order, complex if MATLAB stored it complex."""
    d = f[name][()]
    if d.dtype.names and 'real' in d.dtype.names and 'imag' in d.dtype.names:
        d = d['real'] + 1j * d['imag']
    return np.asarray(d).transpose()


def load(folder=None, patterns=True):
    """Load the dataset from folder (default: the folder above this file).

    With patterns=False, the embedded patterns (44 MB) are not read.
    """
    folder = folder or os.path.join(_HERE, '..')
    ddir = os.path.join(folder, 'data')
    D = SimpleNamespace()
    with h5py.File(os.path.join(ddir, 's_matrix.mat'), 'r') as f:
        D.f_norm = _read(f, 'f_norm').ravel()
        D.S = _read(f, 'S').astype(np.complex128)
        D.ix = _read(f, 'ix').ravel().astype(int)
        D.iy = _read(f, 'iy').ravel().astype(int)
        D.x = _read(f, 'x_lambda').ravel()
        D.y = _read(f, 'y_lambda').ravel()
        D.sigma_max = _read(f, 'sigma_max').ravel()
    if patterns:
        with h5py.File(os.path.join(ddir, 'embedded_patterns.mat'), 'r') as f:
            fn = _read(f, 'f_norm').ravel()
            if fn.shape != D.f_norm.shape or np.max(np.abs(fn - D.f_norm)) > 1e-12:
                raise ValueError('the pattern and S-matrix frequencies differ')
            D.theta_deg = _read(f, 'theta_deg').ravel()
            D.phi_deg = _read(f, 'phi_deg').ravel()
            D.Gtheta = _read(f, 'Gtheta').astype(np.complex64)
            D.Gphi = _read(f, 'Gphi').astype(np.complex64)
    with open(os.path.join(ddir, 'element.json'), encoding='utf-8') as fh:
        D.element = json.load(fh)
    return D


def copol_pattern(Gtheta, Gphi, phi_deg):
    """Co-polar and cross-polar components (Ludwig 3, x-polarized reference).

    Eco = Gtheta cos(phi) - Gphi sin(phi), Ecr = Gtheta sin(phi) + Gphi cos(phi).
    The direction index must be the last axis of Gtheta and Gphi.
    """
    c = np.cos(np.deg2rad(phi_deg))
    s = np.sin(np.deg2rad(phi_deg))
    return Gtheta * c - Gphi * s, Gtheta * s + Gphi * c


def array_gain(D, q, a):
    """Realized gain (linear) of the array in every pattern direction, frequency index q.

    a: incident port waves (N,), normalized here to unit total incident power. All embedded
    patterns share the array centre as phase reference, so the array field is the weighted sum
    of the element fields. Returns (G, Gco, Gcr): total, co-polar and cross-polar gain.
    """
    a = np.asarray(a, dtype=complex).ravel()
    a = a / np.linalg.norm(a)
    Et = a @ D.Gtheta[q].astype(complex)
    Ep = a @ D.Gphi[q].astype(complex)
    Eco, Ecr = copol_pattern(Et, Ep, D.phi_deg)
    return np.abs(Et) ** 2 + np.abs(Ep) ** 2, np.abs(Eco) ** 2, np.abs(Ecr) ** 2
