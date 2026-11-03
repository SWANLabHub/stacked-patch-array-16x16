# Full-wave S-matrix and embedded element patterns of a 16 × 16 stacked-patch array

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23912674.svg)](https://doi.org/10.5281/zenodo.23912674)

This dataset gives the full-wave scattering matrix and the embedded element patterns of a 16 × 16 phased array of
probe-fed stacked patches over an infinite ground plane. It includes MATLAB and Python code that reads the data and
checks it. The paper that uses these data will be referenced here when it is published.

All frequencies are normalized to the centre frequency f0. All lengths are in free-space wavelengths at f0 (λ0).
The data can therefore be scaled to any band in which the material values hold.

## Files
| File | Content |
|---|---|
| `data/s_matrix.mat` | `S` (2 × 256 × 256, complex): scattering matrix at the two frequencies in `f_norm`; `f_norm` = f/f0 = [0.9685; 1] (lower edge of a 6.3 % fractional band, and the centre); `ix`, `iy` (256 × 1): element column and row; `x_lambda`, `y_lambda`: element positions relative to the array centre (λ0); `sigma_max`: largest singular value per frequency |
| `data/embedded_patterns.mat` | `Gtheta`, `Gphi` (2 × 256 × 5941, complex single): embedded element patterns; `theta_deg`, `phi_deg` (5941 × 1): directions; `f_norm` |
| `data/element.json` | element geometry and materials, lengths in λ0 |
| `matlab/load_array_data.m` | reads all data into one structure |
| `matlab/copol_pattern.m`, `matlab/array_gain.m` | co-polar and cross-polar components; realized gain of the array for given port excitations |
| `matlab/example_array_data.m` | prints the checks below and plots two pattern cuts |
| `python/array_data.py` | the same three functions in Python (`load`, `copol_pattern`, `array_gain`) |
| `python/example_array_data.py` | the same checks and plots in Python |
| `python/requirements.txt` | numpy, h5py; matplotlib for the plots |
| `SHA256SUMS` | checksums of the data files |

Both MAT files are MATLAB version 7.3 files, which are HDF5 files. Any HDF5 reader can open them. MATLAB stores a
complex array as a compound type with the fields `real` and `imag`, and an HDF5 reader sees the axes in reverse
order. `python/array_data.py` undoes both.

## Quick start
MATLAB (R2019b or newer):
```matlab
addpath('matlab');
D = load_array_data();                     % D.S(q, m, n), D.Gtheta(q, k, d), ...
run('matlab/example_array_data.m');        % checks and plots
```
Python 3 with numpy and h5py:
```python
import sys; sys.path.insert(0, 'python')
from array_data import load, array_gain
D = load()                                 # D.S[q, m, n], D.Gtheta[q, k, d], ... (zero-based indices)
```
or run `python python/example_array_data.py` (add `--no-plot` or `--save <folder>`).

## Conventions
- **Elements and ports.** Element k has port k, with k = (ix − 1)·16 + iy (k = 1 … 256). The index ix runs along x, iy along y. The pitch is 0.5 λ0 in both directions. In Python, element k has the index k − 1.
- **Polarization.** The E-field of the element is along x. The rows of the array (x) run along the E-plane, the columns (y) along the H-plane.
- **Ports.** 50 Ω coaxial lumped ports at the probe base in the ground plane. All ports are terminated in every solve.
- **Time convention.** exp(+jωt). The far field of an element at position (x, y) carries the phase factor exp(+jk0(x·u + y·v)), with u = sin θ cos φ and v = sin θ sin φ. The examples check this: the steered beam peaks at the commanded direction.
- **S-matrix.** S(q, m, n) = b_m / a_n at frequency q: port n is excited and port m receives. The S-matrix is used as computed (no passivation, no renormalization).
  - Reciprocity ‖S − Sᵀ‖/‖S‖ is below 5 × 10⁻⁶ at both frequencies.
  - The largest singular value is 0.9993 at f/f0 = 1 and 1.0112 at f/f0 = 0.9685. The small excess at the band edge comes from the truncation of the computational domain.
- **Patterns.** `Gtheta` and `Gphi` are the θ and φ components of the far field of one excited port, with all other ports terminated in 50 Ω.
  - |Gtheta|² + |Gphi|² is the realized gain of that port for unit incident power.
  - θ is measured from +z (broadside); φ from +x toward +y.
  - The grid is 2° × 2° for θ ≤ 66° (θ = 0 once).
  - The co-polar component (Ludwig 3, x-polarized) is Gtheta·cos φ − Gphi·sin φ; the cross-polar component is Gtheta·sin φ + Gphi·cos φ.
- **Array patterns.** All element patterns use the array centre as phase reference. For incident port waves a_n with Σ|a_n|² = 1, the realized gain of the array is |Σ a_n Gtheta_n|² + |Σ a_n Gphi_n|² (`array_gain`). To steer the beam to (u0, v0) at frequency f, use a_n ∝ exp(−j2π(f/f0)(x_n·u0 + y_n·v0)) with x_n, y_n in λ0.

## Checks
`example_array_data.m` and `example_array_data.py` print the following values. Use them to check a download or
a new reader.

| Quantity | Value |
|---|---|
| reciprocity ‖S − Sᵀ‖/‖S‖ (f/f0 = 0.9685 / 1) | 2.9 × 10⁻⁶ / 4.7 × 10⁻⁶ |
| largest singular value (f/f0 = 0.9685 / 1) | 1.0112 / 0.9993 |
| median coupling at f/f0 = 1: E-plane neighbour (1, 0) / H-plane neighbour (0, 1) / diagonal (1, 1) | −15.4 / −17.3 / −22.0 dB |
| median coupling at f/f0 = 1: second E-plane (2, 0) / second H-plane (0, 2) | −32.2 / −27.7 dB |
| active reflection, uniform broadside excitation, f/f0 = 1: median / maximum | −16.0 / −12.2 dB |
| realized gain, uniform excitation, broadside, f/f0 = 1 | 29.31 dBi |
| realized gain, uniform excitation steered to θ = 30°, φ = 0, f/f0 = 1 | 29.00 dBi, peak at θ = 30°, φ = 0 |

## How the data were computed
- **Solver:** finite-element solver Palace, version 0.16.0, with second-order Nédélec elements and GMRES (tolerance 10⁻⁶).
- **Domain:** the laminate and spacer end 0.25 λ0 beyond the outer unit cells. The air box extends 1 λ0 beyond the laminate and above the patches, with second-order absorbing walls.
- **Symmetry:** the array is solved as two half models with a mirror wall at y = 0 (between iy = 8 and iy = 9): a magnetic wall gives the even solution and an electric wall the odd solution. The full matrix is assembled from the two. On a 2 × 2 array this agrees with a full solve to 0.6 % (relative norm).
- **Patterns:** the far field is computed on a 2° × 2° grid from the outer absorbing surface. The surface is closed with images of the ground plane and of the symmetry wall.

## Accuracy
- **Power balance:** the power that radiates from the patterns is 1.04 to 1.06 times the accepted power (median over the ports; range 1.00 to 1.14). The pattern gains are therefore high by about 0.2 dB (median) and at most about 0.6 dB.
- **Nearest-neighbour couplings:** accurate to about ±1 dB (air-box convergence) to ±1.5 dB (comparison with Ansys HFSS, with a perfectly matched layer, on 2 × 2 and 4 × 4 arrays).
- **Weak couplings:** couplings below about −25 dB are indicative only. This includes the second neighbours in both planes.
- **Patterns:** on a 2 × 2 array, the embedded patterns agree with Ansys HFSS to 0.15 dB RMS (θ ≤ 66°).
- **Material values:** these are model values, not data-sheet values of a specific laminate or foam.

## Citation
Cite the dataset with its DOI and the paper that uses it.
- All versions (resolves to the latest version): https://doi.org/10.5281/zenodo.23912674
- Version 1.0.0: https://doi.org/10.5281/zenodo.23912675

Kalfa M. Full-wave S-matrix and embedded element patterns of a 16 × 16 stacked-patch array, version 1.0.0 [dataset].
Zenodo, 2026. https://doi.org/10.5281/zenodo.23912675

## Licence
The data and this documentation: Creative Commons Attribution 4.0 International (CC BY 4.0), see `LICENSE`.
The code in `matlab/` and `python/`: MIT licence, see `LICENSE-CODE`.
