# Full-wave S-matrix and embedded element patterns of a 16 × 16 stacked-patch array

This dataset gives the full-wave scattering matrix and the embedded element patterns of a 16 × 16 phased array of
probe-fed stacked patches over an infinite ground plane. The paper that uses these data will be referenced here when it is published.

All frequencies are normalized to the centre frequency f0. All lengths are in free-space wavelengths at f0 (λ0).
The data can therefore be scaled to any band in which the material values hold.

## Files
| File | Content |
|---|---|
| `data/s_matrix.mat` | `S` (2 × 256 × 256, complex): scattering matrix at the two frequencies in `f_norm`; `f_norm` = f/f0 = [0.9685; 1] (lower edge of a 6.3 % fractional band, and the centre); `ix`, `iy` (256 × 1): element column and row; `x_lambda`, `y_lambda`: element positions relative to the array centre (λ0); `sigma_max`: largest singular value per frequency |
| `data/embedded_patterns.mat` | `Gtheta`, `Gphi` (2 × 256 × 5941, complex single): embedded element patterns; `theta_deg`, `phi_deg` (5941 × 1): directions; `f_norm` |
| `data/element.json` | element geometry and materials, lengths in λ0 |

Both MAT files are MATLAB version 7.3 (HDF5).

## Conventions
- **Elements and ports.** Element k has port k, with k = (ix − 1)·16 + iy. The index ix runs along x, iy along y. The pitch is 0.5 λ0 in both directions.
- **Polarization.** The E-field of the element is along x. The rows of the array (x) run along the E-plane, the columns (y) along the H-plane.
- **Ports.** 50 Ω coaxial lumped ports at the probe base in the ground plane. All ports are terminated in every solve.
- **S-matrix.** The S-matrix is used as computed (no passivation, no renormalization).
  - Reciprocity ‖S − Sᵀ‖/‖S‖ is below 5 × 10⁻⁶ at both frequencies.
  - The largest singular value is 0.9993 at f/f0 = 1 and 1.0112 at f/f0 = 0.9685; the small excess at the band edge comes from the truncation of the computational domain.
- **Patterns.** `Gtheta` and `Gphi` are the θ and φ components of the far field of one excited port, with all other ports terminated.
  - |Gtheta|² + |Gphi|² is the realized gain of that port.
  - θ is measured from +z (broadside); φ from +x toward +y.
  - The phase reference is the array centre.
  - The grid is 2° × 2° for θ ≤ 66° (θ = 0 once).
  - The co-polar component (Ludwig 3, x-polarized) is Gtheta·cos φ − Gphi·sin φ.

## How the data were computed
- **Solver:** finite-element solver Palace, version 0.16.0, with second-order Nédélec elements and GMRES (tolerance 10⁻⁶).
- **Domain:** the laminate and spacer end 0.25 λ0 beyond the outer unit cells. The air box extends 1 λ0 beyond the laminate and above the patches, with second-order absorbing walls.
- **Symmetry:** the array is solved as two half models with a mirror wall at y = 0 (between iy = 8 and iy = 9): a magnetic wall gives the even solution and an electric wall the odd solution. The full matrix is assembled from the two. On a 2 × 2 array this agrees with a full solve to 0.6 % (relative norm).
- **Patterns:** the far field is computed on a 2° × 2° grid from the outer absorbing surface. The surface is closed with images of the ground plane and of the symmetry wall.

## Accuracy
- **Power balance:** the power that radiates from the patterns is 1.04 to 1.06 times the accepted power (median over the ports; range 1.00 to 1.14). The pattern gains are therefore high by about 0.2 dB (median) and at most about 0.6 dB.
- **Nearest-neighbour couplings:** accurate to about ±1 dB (air-box convergence) to ±1.5 dB (comparison with Ansys HFSS, with a perfectly matched layer, on 2 × 2 and 4 × 4 arrays).
- **Weak couplings:** couplings below about −25 dB are indicative only. This includes the second E-plane neighbour.
- **Patterns:** on a 2 × 2 array, the embedded patterns agree with Ansys HFSS to 0.15 dB RMS (θ ≤ 66°).
- **Material values:** these are model values, not data-sheet values of a specific laminate or foam.

## Reading the data
MATLAB:
```matlab
s = load('data/s_matrix.mat');          % s.S(q, m, n): port n excited, port m received, frequency q
p = load('data/embedded_patterns.mat'); % p.Gtheta(q, k, d)
```
Python (h5py; MATLAB version 7.3 stores complex numbers as compound `real`/`imag` and reverses the axis order):
```python
import h5py, numpy as np
with h5py.File('data/s_matrix.mat', 'r') as f:
    S = f['S'][()]; S = (S['real'] + 1j * S['imag']).transpose(2, 1, 0)   # S[q, m, n]
    f_norm = f['f_norm'][()].ravel()
```

## Citation
Cite the dataset with its DOI (see the repository's release page or `CITATION.cff`) and the paper that uses it.

## Licence
See `LICENSE`.
