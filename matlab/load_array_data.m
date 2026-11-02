function D = load_array_data(folder, varargin)
%LOAD_ARRAY_DATA Load the S-matrix and embedded element patterns of the 16 x 16 stacked-patch array.
%   D = LOAD_ARRAY_DATA() reads the dataset from the folder above this file.
%   D = LOAD_ARRAY_DATA(folder) reads folder/data/s_matrix.mat, folder/data/embedded_patterns.mat
%   and folder/data/element.json.
%   D = LOAD_ARRAY_DATA(folder, 'patterns', false) skips the embedded patterns (44 MB).
%
%   Fields of D (q: frequency, m and n: ports, d: direction):
%     f_norm      nf x 1          frequency divided by the centre frequency f0
%     S           nf x N x N      complex S-matrix, S(q, m, n) = b_m / a_n (port n excited, port m received)
%     ix, iy      N x 1           element column (along x) and row (along y); element k = (ix - 1)*16 + iy
%     x, y        N x 1           element position relative to the array centre, in wavelengths at f0
%     sigma_max   nf x 1          largest singular value of S at each frequency
%     theta_deg   Nd x 1          pattern direction: theta from +z (broadside)
%     phi_deg     Nd x 1          pattern direction: phi from +x toward +y
%     Gtheta      nf x N x Nd     embedded element pattern, theta component (complex single)
%     Gphi        nf x N x Nd     embedded element pattern, phi component (complex single)
%     element     struct          element geometry and materials (lengths in wavelengths at f0)
%
%   |Gtheta|^2 + |Gphi|^2 is the realized gain of one port for unit incident power, with all other
%   ports terminated in 50 ohm. All patterns use the array centre as phase reference and the time
%   convention exp(+j*omega*t). See README.md for the conventions and the accuracy.
%
%   Example:
%     D = load_array_data();
%     S0 = squeeze(D.S(D.f_norm == 1, :, :));   % S-matrix at the centre frequency
%
%   See also COPOL_PATTERN, ARRAY_GAIN, EXAMPLE_ARRAY_DATA.

if nargin < 1 || isempty(folder)
    folder = fullfile(fileparts(mfilename('fullpath')), '..');
end
p = inputParser;
p.addParameter('patterns', true, @(v) islogical(v) || isnumeric(v));
p.parse(varargin{:});

ddir = fullfile(folder, 'data');
s = load(fullfile(ddir, 's_matrix.mat'));
D = struct();
D.f_norm = s.f_norm(:);
D.S = s.S;
D.ix = s.ix(:);
D.iy = s.iy(:);
D.x = s.x_lambda(:);
D.y = s.y_lambda(:);
D.sigma_max = s.sigma_max(:);

if p.Results.patterns
    e = load(fullfile(ddir, 'embedded_patterns.mat'));
    assert(isequal(size(e.f_norm(:)), size(D.f_norm)) && max(abs(e.f_norm(:) - D.f_norm)) < 1e-12, ...
        'load_array_data: the pattern and S-matrix frequencies differ');
    D.theta_deg = e.theta_deg(:);
    D.phi_deg = e.phi_deg(:);
    D.Gtheta = e.Gtheta;
    D.Gphi = e.Gphi;
end

D.element = jsondecode(fileread(fullfile(ddir, 'element.json')));
end
