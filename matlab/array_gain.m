function [G, Gco, Gcr] = array_gain(D, q, a)
%ARRAY_GAIN Realized gain of the array for given incident port waves.
%   G = ARRAY_GAIN(D, q, a) returns the realized gain (linear) in every direction of
%   D.theta_deg, D.phi_deg at frequency index q, for incident port waves a (N x 1).
%   The waves are normalized to unit total incident power, sum(abs(a).^2) = 1.
%   Because all embedded patterns share the array centre as phase reference, the array field
%   is the weighted sum of the element fields:
%     G = |sum_n a_n Gtheta_n|^2 + |sum_n a_n Gphi_n|^2
%   [G, Gco, Gcr] = ARRAY_GAIN(...) also returns the co-polar and cross-polar parts (Ludwig 3).
%
%   Example (uniform excitation steered to theta = 30 deg, phi = 0, at the centre frequency):
%     D = load_array_data();
%     q = find(D.f_norm == 1);
%     u0 = sind(30); v0 = 0;
%     a = exp(-2j*pi*D.f_norm(q)*(D.x*u0 + D.y*v0));     % time convention exp(+j*omega*t)
%     G = array_gain(D, q, a);
%
%   See also LOAD_ARRAY_DATA, COPOL_PATTERN.

assert(isfield(D, 'Gtheta'), 'array_gain: load the dataset with the embedded patterns');
n = size(D.Gtheta, 2);
a = double(a(:));
assert(numel(a) == n, 'array_gain: one incident wave per port expected');
a = a / norm(a);
Et = (a.' * double(reshape(D.Gtheta(q, :, :), n, []))).';    % Nd x 1
Ep = (a.' * double(reshape(D.Gphi(q, :, :), n, []))).';
G = abs(Et).^2 + abs(Ep).^2;
if nargout > 1
    [Eco, Ecr] = copol_pattern(Et, Ep, D.phi_deg);
    Gco = abs(Eco).^2;
    Gcr = abs(Ecr).^2;
end
end
