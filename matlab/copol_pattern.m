function [Eco, Ecr] = copol_pattern(Gtheta, Gphi, phi_deg)
%COPOL_PATTERN Co-polar and cross-polar components (Ludwig 3, x-polarized reference).
%   [Eco, Ecr] = COPOL_PATTERN(Gtheta, Gphi, phi_deg) returns
%     Eco = Gtheta .* cos(phi) - Gphi .* sin(phi)
%     Ecr = Gtheta .* sin(phi) + Gphi .* cos(phi)
%   The direction index must be the last dimension of Gtheta and Gphi (as in LOAD_ARRAY_DATA),
%   or Gtheta and Gphi are vectors with one entry per direction. |Eco|^2 + |Ecr|^2 equals
%   |Gtheta|^2 + |Gphi|^2.
%
%   Example (co-polar pattern of element 136 at the centre frequency):
%     D = load_array_data();
%     q = find(D.f_norm == 1);
%     Eco = copol_pattern(squeeze(D.Gtheta(q, 136, :)), squeeze(D.Gphi(q, 136, :)), D.phi_deg);
%
%   See also LOAD_ARRAY_DATA, ARRAY_GAIN.

assert(isequal(size(Gtheta), size(Gphi)), 'copol_pattern: Gtheta and Gphi differ in size');
nd = numel(phi_deg);
if isvector(Gtheta)
    assert(numel(Gtheta) == nd, 'copol_pattern: one pattern value per direction expected');
    shp = size(Gtheta);
else
    assert(size(Gtheta, ndims(Gtheta)) == nd, 'copol_pattern: the direction index must be the last dimension');
    shp = [ones(1, ndims(Gtheta) - 1), nd];
end
c = reshape(cosd(double(phi_deg(:))), shp);
s = reshape(sind(double(phi_deg(:))), shp);
Eco = Gtheta .* c - Gphi .* s;
Ecr = Gtheta .* s + Gphi .* c;
end
