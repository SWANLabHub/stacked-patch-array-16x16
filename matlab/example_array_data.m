%EXAMPLE_ARRAY_DATA Load the dataset, print basic checks and plot two pattern cuts.
%   Run from any folder: run('<dataset>/matlab/example_array_data.m').
%   The printed values are the same as those of python/example_array_data.py.

here = fileparts(mfilename('fullpath'));
addpath(here);
D = load_array_data(fullfile(here, '..'));
[~, q0] = min(abs(D.f_norm - 1));                 % centre frequency
nf = numel(D.f_norm); N = numel(D.ix);
fprintf('%d ports, %d frequencies (f/f0 = %s), %d pattern directions\n', N, nf, mat2str(D.f_norm.', 6), numel(D.theta_deg));

% --- S-matrix checks ---------------------------------------------------------------
for q = 1:nf
    Sq = squeeze(D.S(q, :, :));
    fprintf('f/f0 = %.4f: reciprocity ||S - S.''||/||S|| = %.1e, largest singular value = %.4f\n', ...
        D.f_norm(q), norm(Sq - Sq.', 'fro') / norm(Sq, 'fro'), max(svd(Sq)));
end

% median coupling |S_mn| (dB) over all element pairs with a given lattice offset, centre frequency
S0 = squeeze(D.S(q0, :, :));
offsets = [1 0; 0 1; 1 1; 2 0; 0 2];
names = {'E-plane neighbour (1, 0)', 'H-plane neighbour (0, 1)', 'diagonal (1, 1)', ...
         'second E-plane (2, 0)', 'second H-plane (0, 2)'};
dix = D.ix - D.ix.'; diy = D.iy - D.iy.';
for r = 1:size(offsets, 1)
    sel = abs(dix) == offsets(r, 1) & abs(diy) == offsets(r, 2);
    fprintf('  %-26s median %6.1f dB\n', names{r}, median(20*log10(abs(S0(sel)))));
end

% active reflection coefficient for uniform broadside excitation
a = ones(N, 1) / sqrt(N);
Gact = (S0 * a) ./ a;
fprintf('active reflection, uniform broadside excitation: median %.1f dB, max %.1f dB\n', ...
    median(20*log10(abs(Gact))), max(20*log10(abs(Gact))));

% --- array patterns ---------------------------------------------------------------
G = array_gain(D, q0, ones(N, 1));
fprintf('uniform excitation, realized gain at broadside: %.2f dBi\n', 10*log10(G(D.theta_deg == 0)));

th0 = 30; ph0 = 0;                                    % steering direction
u0 = sind(th0)*cosd(ph0); v0 = sind(th0)*sind(ph0);
a = exp(-2j*pi*D.f_norm(q0)*(D.x*u0 + D.y*v0));       % time convention exp(+j*omega*t)
[G, Gco] = array_gain(D, q0, a);
[~, k] = max(G);
fprintf('steered to theta = %d, phi = %d deg: peak %.2f dBi at theta = %d, phi = %d deg (co-polar %.2f dBi)\n', ...
    th0, ph0, 10*log10(G(k)), D.theta_deg(k), D.phi_deg(k), 10*log10(Gco(k)));

% --- plots: co-polar cuts of a central element and the steered beam --------------------------
kc = (8 - 1)*16 + 8;                                  % element ix = 8, iy = 8
[Eco, ~] = copol_pattern(D.Gtheta(:, kc, :), D.Gphi(:, kc, :), D.phi_deg);
cut = @(phc) cut_index(D, phc);
figure('Name', 'Embedded element pattern, element (8, 8)');
tiledlayout(1, 2);
for pc = [0 90]
    nexttile; hold on; grid on;
    [i, t] = cut(pc);
    for q = 1:nf
        plot(t, 10*log10(abs(squeeze(Eco(q, 1, i))).^2), 'DisplayName', sprintf('f/f0 = %.4f', D.f_norm(q)));
    end
    xlabel('\theta (deg)'); ylabel('Co-polar realized gain (dBi)');
    title(sprintf('\\phi = %d deg', pc)); legend('Location', 'south');
end
figure('Name', 'Steered array pattern');
[i, t] = cut(0);
plot(t, 10*log10(Gco(i))); grid on;
xlabel('\theta (deg)'); ylabel('Co-polar realized gain (dBi)');
title(sprintf('Uniform excitation steered to \\theta = %d deg, \\phi = 0 (f/f0 = 1)', th0));

function [i, t] = cut_index(D, phc)
% directions of the cut phi = phc (theta >= 0) and phi = phc + 180 (shown as theta < 0)
pos = find(abs(D.phi_deg - phc) < 1e-6 | D.theta_deg == 0);
neg = find(abs(D.phi_deg - mod(phc + 180, 360)) < 1e-6 & D.theta_deg > 0);
t = [-D.theta_deg(neg); D.theta_deg(pos)];
i = [neg; pos];
[t, o] = sort(t); i = i(o);
end
