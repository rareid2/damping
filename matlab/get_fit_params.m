function [n, An] = get_fit_params(L, MLT, AE_level, b_debug, h_debug)
% [n, An] = get_fit_params(L, MLT, AE_level)
% 
% Determine fit params, e.g., from Bortnik [2007] Sec 2.2 for CRRES
% suprathermal fluxes

% L in R_E
% MLT in hours (0-24)
% AE_level either 1 (least active), 2, or 3 (most active)
% 
% J_perp units: cm^-2s^-1str^-1keV^-1
% 
% set b_debug to true for debugging info and plots
% give axis handles h_debug to plot debugging plots in a given axis handle

% By Daniel Golden (dgolden1 at stanford dot edu) May 2010
% $Id: get_fit_params.m 1301 2011-05-10 23:59:48Z dgolden $

%% Setup
% Raytracer root: Defining this locally because passing it around is annoying.

input_dir = 'crres_data';
% input_dir = fullfile(raytracer_root,'matlab','damping','crres_suprathermal_fluxes');

if ~exist('AE_level', 'var') || isempty(AE_level)
	AE_level = 3;
end
if ~isnumeric(AE_level) || ~isscalar(AE_level) || AE_level < 1 || AE_level > 3
	error('AE_level must be between 1 (least active) and 3 (most active)');
end
if ~exist('b_debug', 'var') || isempty(b_debug)
	b_debug = false;
end
if b_debug && (~exist('h_debug', 'var') || isempty(h_debug))
	figure;
	h_debug(1) = axes;
	figure;
	h_debug(2) = axes;
end

me = 9.10938188e-31; % electron mass (kg)

%% Load data
% Only load the file if it's not already in memory
persistent crres last_L last_MLT last_AE_level
if isempty(crres) || (last_AE_level ~= AE_level)
	d = dir(fullfile(input_dir, sprintf('crres*ae%d.mat', AE_level)));
	for kk = 1:length(d)
		[crres(kk).L, crres(kk).MLT, crres(kk).J_int, crres(kk).t_b] = ...
			get_fluxes(fullfile(input_dir, d(kk).name), 'gridfit');
		crres(kk).E = str2double(d(kk).name(8:12))/1e3; % keV
	end
	
	last_AE_level = AE_level;
end

%% Find this bin
for kk = 1:length(crres)
	row = nearest(L, crres(kk).L(:,1));
	col = nearest(MLT, crres(kk).MLT(1,:));
	J_perp(kk) = crres(kk).J_int(row, col);
	E(kk) = crres(kk).E;
end


% disp('Log E:'); disp(log(E));
% disp('Log Jp:');disp(log(J_perp));

%% Fit to flux
% Flux is of the form J_perp = J0/(E^m)
p = polyfit(log(E), log(J_perp), 1);
% disp('polyfit coeffs:'); disp(p);

J0 = exp(p(2));
m = -p(1);
% disp('J0: '); disp(J0); disp('m: '); disp(m);
%% Plot flux vs. E
if b_debug
	axes(h_debug(1)); cla;
	
	h = loglog(E, J_perp, 'kd');
	set(h, 'markersize', 10, 'color', 'k', 'markerfacecolor', 'k');
	hold on
	E_fit = logspace(log10(E(1)), log10(E(end)));
	% h = loglog(E_fit, exp(p(2))./E_fit.^(-p(1)));
	h = loglog(E_fit, J0./(E_fit.^m));
	set(h, 'linewidth', 2, 'color', 'k');

	grid on;
	xlabel('Energy (keV)');
	ylabel('Flux J_{perp} (cm^{-2}s^{-1}str^{-1}keV^{-1})');
    
    xlim([1e-1 2e1])
	increase_font;
end

%% Plot phase space density (PSD)
n = 2*m + 2;
An = 2*J0/(0.5*6.25e11*me).^(m-1);

if b_debug
	f = J_perp*6.25e11*me./(2*E/(6.25e11*me));
	v = sqrt(2*E/(6.25e11*me));

	axes(h_debug(2)); cla;
	
	h = loglog(v, f, 's');
	set(h, 'markersize', 10, 'color', 'k', 'markerfacecolor', 'k');

	% Plot PSD fit
	v_fit = logspace(log10(v(1)), log10(v(end)));
	hold on;
	h = loglog(v_fit, An./(v_fit.^n), '-');
	set(h, 'linewidth', 2, 'color', 'k');
	grid on
	xlabel('Velocity (cm/s)');
	ylabel('PSD f(v) (s^3/cm^6)');
    
    xlim([8e8 1e10])
	increase_font;
end
