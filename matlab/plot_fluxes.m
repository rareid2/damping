function plot_fluxes(filename, str_type, interp_type, ax)
% plot_fluxes(filename, str_type, interp_type, ax)
% 
% Plot suprathermal fluxes a la Bortnik [2007] (Modeling the propagation
% characteristics...)
% 
% str_type is one of 'fluxes' or 'time'
% 
% interp_type is one of 'none', 'linear', or 'gridfit'
% 
% Use axis ax if provided

% By Daniel Golden (dgolden1 at stanford dot edu) May 2010
% $Id: plot_fluxes.m 1205 2011-03-04 23:15:40Z dgolden $

%% Setup
if ~exist('ax', 'var') || isempty(ax)
	figure;
	ax = axes;
	b_own_axes = true;
else
	b_own_axes = false;
end
if ~exist('str_type', 'var') || isempty(str_type)
	str_type = 'fluxes';
end
if ~exist('interp_type', 'var') || isempty(interp_type)
	interp_type = 'none';
end


%% Get fluxes
[L, MLT, J_plot, t_b] = get_fluxes(filename, interp_type);

%% Massage so pcolor works better
J_plot = [J_plot, zeros(size(J_plot, 1), 1)];
L = [L, L(:,1)];
MLT = [MLT, MLT(:,1)];
t_b = [t_b, zeros(size(t_b, 1), 1)];

%% Plot
saxes(ax);

if strcmp(str_type, 'fluxes')
	p = pcolor(L.*cos((MLT-0.5)/24*2*pi - pi/2), L.*sin((MLT-0.5)/24*2*pi - pi/2), log10(J_plot));
	set(p, 'linestyle', 'none');
	caxis(log10([2e5 2e8]))
	axis square;
elseif strcmp(str_type, 'time')
	p = pcolor(L.*cos((MLT-0.5)/24*2*pi - pi/2), L.*sin((MLT-0.5)/24*2*pi - pi/2), log10(t_b));
	set(p, 'linestyle', 'none');
	caxis(log10([1 200]));
	axis square;
else
	error('Weird value for str_type: %s', str_type);
end

% Plot day/night terminator
hold on;
patch(cos(linspace(0, pi, 25)), sin(linspace(0, pi, 25)), 'w');
patch(cos(linspace(pi, 2*pi, 25)), sin(linspace(pi, 2*pi, 25)), 'k');

% Title
[~, fname] = fileparts(filename);
E = str2double(fname(8:12));
AE = str2double(fname(16));
if E < 1e3
% 	title(sprintf('E = %0.0f eV, AE = %d', E, AE));
	title(sprintf('%0.0f eV', E));
else
% 	title(sprintf('E = %0.2f keV, AE = %d', E/1e3, AE));
	title(sprintf('%0.2f keV', E/1e3));
end

% Ditch x and y labels
set(ax, 'xtick', [], 'ytick', []);

if b_own_axes
	increase_font;
end

%% Overlay grids
hold on;
theta = (0:24)/24*2*pi;
L_grid = 2:8;

% L grid
for kk = 1:length(L_grid)
	plot(L_grid(kk)*cos(theta), L_grid(kk)*sin(theta), 'k--');
end

% MLT grid
for kk = 1:length(theta(1:end-1))
	plot([1 8]*cos(theta(kk)), [1 8]*sin(theta(kk)), 'k--');
end

% Rectangular grid at 0, 6, 12, 18 MLT
plot([-8 8], [0 0], 'k-');
plot([0 0], [-8 8], 'k-');
