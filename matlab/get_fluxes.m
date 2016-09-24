function [L, MLT, J_perp_out, t_b] = get_fluxes(fluxfilename, interp_type)
% [L, MLT, J_perp_out, t_b] = get_fluxes(fluxfilename, interp_type)
% 
% Interpolate crres fluxes a few different ways
%
% interp_type can be one of:
%  'none' (default)
%  'linear' for interpolation using TriScatteredInterp, or
%  'gridfit' to use the FEX gridfit function (a "gradient regularizer")

% By Daniel Golden (dgolden1 at stanford dot edu) May 2010
% $Id: get_fluxes.m 1205 2011-03-04 23:15:40Z dgolden $

%% Setup
if ~exist('interp_type', 'var') || isempty(interp_type)
	interp_type = 'none';
end

%% Load data
load(fluxfilename);
J_perp(J_perp == 1e12 | t_b < 10) = nan; % 1e12 is actually nan
t_b(t_b == 1e12 | t_b < 10) = nan;

%% Interpolate
switch interp_type
	case 'none'
		J_perp_out = J_perp;
% 		t_b = t_b;
	case 'linear'		
		b_valid = ~isnan(J_perp);
		F_J = TriScatteredInterp(L(b_valid), MLT(b_valid), J_perp(b_valid));
		J_int = F_J(L, MLT);

		% Some interpolated values get NaNs (in the innermost or outermost
        % L) because there are no valid values around them.  Turn them into
        % nearest neighbor interpolants
        F_J_nearest = TriScatteredInterp(L(b_valid), MLT(b_valid), J_perp(b_valid), 'nearest');
        J_int_nearest = F_J_nearest(L, MLT);
        
        J_int(isnan(J_int)) = J_int_nearest(isnan(J_int));
		
		J_perp_out = J_int;
	case 'gridfit'
		% addpath(fullfile(danmatlabroot, 'gridfit'));

		% Extend MLT to go from -12 to +36 for the fit
		MLT_fit = [MLT(:, 13:24)-24, MLT, MLT(:, 1:12)+24];
		L_fit = [L(:, 13:24), L, L(:, 1:12)];
		J_perp_fit = [J_perp(:, 13:24), J_perp, J_perp(:, 1:12)];
		b_valid = ~isnan(J_perp_fit);
		
		% Smooth 'n' fit
		J_int = gridfit(MLT_fit(b_valid), L_fit(b_valid), log10(J_perp_fit(b_valid)), ...
			sort(unique(MLT_fit(:))), sort(unique(L_fit(:))), 'smoothness', 3);
		J_perp_out = 10.^(J_int(:, 13:36));
	otherwise
		error('Weird interp_type: %s', interp_type);
end
