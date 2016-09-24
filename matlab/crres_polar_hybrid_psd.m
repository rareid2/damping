function f = crres_polar_hybrid_psd(vperp, vpar, n, An, L, L_pp)
% A hybrid psd model that smoothly varies between CRRES fluxes outside the
% plasmasphere and POLAR fluxes inside the plasmasphere

% By Daniel Golden (dgolden1 at stanford dot edu) May 2010
% $Id: crres_polar_hybrid_psd.m 1205 2011-03-04 23:15:40Z dgolden $

% Using this scheme, the two weights are equal when L_pp == L, and when
% |L_pp - L| == 0.5, the dominant psd is weighted 12x more than the less
% dominant psd.  If the L value is huge (at high latitudes), the exp
% function overflows.
if L_pp - L > 1 % way inside plasmasphere
	f = suprathermal(vperp,vpar);
elseif L - L_pp > 1 % way outside plasmasphere
	f = crres_psd(n, An, vperp, vpar);
else
	f_polar = suprathermal(vperp,vpar);
	f_crres = crres_psd(n, An, vperp, vpar);

    % disp('f_polar: '); disp(f_polar);
    % disp('f_crres: '); disp(f_crres);
    
	w_polar = exp(5*(L_pp - L))./(1 + exp(5*(L_pp - L))); % higher weight inside plasmasphere
	w_crres = exp(5*(L - L_pp))./(1 + exp(5*(L - L_pp))); % higher weight outside plasmasphere

	% Find the centroid in log-space
	f = exp(centroid(log([f_polar, f_crres]), [w_polar, w_crres]));
end

% fprintf('L - L_pp = %0.1f, vperp = %0.1e, vpar = %+0.1e, fpolar = %0.1e, fcrres = %0.1e, f = %0.1e\n', ...
%     L - L_pp, vperp, vpar, f_polar, f_crres, f);

disp('')
