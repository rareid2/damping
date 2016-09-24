function f = crres_psd(n, An, vperp, vpar)
% Get phase space density using CRRES suprathermal fluxes from Bortnik
% [2007]
% 
% vperp and vpar in m/s

% By Daniel Golden (dgolden1 at stanford dot edu) May 2010
% $Id: crres_psd.m 1205 2011-03-04 23:15:40Z dgolden $

% Just a crutch to avoid the singularity
v0=1;

% Convert to cm/s
v = 100*sqrt(vperp.^2 + vpar.^2 + v0);

f = An./(v.^n);

% Convert to s^3/m^6 from s^3/cm^6 
f = f*100^6;