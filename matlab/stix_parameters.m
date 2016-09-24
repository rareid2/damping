function [S,D,P,R,L] = stix_parameters(w, qs, Ns, ms, nus, B0mag)
% Compute the stix parameters for a multicomponent plasma
% w = frequency
% qs = vector of charges
% Ns = vector of number densities in m^-3
% ms = vector of masses
% nus = vector of collision frequencies
% B0mag = the magnetic field magnitude

  [~, ~, EPS0, ~] = physconstf;

  wps2 = (Ns.*qs.^2./ms./EPS0).*(w./(w+j*nus));
  wcs = ((qs*B0mag)./ms).*(w./(w+j*nus));

  % Evaluate the stix parameters given the multicomponent plasma relations
  R = 1-sum(wps2./(w*(w+wcs)));
  L = 1-sum(wps2./(w*(w-wcs)));
  P = 1-sum(wps2./w^2);
  S = 1/2*(R+L);
  D = 1/2*(R-L);
  % fprintf('w: %0.2f qs: %e, Ns: %e, ms: %e, nus: %e, B0: %e\n',w,qs,Ns,ms,nus,B0mag);
  % fprintf('R: %0.2f L: %0.2f P: %0.2f S: %0.2f D: %0.2f\n',R,L,P,S,D);