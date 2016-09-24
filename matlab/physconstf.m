function [R_E, c, EPS0, MU0, q, me, mp] = physconstf
% [R_E, c, EPS0, MU0, q, me] = physconstf
% Physical constants

% Originally by F. R. Foust
% Modified by Daniel Golden (dgolden1 at stanford dot edu) Feburary 2010
% $Id: physconstf.m 863 2010-04-01 23:56:05Z dgolden $


EPS0 = 8.854187817e-12;
MU0 = pi* 4e-7;
c = sqrt(1/EPS0/MU0);
R_E = 6370e3; % radius of earth in m
q = 1.60217646e-19; % elementary charge
me = 9.10938188e-31; % electron mass
mp = 1.67262158e-27; % proton mass
