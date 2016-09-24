function pp = kp_to_pp(kp)
% pp = kp_to_pp(kp)
% 
% Function to convert from kp to PP_L of GCPM model
% I had to actually look with my eyeballs to calculate these values, but
% this function interpolates

% By Daniel Golden (dgolden1 at stanford dot edu) May 2010
% $Id: kp_to_pp.m 969 2010-05-26 22:04:27Z dgolden $

kp_map = [30  33  37  40  43  47  50  53  57  60  63  67  70  73  77  80]/10;
pp_map = [4.4 4.3 4.1 3.9 3.8 3.6 3.5 3.4 3.2 3.1 2.9 2.7 2.6 2.4 2.3 2.1];

% The easy way: interpolate
% pp = interp1(kp_map, pp_map, kp);

% The more complicated, and maybe slightly more correct way: fit a linear
% line.  My samples are completely linear anyway, so this isn't just a
% total kludge
p = polyfit(kp_map, pp_map, 1);
pp = polyval(p, kp);