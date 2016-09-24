function dec = fpart(float)
% Return the floating part of the number
% 
% For a negative number, the result is negative

% By Daniel Golden (dgolden1 at stanford dot edu) June, 2008
% $Id: fpart.m 1469 2012-02-24 01:39:12Z dgolden $

idx_pos = float >= 0;
dec(idx_pos) = float(idx_pos) - floor(float(idx_pos));
dec(~idx_pos) = float(~idx_pos) - ceil(float(~idx_pos));