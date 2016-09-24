function r_center = centroid(r, m)
% centr = centroid(r, m)
% 
% Function to find the 1-d centroid of a vector
% 
% INPUTS
% r: distance of points from reference points (x-values)
% m: mass of points
% 
% OUTPUTS
% r_center: r value of the center of mass
% 
% If m is a matrix, a centroid will be returned for each column
% 
% This is the same as a weighted average, with values r and weights m

% By Daniel Golden (dgolden1 at stanford dot edu)
% $Id: centroid.m 973 2010-05-28 01:06:38Z dgolden $

if isvector(r)
	r = r(:);
end
if isvector(m)
	m = m(:);
end

if ~all(size(r) == size(m))
    r = repmat(r, 1 + size(m, 1) - size(r, 1), 1 + size(m, 2) - size(r, 2));
    if ~all(size(r) == size(m))
        error('size r [%d %d] ~= size m [%d %d]', size(r), size(m));
    end
end

r_center = sum(r.*m)./(sum(m));
