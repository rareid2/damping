function varargout = sph2cartd(varargin)
% [x, y, z] = sph2cartd([theta phi rho])
% Convert spherical to cartesian coordinates, with phi measured from +z
% axis.
% 
% Matlab's sph2cart function measures phi from the x-y plane
% Inverse: cart2sphd

% By Daniel Golden (dgolden1 at stanford dot edu) February 2010
% $Id: sph2cartd.m 834 2010-03-22 23:38:03Z dgolden $
if nargin == 1
	pos_sph = varargin{1};
	t = pos_sph(:, 1);
	p = pos_sph(:, 2);
	r = pos_sph(:, 3);
elseif nargin == 3
	t = varargin{1};
	p = varargin{2};
	r = varargin{3};
end

[x, y, z] = sph2cart(t, pi/2 - p, r);

if nargout <= 1
	pos_cart = [x y z];
	varargout{1} = pos_cart;
else
	varargout{1} = x;
	varargout{2} = y;
	varargout{3} = z;
end