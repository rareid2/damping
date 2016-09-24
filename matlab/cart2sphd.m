function varargout = cart2sphd(varargin)
% [theta, phi, rho] = cart2sphd([x y z])
% Convert cartesian to spherical coordinates, with phi measured from +z
% axis.
% 
% Matlab's cart2sph function measures phi from the x-y plane
% Inverse: sph2cartd

% By Daniel Golden (dgolden1 at stanford dot edu) February 2010
% $Id: cart2sphd.m 838 2010-03-23 21:47:11Z dgolden $

if nargin == 1
	pos_cart = varargin{1};
	x = pos_cart(:, 1);
	y = pos_cart(:, 2);
	z = pos_cart(:, 3);
elseif nargin == 3
	x = varargin{1};
	y = varargin{2};
	z = varargin{3};
end

[t, p, r] = cart2sph(x, y, z);
p = pi/2 - p; % Measure phi from +z axis instead of x-y plane

if nargout <= 1
	pos_sph = [t, p, r];
	varargout{1} = pos_sph;
else
	varargout{1} = t;
	varargout{2} = p;
	varargout{3} = r;
end