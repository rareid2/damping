function ray_int = interpolate_ray(ray_orig, new_time, str_int_idx)
% ray_int = interpolate_ray(ray_orig, new_time, str_int_idx)
% 
% Interpolate ray over new time vector
%
% If str_int_idx is 'time' (default), the new_time vector is
% intepreted as time points.
% If str_int_idx is 'index', the new_time vector is interpreted
% as indices into ray_orig; otherwise it's 'time' (default)

% By Daniel Golden (dgolden1 at stanford dot edu) March 2010
% $Id$

if ~exist('str_int_idx', 'var') || isempty(str_int_idx)
	str_int_idx = 'time';
end

if strcmp(str_int_idx, 'time')
	if new_time(end) > ray_orig.time(end)
	   error('End time (%f sec) must be less than last ray time (%f sec)', new_time(end), ray_orig.time(end)); 
	end
	
	names = fieldnames(ray_orig);
	for kk = 1:length(names)
	    if size(ray_orig.(names{kk}), 1) == 1
	        ray_int.(names{kk}) = ray_orig.(names{kk});
	    else
	        ray_int.(names{kk}) = interp1(ray_orig.time, ray_orig.(names{kk}), new_time);
	    end
	end
elseif strcmp(str_int_idx, 'index')
	if max(new_time) > length(ray_orig.time)
		error('Max index %d must be less than or equal to length of ray %d', max(new_time), length(ray_orig.time));
	end
	assert(all(islogical(new_time)) || all(fpart(new_time) == 0));

	names = fieldnames(ray_orig);
	for kk = 1:length(names)
	    if size(ray_orig.(names{kk}), 1) == 1
	        ray_int.(names{kk}) = ray_orig.(names{kk});
	    else
	        ray_int.(names{kk}) = ray_orig.(names{kk})(new_time, :);
	    end
	end

end
