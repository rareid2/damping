function convert_ascii_to_mat
% Convert an ASCII CRRES database file, with columns
%         MLT           L       J_perp           t_b
% To a .mat file

% By Daniel Golden (dgolden1 at stanford dot edu)
% $Id: convert_ascii_to_mat.m 1205 2011-03-04 23:15:40Z dgolden $

input_dir = '/home/dgolden/vlf/case_studies/crres_suprathermal_fluxes';
output_dir = fullfile(danmatlabroot, 'crres_suprathermal_fluxes');

d = dir(fullfile(input_dir, 'crres*.txt'));

for kk = 1:length(d)
	blah = load(fullfile(input_dir, d(kk).name));
	
	MLT = reshape(blah(:,1), length(unique(blah(:,2))), length(unique(blah(:,1))));
	L = reshape(blah(:,2), length(unique(blah(:,2))), length(unique(blah(:,1))));
	J_perp = reshape(blah(:,3), length(unique(blah(:,2))), length(unique(blah(:,1))));
	t_b = reshape(blah(:,4), length(unique(blah(:,2))), length(unique(blah(:,1))));
	
	[~, fname, ext] = fileparts(d(kk).name);
	output_filename = fullfile(output_dir, [fname '.mat']);
	save(output_filename, 'MLT', 'L', 'J_perp', 't_b');
	
	fprintf('Wrote %s\n', output_filename);
end