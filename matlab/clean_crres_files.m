% Clean up the CRRES data files, and fill in NaNs.

input_dir = 'crres_data';
output_dir = 'crres_data_clean';


%% Load data
% Only load the file if it's not already in memory


d = dir(fullfile(input_dir, sprintf('crres_*.mat')));

crres = {};
for kk = 1:length(d)
        [crres(kk).L, crres(kk).MLT, crres(kk).J_int, crres(kk).t_b] = ...
            get_fluxes(fullfile(input_dir, d(kk).name), 'gridfit');
        crres(kk).E = str2double(d(kk).name(8:12))/1e3; % keV
        crres(kk).AE = str2double(d(kk).name(16));
end    

squeeze(cress);

save(fullfile(output_dir, 'crres_clean.mat'), 'crres');

