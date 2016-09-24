% script to plot all fluxes in original and gridfitted forms

% By Daniel Golden (dgolden1 at stanford dot edu) May 2010
% $Id: plot_all_fluxes.m 1205 2011-03-04 23:15:40Z dgolden $

fnames = {'crres_e00213_ae3', 'crres_e01090_ae3', 'crres_e04250_ae3', 'crres_e16500_ae3'};
figure;

%% Vertical
% for kk = 1:4
%   plot_fluxes(fnames{kk}, 'fluxes', 'none', subplot(4, 2, kk*2-1));
%   plot_fluxes(fnames{kk}, 'fluxes', 'gridfit', subplot(4, 2, kk*2));
% end
% 
% figure_grow(gcf, 0.7, 2);
% increase_font(gcf, 12);

%% Horizontal
for kk = 1:4
  plot_fluxes(fnames{kk}, 'fluxes', 'none', subplot(2, 4, kk));
  if kk == 1
    ylabel('Original');
  end
  plot_fluxes(fnames{kk}, 'fluxes', 'gridfit', subplot(2, 4, 4+kk)); title('')
  if kk == 1
    ylabel('Smoothed');
  end
  
end

figure_grow(gcf, 1.6, 1);
increase_font(gcf, 14);
