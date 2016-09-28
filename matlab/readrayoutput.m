function [out] = readrayoutput(filename)
% Read the raytracer output file "filename"
% Return as an array of structures with the following members:
% 
% For a ray index i:
%
% out{i}.stopcond - Reason for stopping (see raytracer.f95 for an explanation)
% out{i}.time - n x 1 vector of times for each timestep
% out{i}.pos - n x 3 vector of positions, in meters for each timestep
% out{i}.vprel - n x 3 vector of relative (to c) phase velocities
% out{i}.vgrel - n x 3 vector of relative (to c) group velocities
% out{i}.n - n x 3 vector of refractive index vectors
% out{i}.B0 - n x 3 vector of background magnetic field (Tesla)
% out{i}.Nspec - Number of plasma species
% out{i}.qs - n x Nspec vector of charges in Coulombs
% out{i}.ms - n x Nspec vector of masses in kg
% out{i}.Ns - n x Nspec vector of number densities in m^-3
% out{i}.nus - n x Nspec vector of collision frequencies in s^-1
% 
% Updated to v.1.16 , APS 8.12.2016 (added frequency in ray log)

  %out{1}.pos 
  % Read the file
  x = load(filename,'-ascii');
  % fields: raynum, stopcond, time, pos(3), vprel(3), vgrel(3), n(3),
  % B0(3), Nspec, qs(Nspec), ms(Nspec), Ns(Nspec), nus(Nspec)
  
  
  rayindexes = x(:,1);
  numrays = length(unique(rayindexes)); %max(rayindexes);

  for( ii=1:numrays )
    tmp = x(rayindexes==ii,:);
    out{ii}.stopcond = tmp(1,2);
    out{ii}.time = tmp(:,3);
    out{ii}.pos = tmp(:,4:6);
    out{ii}.vprel = tmp(:,7:9);
    out{ii}.vgrel = tmp(:,10:12);
    out{ii}.n = tmp(:,13:15);
    out{ii}.B0 = tmp(:,16:18);
    % The number of species is the same for each timestep
    out{ii}.w     = tmp(1,19);
    out{ii}.Nspec = tmp(1,20);
    out{ii}.qs =  tmp(:,(21+0*out{ii}.Nspec):(21+1*out{ii}.Nspec-1));
    out{ii}.ms =  tmp(:,(21+1*out{ii}.Nspec):(21+2*out{ii}.Nspec-1));
    out{ii}.Ns =  tmp(:,(21+2*out{ii}.Nspec):(21+3*out{ii}.Nspec-1));
    out{ii}.nus = tmp(:,(21+3*out{ii}.Nspec):(21+4*out{ii}.Nspec-1));

    if (size(tmp) > (21 + 4*out{ii}.Nspec - 1))
      disp('loading damping...');
      out{ii}.damping = tmp(:, end);
    end;

    % Turn vectors of constants into scalars
    out{ii}.qs = out{ii}.qs(1,:);
    out{ii}.ms = out{ii}.ms(1,:);
  end;
