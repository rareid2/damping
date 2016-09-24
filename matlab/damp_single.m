function ray_power = damp_single(varargin)
    % Calculate damping for a single ray.
    % (this version for standalone compilation)
    t_start = now;

    p = inputParser;

    % Get paths
    % setup_for_raytracing;
    p.addParamValue('rayfile','/shared/users/asousa/WIPP/3dWIPP/outputs/testoutputfile.mat');
    p.addParamValue('AE_level',3);
    p.addParamValue('kp',4);
    % p.addParamValue('f',1000);  % Hz
    p.addParamValue('steplength',0.1);  % Seconds; the average step size to use
    p.addParamValue('debugging',false);
    p.parse(varargin{:});


    rayfile = p.Results.rayfile;
    AE_level = p.Results.AE_level;
    kp = p.Results.kp;
    steplength = p.Results.steplength;
    DEBUG = p.Results.debugging;



    % % Run the damping code for a single ray.
    % close all; clear all;

    % Get paths
    % setup_for_raytracing;
    % Get constants
    [R_E, c, EPS0, MU0, q, me] = physconstf;

    % Inputs
    % rayfile = fullfile(project_root,'outputs','testoutputfile.mat');
    % f = 1000;
    % AE_level = 3;
    % kp = 4;

    L_pp = kp_to_pp(kp);


    % Solver params
    TOL = 1e-3;
    % INTEGMETHOD = 'fast';
    INTEGMETHOD = 'accurate';

    all_rays = readrayoutput(rayfile);

    % ray = ray.ray_out;
    % ray.f = f; % Hz

    % w = 2*pi*f;


    % Loop over each ray in the rayfile:

    for (rr = 1:length(all_rays))
        % Select current ray in file
        ray = all_rays{rr};

        w = ray.w;

        % % ------------------- Downsample and interpolate for speed -------------------------------
        t_step = round(length(ray.time)/ray.time(end)*steplength); % Average of one step every 0.1 seconds
        t = ray.time(1:t_step:end);
        % t = ray.time;
        % disp(t_step);
        % Make sure we include the end point in the interpolation, or we'll end up
        % with NaNs in the power at the end of the vector
        if t(end) ~= ray.time(end)
            t(end+1) = ray.time(end);
        end

        % % ray_int = interpolate_ray(ray, t);

        % pos   = ray.pos(1:t_step:end,:);
        % vgrel = ray.vgrel(1:t_step:end,:);
        % n     = ray.n(1:t_step:end,:);
        % B0    = ray.B0(1:t_step:end,:);
        % qs    = ray.qs(1:t_step:end,:);
        % ms    = ray.ms(1:t_step:end,:);
        % Ns    = ray.Ns(1:t_step:end,:);
        % nus   = ray.nus(1:t_step:end,:);


        % t_step = 0.1;
        % t = 0:t_step:ray.time(end);

        %% Interpolate onto working time axis
        ray_int = interpolate_ray(ray, t);
        pos   = ray_int.pos;
        vgrel = ray_int.vgrel;
        n     = ray_int.n;
        B0    = ray_int.B0;
        qs    = ray_int.qs;
        ms    = ray_int.ms;
        Ns    = ray_int.Ns;
        nus   = ray_int.nus;


        % ignore collisions
        nus= 0*nus;

        % set up electron distribution
        % temperature
        % kT = 1e3*q; % 1 keV

        qe = -q;

        % Which resonances to do (0 = landau, ±1 = cyclotron)
        m = [0];

        magnitude = zeros(size(t));
        magnitude(1) = 1;
        kis = zeros(size(t));

        for( ii=2:length(t) )
            t_iter_start = now;

            % No one will ever notice if I use a dipole model here... (as long as I
            % calculate the plasmapause the same way)
            mag_lat = atan(pos(ii,3)/sqrt(pos(ii,1)^2 + pos(ii,2)^2));
            r = sqrt(sum(pos(ii,:).^2));
            L = r/cos(mag_lat)^2/R_E;
            MLT = mod((atan2(pos(ii,2), pos(ii,1)) + pi)/(2*pi)*24, 24); % MLT in hours; 0 MLT is in -x direction
            [n_fit, An_fit] = get_fit_params(L, MLT, AE_level, false);

            if ~isfinite(n_fit) || ~isfinite(An_fit)
                error('invalid n_fit %g or An_fit %g (L=%0.2f,MLT=%0.2f)', n_fit, An_fit, L, MLT);
            end
            fe = @(vperp, vpar) crres_polar_hybrid_psd(vperp, vpar, n_fit, An_fit, L, L_pp);
            fs = {fe};

            wce_h = ((qe*norm(B0(ii,:)))./me);
            % Hot charge
            qe_h = qe;
            % Hot mass
            me_h = me;
            % vector of hot plasma properties
            wchs = [wce_h];
            qhs = [qe_h];
            mhs = [me_h];

            % Just trust that the n in the file is the one we want to use.
            % They might in fact be off a little due to the interpolation.
            k = n(ii,:)*w/c;
            kmag = norm(k);
            % Signed component along B
            Bhat = B0(ii,:)/norm(B0(ii,:));

            kpar = k*Bhat';
            % Perpendicular vector
            kperp = k - kpar*Bhat;
            % Convert kperp into a magnitude
            kperp = norm(kperp);

            
            % fprintf('Kperp: %f, Kpar: %f\n',kperp, kpar);
            if( kmag ~= 0  )
            % Compute the spatial damping rate

            ki =  spatialdamping( fs, kperp, kpar, w, m, wchs, qhs, mhs, ...
                                  qs, Ns(ii,:), ms, nus(ii,:), ...
                                  norm(B0(ii,:)), TOL, INTEGMETHOD );
            % Take the component of ki along vg (yes, this is right!).
            ki_along_vg = ki*(k*vgrel(ii,:)')/(norm(k)*norm(vgrel(ii,:)));

            dist = norm(pos(ii,:)-pos(ii-1,:));
            kis(ii) = ki_along_vg;
            magnitude(ii) = magnitude(ii-1)*exp(-dist*ki_along_vg);

            if DEBUG
                [theta,phi,radius] = cart2sph(pos(ii,1), pos(ii,2), pos(ii,3));
                fprintf('t=%0.2f @(%0.2f, %0.2f, %0.2f), dt=%0.03f, damping=%0.2f dB, current power=%0.2f dB (iteration took %s)\n', ...
                  t(ii), phi*180/pi, theta*180/pi, radius/R_E, t(ii)-t(ii-1), db(exp(-dist*ki_along_vg)), db(magnitude(ii)), time_elapsed(t_iter_start, now));
                % fprintf('\tlat: %0.2f lon: %0.2f r: %0.2f\n',phi*180/pi, theta*180/pi, radius/R_E);
                else
                %    disp('Re{n} = 0, not solving evanescent mode');
                end;
            end
        end

        % Interpolate back to original time axis:

        % kis_interp = interp1(t, kis, ray.time);
        % mag_interp = zeros(size(kis_interp));
        % mag_interp(1) = 1;

        % for (ii=2:size(ray.time))
        %     dd = norm(ray.pos(ii,:) - ray.pos(ii-1,:));
        %     mag_interp(ii) = mag_interp(ii-1)*exp(-dd*kis_interp(ii));
        % end


        % ray.magnitude = mag_interp;

        %% Interpolate calculated power over true ray time
        % Assume that the magnitude at the true first and last time points is equal
        % to the magnitude at the calculated first and last time points (to avoid
        % NaNs when interpolating)
        if t(1) ~= ray.time(1)
            t = [ray.time(1); t];
            magnitude = [magnitude(1); magnitude];
        end
        if t(end) ~= ray.time(end)
            t = [t; ray.time(end)];
            magnitude = [magnitude; magnitude(end)];
        end

        % Interpolate in log space; power = magnitude^2!!
        ray_power = exp(interp1(t, log(magnitude), ray.time)).^2;



        % dd = sqrt(sum(ray.pos(2:end,:) - ray.pos(1:end-1,:).^2,2));
        % mag_interp = cumprod(exp(-1*dd*kis_interp(2:end)));

        if DEBUG
            figure();
            plot(t, magnitude.^2,'r',ray.time, ray_power,'b');
        end

        fprintf('damping complete -- %s\n',time_elapsed(t_start,now));


    % disp(ii);
    % fprintf('damping complete... %g seconds\n',time_elapsed(t_start, now));
    end % ray loop

end % function def