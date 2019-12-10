% Compare Forrest's Matlab-based implementation vs my C++ code:

% clear all; close all;
rayfile = "/Users/austin/Desktop/Stanford Raytracer/raytracer_osx_build/damping/test_data/example_ray_mode1.ray";
c_path = "/Users/austin/Desktop/Stanford Raytracer/raytracer_osx_build/bin";
c_outfile = "/Users/austin/Desktop/Stanford Raytracer/raytracer_osx_build/damping/test_data/example_ray_mode1.damp";
c_mode = 1;  % 0 = ngo, 1 = Foust+
mat_path = pwd;
AE_level = 1;
Kp = 8;



% Run Matlab version:

ray_pwr_matlab = damp_single('rayfile',rayfile,'DEBUGGING', 'False','steplength',0.05,'AE_level',AE_level,'kp',Kp);


% Run C++ version:
% disp('Compiling C version...');
cd(c_path);
% system('make');

% c_cmd = sprintf('bin/damping -i %s -o %s -a %d -k %d -m %d',rayfile, c_outfile, AE_level, Kp, c_mode);
c_cmd = sprintf('./damping --inp_file="%s" --out_file="%s" --Kp=%g --AE=%g --mode=%d --geom_factor=0', ...
                               rayfile,     c_outfile,     Kp,   AE_level, c_mode);

% disp(c_cmd)

% ./damping --inp_file "/Users/austin/Desktop/Stanford Raytracer/raytracer_osx_build/example_scripts/test_outputs/example_ray_mode6.ray" --out_file "/Users/austin/Desktop/Stanford Raytracer/raytracer_osx_build/example_scripts/test_outputs/example_ray_mode6.damp"  --Kp 2 --AE 1.6 --yearday 2010001 --msec 0 --geom_factor=1 --mode 0
system(c_cmd);

cd(mat_path);

% Load c results:

% cdata = readrayoutput(c_outfile);
% 12.2019 -- Austin's damping code keeps the ray and damping separate
cdata = readmatrix(c_outfile,'FileType','text');
ctime = cdata(:,2);
cdamp = cdata(:,3);
figure(1);
subplot(211);
% plot(cdata{1}.time, cdata{1}.damping,'r', cdata{1}.time, ray_pwr_matlab, 'b');
plot(ctime, cdamp, 'r', ctime, ray_pwr_matlab, 'b');
legend('C++','Matlab');
ylabel('Relative power');
subplot(212);
% plot(cdata{1}.time, abs(cdata{1}.damping - ray_pwr_matlab));
plot(ctime, abs(cdamp - ray_pwr_matlab));

xlabel('Time (sec)');
ylabel('difference');
