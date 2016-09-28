% Compare Forrest's Matlab-based implementation vs my C++ code:

% clear all; close all;

rayfile = '/shared/users/asousa/WIPP/3dWIPP/outputs/one_ray.ray';
c_path = '/shared/users/asousa/WIPP/3dWIPP/damping/';
c_outfile = '/shared/users/asousa/WIPP/3dWIPP/outputs/c_damping.ray';
c_mode = 1;  % 0 = ngo, 1 = Foust
mat_path = pwd;
AE_level = 1;
Kp = 7.2;


% Run Matlab version:

ray_pwr_matlab = damp_single('rayfile',rayfile,'DEBUGGING', 'True','steplength',0.05,'AE_level',AE_level,'kp',Kp);


% Run C++ version:
disp('Compiling C version...');
cd(c_path);
system('make');

c_cmd = sprintf('bin/damping -i %s -o %s -a %d -k %d -m %d',rayfile, c_outfile, AE_level, Kp, c_mode);
% c_cmd = ['bin/damping -i ',rayfile,' -o ', c_outfile, ' -a ', AE_level, ' -k ', Kp];
% disp(c_cmd)
system(c_cmd);

cd(mat_path);

% Load c results:

cdata = readrayoutput(c_outfile);

figure(1);
subplot(211);
plot(cdata{1}.time, cdata{1}.damping,'r', cdata{1}.time, ray_pwr_matlab, 'b');
legend('C++','Matlab');
ylabel('Relative power');
subplot(212);
plot(cdata{1}.time, abs(cdata{1}.damping - ray_pwr_matlab));
xlabel('Time (sec)');
ylabel('difference');
