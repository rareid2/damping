% Compare Forrest's Matlab-based implementation vs my C++ code:

% clear all; close all;

rayfile = '/shared/users/asousa/WIPP/3dWIPP/outputs/one_ray.ray';
c_path = '/shared/users/asousa/WIPP/3dWIPP/damping/';
c_outfile = '/shared/users/asousa/WIPP/3dWIPP/outputs/c_damping.ray';
mat_path = pwd;
AE_level = [1 3];
Kp = [0, 4, 8];



% Run C++ version:
disp('Compiling C version...');
cd(c_path);
system('make');


for (aa = 1:length(AE_level))
    for (kk = 1:length(Kp))
        c_outfile = sprintf('/shared/users/asousa/WIPP/3dWIPP/outputs/c_damping_kp%d_ae%d_foust.ray',Kp(kk),AE_level(aa));
        c_cmd = sprintf('bin/damping -i %s -o %s -a %d -k %d -m %d',rayfile, c_outfile, AE_level(aa), Kp(kk), 1);
        system(c_cmd);
    end
end

c_outfile = sprintf('/shared/users/asousa/WIPP/3dWIPP/outputs/c_damping_ngo.ray');
c_cmd = sprintf('bin/damping -i %s -o %s -a %d -k %d -m %d',rayfile, c_outfile, AE_level(aa), Kp(kk), 0);
system(c_cmd);


cd(mat_path);

% Load c results:

% cdata = readrayoutput(c_outfile);

figure(1);
hold on;

% mags = zeros(length(cdata{1}.time),length(Kp)*length(AE_level));
for (aa = 1:length(AE_level))
    for (kk = 1:length(Kp))
        c_outfile = sprintf('/shared/users/asousa/WIPP/3dWIPP/outputs/c_damping_kp%d_ae%d.ray',Kp(kk),AE_level(aa));
        cdata = readrayoutput(c_outfile);
        % mags(:, kk) = cdata{1}.damping;
        plot(cdata{1}.time, cdata{1}.damping, 'color',rand(1,3),'lineWidth',3);
        text(cdata{1}.time(end), cdata{1}.damping(end),sprintf('Kp= %d, Ae= %d',Kp(kk), AE_level(aa)));
    end;
end

c_outfile = sprintf('/shared/users/asousa/WIPP/3dWIPP/outputs/c_damping_ngo.ray');
cdata = readrayoutput(c_outfile);
plot(cdata{1}.time, cdata{1}.damping, 'g','lineWidth',3);
text(cdata{1}.time(end), cdata{1}.damping(end),'Ngo');

grid on;
xlabel('Time (sec)');
ylabel('Relative power');
ylim([0 1]);