% Compare the effect of cyclotron and higher-order resonance modes:
clear all; close all;

f0 = '/shared/users/asousa/WIPP/3dWIPP/outputs/foust_0_0.ray';
f1 = '/shared/users/asousa/WIPP/3dWIPP/outputs/foust_1_1.ray';
f2 = '/shared/users/asousa/WIPP/3dWIPP/outputs/foust_10_10.ray';


d0 = readrayoutput(f0);
d1 = readrayoutput(f1);
d2 = readrayoutput(f2);


figure(1);
hold on;
plot(d0{1}.time, d0{1}.damping,'r');
plot(d1{1}.time, d1{1}.damping,'g');
plot(d2{1}.time, d2{1}.damping,'b');

legend('landau','landau + cyclotron','-10..10');
xlabel('time (sec)');
% ylim([0 1]);