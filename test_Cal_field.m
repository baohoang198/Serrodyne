clc, clear all, close all
%% Cal field
c = 3e8;
f0 = 10e9;
w0 = 2*pi*f0;
ws = 0.0001*w0;
fs = 0.0001*f0;
n = -5:1:5;
w = w0+n*ws;
f = w./(2*pi);
%
Mod_depth =2;
syms t m;
Time = linspace(0,1/fs,1000);
refract_func = -1i*Mod_depth/(2*m*pi).*exp(1j*m*ws*t);
n = symsum(refract_func,m,-60,-1)+symsum(refract_func,m,1,60)+Mod_depth/2;
nt = (double(vpa(subs(n,t,Time)))');
% figure,
% plot(Time,nt)
% xlabel('Time')
% ylabel('Refractive constant')
% axis([-inf inf 0 2])
for i = 1:length(n)
    fx(i,:) = f-(f0+n(i)*fs);
    y(i,:) = dirac(fx(i,:));
end
idx =y==Inf;
y(idx) = 1;
%
% figure,
% for i = 1:length(n)
% stem
% end
%
harm = -5:1:5;
d = 0.01;
E_out = exp(-1*w0*d/c*real(nt));
syms t;
for i = 1:length(harm)
 fnc = double(vpa(int(exp(-1j*(harm(i)*ws*t + w0*(d/c)*n)),[-0.5/fs,0.5/fs])));
%  fnc = double(vpa(int(exp(-1j*(harm(i)*ws*t + ws*t)),[-0.5/fs,0.5/fs])));
%  fnc = double(vpa(int(exp(-1j*(w0*d/c*n)),[-0.5/fs,0.5/fs])))
 E(i) = fnc;
end
figure,
stem(harm,abs(E)./max(abs(E)))
ylabel('|Eout/Ein|^2')
xlabel('Harmonic')
title(['M = ',Mod_depth])
% E = sum(E);