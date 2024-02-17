% clc
% clear all,
% close all,
% c0 = 3e8;
% f=0:0.1:12;
% f = f*1e9;
% n = 0.2:0.1:2;
% n = flip(n);
% for i = 1:length(n)
% y1=(2*pi/c0)*n(i)*f;
% plot(f,y1)
% hold on
% y2=-(2*pi/c0)*n(i)*f;
% plot(f,y2)
% end
% hold off

clc
clear all,
close all,
c0 = 3e8;
f=8:0.1:12;
f = f*1e9;
fs = 1e9;
n1 = [-10:1:-5]
n2 = [5:1:10]
n = [n1,n2]
% n = [-10:1:10];
period = c0/fs;
for i = 1:length(n)
    if n(i) <= 0
        y1=(2*pi/c0)*(f)+2*n(i)*pi/period;
        plot(f*1e-9,y1)
        hold on
    elseif n(i)>= 0
        y2=-(2*pi/c0)*(f)+2*n(i)*pi/period;
        plot(f*1e-9,y2)
        hold on
    end 
end
hold off