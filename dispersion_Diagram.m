clc, clear, close all,
M = -100:1:100;
N = -100:1:100;
Mod_depth = 2;
f0 = 10e9;
c = 3e8;
w0 = 2*pi*f0;
fs1 = 0.0001*f0;
ws1 = 0.1*w0;
f = (8.5:0.1:11.5)*1e9;
w = 2*pi*f;
fs = 0.1*f;
ws = 0.1*w;
%% refractive
% syms t m;
% Time = linspace(0,2/fs1,1000);
% refract_func = -1i*Mod_depth/(2*m*pi).*exp(1j*m*ws1*t);
% n = symsum(refract_func,m,-60,-1)+symsum(refract_func,m,1,60)+Mod_depth/2;
% nt = (double(vpa(subs(n,t,Time)))');
% figure,
% plot(Time,nt)
% xlabel('Time')
% ylabel('Refractive constant')
% axis([-inf inf 0 2])
%%
for in = 1: length(N)
    for im = 1 : length(N)
        if N(in)-M(im) == 0
            qp_nm = Mod_depth/2;
        else
            qp_nm = -1j*(Mod_depth/((N(in)-M(im))*pi*2));
        end
        Q(in,im) = ((qp_nm).^2)*(N(in)*ws1+w0)*(0.5*(N(in)+M(im))*ws1+w0);
    end
end
[An,Kappa_2] = eig(Q);
Kappa = (1/c)*sqrt(eig(Q));
k = real(Kappa);
%
figure,
for mode = [-2:2]
    k1 = ((w-w0)/c-mode*ws1/c);
    k2 = (-(w-w0)/c-mode*ws1/c);
    plot(f*1e-9,k1);
    hold on
    plot(f*1e-9,k2);
end
xline(f0*1e-9,'--',Color='r')
plot(f0*1e-9, -k(120)/(ws1/c),'.', 'MarkerSize', 30)
plot(f0*1e-9, -k(160)/(ws1/c),'.', 'MarkerSize', 30)
plot(f0*1e-9, k(200)/(ws1/c),'.', 'MarkerSize', 30)
plot(f0*1e-9, k(120)/(ws1/c),'.', 'MarkerSize', 30)
plot(f0*1e-9, k(160)/(ws1/c),'.', 'MarkerSize', 30)
xlabel('Freq')
ylabel('Wavenumber')
%
% alpha = [];
syms t n;
Time = linspace(0,2/fs1,1000);
alpha = symmatrix('alpha',[1 10]);
func1 = exp(1j*n*ws*t).*exp(1j*w0*t);
for i = 1 : 10
    if i == 1
        func = alpha(i).*symsum(func1,n,-100,100);
    else
        func = alpha(i).*symsum(func1,n,-100,100) + alpha(i-1).*symsum(func1,n,-100,100);
    end
end
% fx = (double(vpa(subs(func1,t,Time)))');
display(func)
%% Cal field
% alp = sym('alpha',[1,201])
for i = 1: length(M)
    An1(i,:) = sum(An(i));
end
alp = sym('alpha',[1,201]);
func1 = alp*An1 == 1;
x = solve(func1);