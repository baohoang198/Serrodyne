clc,
close all,
clear all,
%% Dispersion diagram (Fig. 2)
N = -60:60;
M = -60:60;
% f0 = 10e9;
f0=8:0.1:12;
f0 = f0*1e9;
w0 = 2*pi*f0;
ws = 0.1*w0;
c = 3e8;
ModulDepth = 2;
Qnm = zeros(length(M),length(N),length(f0));
for in = 1:length(N)
    for im = 1:length(M)
        for i_f = 1:length(f0)
            Mode_nm = N(in)-M(im);
            if Mode_nm == 0
                q = ModulDepth/2;
            else
                q = -1i*ModulDepth/(Mode_nm*pi);
            end
            Qnm(in,im,i_f) = (q.^2).*(N(in)*ws(i_f)+w0(i_f)).*(0.5*(N(in)+M(im)).*ws(i_f)+w0(i_f));
        end
    end
end
for i_f = 1:length(f0)
    Kappa(:,:,i_f) = (sqrt((1/c^2).*eig(Qnm(:,:,i_f))));
%     Kappa(:,:,i_f) = real(Kappa);
end
Kappa = real(Kappa);
    % c0 = 3e8;
    for i = floor(length(Kappa)/2)-2:floor(length(Kappa)/2)+2
        %     if n(i) <= 0
%         y1=((2*pi/c)*(f)-Kappa(i))/(ws/c);
        y1=reshape((Kappa(i,1,:)),1,41,1);
        plot(f0*1e-9,y1)
        hold on
        %     elseif n(i)>= 0
        %         y2=-(2*pi/c)*(f)+Kappa(i);
        %         plot(f*1e-9,y2)
        %         hold on
        %     end
    end
    xline(10,'--')
    hold off