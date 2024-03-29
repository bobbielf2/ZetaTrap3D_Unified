ka = 2; % wavenumber
use_unified = 1; % use the "unified" quadrature for higher-order correction
        ord = 5; % convergence order = 3, 5, 7 (unified quad only)

% cylinder parameters
cylax = [0.5,0.5,1];    % cylinder axis
r = 2.5;                % cylinder radius
%n = 20;                % #pts on circle
L = 30; a = 1/2; b = 3/4; L0 = 10;      % complexification parameters

% convergence
N = 10*round(2.^(2:0.5:6)/2);
Us = zeros(size(N));
Ud = zeros(size(N));
hh = zeros(size(N));
for ii = 1:numel(N)
    n = N(ii);
    s = cylinder_cmpl(cylax,r,L,n,a,b,L0);  % complexified cylinder

    % fix a target point
    i = 1+6*s.Nu/10;
    j = 1+34*s.Nu/10;
    ind = sub2ind([s.Nu,s.Nv],i,j);
    t.x = s.x(:,ind);

    % precompute correction weight
    if use_unified
        lptypes = {'s','d'};
        [zweis, indstens, stens] = Helm3dPointZeta_multi_cmpl(ka,ord,lptypes,s,ind);
        ind_s = indstens{1};
        Ds = zweis{1};
        ind_d = indstens{2};
        Dd = zweis{2};
    else
        ord = 3;
        [Zs,Zd] = epstein_zeta_cmpl(1,s.E(ind),s.F(ind),s.G(ind),s.e(ind),s.f(ind),s.g(ind));
        Ds = s.w(ind)/(4*pi).*(-Zs/s.h + 1i*ka);
        Dd = Zd.*s.w(ind)/(4*pi*s.h);
    end

    % evaluate potential
    sig = (sin(3*s.u+2) - 3*cos(0.7*s.v-pi)+1i*cos(2*s.u)*0.3.*s.v);% complex density
    As = Helm3dSLP_cmpl(ka,t,s); % slp
    Ad = Helm3dDLP_cmpl(ka,t,s); % dlp
    if use_unified
        As(ind) = 0;
        As(ind_s) = As(ind_s)+Ds.';
        Ad(ind) = 0;
        Ad(ind_d) = Ad(ind_d)+Dd.';
    else
        As(ind) = Ds;
        Ad(ind) = Dd;
    end
    Us(ii) = As*sig(:);
    Ud(ii) = Ad*sig(:);
    hh(ii) = s.h; % mesh size

    % plot
    if n < 40
        clf
        % plot surface
        subplot(1,4,1)
        x = real(t.x);
        plot3(x(1,:),x(2,:),x(3,:),'.r','markersize',20)
        hold on, plot_surf(s,1), hold off
        axis([-10,5,-10,5,-15,5])
        legend('= target'), title('surface (real part)')
        
        % plot density
        subplot(1,4,2)
        plot_surf(struct('x',[s.u(:),s.v(:),sig(:)]','Nu',s.Nu,'Nv',s.Nv))
        title('density (real part)')
    end
end
% convergence plot 
err_s = abs(Us-Us(end));
err_d = abs(Ud-Ud(end));
% err vs h
subplot(1,4,3) % SLP
loglog(hh,err_s,'o',hh,hh.^ord,'--'), xlabel('h')
legend({'SLP',['$O(h^',num2str(ord),')$']},'interpreter','latex')
title('SLP')
subplot(1,4,4) % DLP
loglog(hh,err_d,'o',hh,hh.^ord,'--'), xlabel('h')
legend({'DLP',['$O(h^',num2str(ord),')$']},'interpreter','latex')
title('DLP')


function plot_surf(s,wrap)
if nargin < 2, wrap = 0; end
Nr = [s.Nu,s.Nv];   % reshaping sizes
% note: assume grid generated by 'ndgrid'. If used 'meshgrid' then need to swap dim(x) and dim(y)!
x = reshape(real(s.x(1,:)),Nr);
y = reshape(real(s.x(2,:)),Nr);
z = reshape(real(s.x(3,:)),Nr);
if wrap
    x = x([1:end,1],:);
    y = y([1:end,1],:);
    z = z([1:end,1],:);
end
surf(x,y,z)
end
