ka = 2;

%N = 32*round(2.^(1.5:0.5:6)/2);
N = 32*round(2.^(2:0.5:6)/2);
Us = zeros(size(N));
Ud = zeros(size(N));
hh = zeros(size(N));
for ii = 1:numel(N)
    % complex path
    % n = 64;
    n = N(ii);
    L = 30;
    a = 1/2; b = 3/4; L0 = 10;
    s = patch_cmpl(L,n,a,b,L0);

    % fix a target point
    i = 1+12*n/32;
    j = 1+11*n/32;
    ind = sub2ind(n+1*[1,1],i,j);
    t.x = s.x(:,ind);

    % precompute correction weight (1-point)
    [Zs,Zd] = epstein_zeta_cmpl(1,s.E(ind),s.F(ind),s.G(ind),s.e(ind),s.f(ind),s.g(ind));
    Ds = s.w(ind)/(4*pi).*(-Zs/s.h + 1i*ka);
    Dd = Zd.*s.w(ind)/(4*pi*s.h);

    % evaluate potential
    sig = (sin(0.6*s.u+2) - 3*cos(0.7*s.v-pi)).*exp(-0.05*(s.u.^2+s.v.^2)); % density
    A = Helm3dSLP_cmpl(ka,t,s); A(ind) = Ds; % slp
    u = A*sig(:);
    Us(ii) = u;
    A = Helm3dDLP_cmpl(ka,t,s); A(ind) = Dd; % dlp
    u = A*sig(:);
    Ud(ii) = u;
    hh(ii) = s.h; % mesh size

    if n < 70
        clf
        % plot surface
        subplot(1,3,1)
        plot3(real(t.x(1)),real(t.x(2)),t.x(3),'.r','markersize',30)
        hold on, plot_surf(s), hold off
        legend('= target'), title('surface (real part)')
        % plot density
        subplot(1,3,2)
        plot_surf(struct('x',[s.u(:),s.v(:),sig(:)]','Nu',n+1,'Nv',n+1))
        title('density')
    end
end
% convergence plot
err_s = abs(Us-Us(end));
err_d = abs(Ud-Ud(end));
subplot(1,3,3)
% err vs N
%loglog(N,err_s,'*',N,err_d,'o',N,1*N.^-3,'--','LineWidth',1), xlabel('$n$','Interpreter','latex')
% err vs h
loglog(hh,err_s,'*',hh,err_d,'o',hh,1e-5*hh.^3,'--','LineWidth',1), xlabel('$h$','Interpreter','latex') % ylim([8e-11,1e-2])
legend('SLP','DLP','$O(h^3)$','interpreter','latex','Location','northwest')
title('error')
%% annotations
annotation(gcf,'textbox',...
    [0.452903242672827 0 0.213023505384887 0.110687022900762],...
    'String','($\sigma = 0.0017$ at target)',...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',14,...
    'FitBoxToText','off');
annotation(gcf,'textbox',...
    [0.0267101239158896 0.0036871669307216 0.359527389957584 0.110687022900762],...
    'String','target $x=[-7.5-0.001i, -9.375-0.108i,-0.004]$',...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',14,...
    'FitBoxToText','off');

function plot_surf(s)
Nr = [s.Nu,s.Nv];   % reshaping sizes
% note: assume grid generated by 'ndgrid'. If used 'meshgrid' then need to swap dim(x) and dim(y)!
x = reshape(real(s.x(1,:)),Nr);
y = reshape(real(s.x(2,:)),Nr);
z = reshape(real(s.x(3,:)),Nr);
surf(x,y,z)
end