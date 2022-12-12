% Quasi-periodic scattering in 3D
% This code solves a 2-layered media transmission problem
% using boundary integral equation formulation
% with a periodizing scheme. See [1] for more details.
% 
% [1] B. Wu and M.H. Cho, arXiv:2211.14656 (2022)

addpath ../multilayermedia/
addpath ../potentials/
addpath ../utils/
addpath ../zetafunc/

%----------------------------------------
% Set up geometry & collocation pts
h = 0; % h=0 is flat interface, h>0 is corrugate inteface
Nu = 40; Nv = 40;   % # of pts on periodic surface
N = Nu*Nv;          % total # of pts
Mw = 20;            % # of gauss pts on L, R, F, B walls
Mu = 20;            % # of pts on U wall is Mu^2
P = 30;     % mesh order on proxy sphere ( (P-1)x(2P) mesh )
Rp = 1.5;     % proxy sphere radius
K = 10;     % use 2K+1 terms in Rayleigh-Bloch expansion

s = periSurf2D(Nu,Nv,h); % scattering surface

% Set up side walls, top & bottom walls, proxy spheres ...
% and plot them
clf, subplot(1,2,1)
[Wl,Wr,Wf,Wb,Wu,Wd,px,px2] = periWall2D(Mw,Mu,P,Rp,s);

%----------------------------------------
% Set scattering params
% Incident wave vector = ka1*(cos(th)*cos(ph), cos(th)*sin(ph), sin(th))
% where incident angle -pi < th < 0, phase angle -pi/2 <= ph <= pi/2
ka1 = 5;               	% wavenumbers
ka2 = 10;
%th = -acos(1-2*pi/ka);	% incident angle (Wood anomaly: th = -acos(1-2*pi*n/s.d/ka) for any integer n)
th = -pi/3;
ph = pi/4;              % phase angle
kv = ka1*[cos(th)*cos(ph);cos(th)*sin(ph);sin(th)]; % wave vector
alph_x = exp(1i*kv(1)*s.du);  % Bloch phases (quasi-periodicity parameter)
alph_y = exp(1i*kv(2)*s.dv);
s.alph_x=alph_x; s.alph_y=alph_y; % allow zeta correction account for Bloch phases

% Transmission BC
f = -exp(1i*sum(kv.*s.x).'); % f=[ -u^inc; -d/dn{u^inc} ]
f = [f;-exp(1i*sum(kv.*s.x).').*(1i*sum(kv.*s.nx).')];

%----------------------------------------
% Precompute zeta weights
fprintf('precomputing zeta weights... ')
tic
ord = 5; % zeta quadr order
ZCs = Helm3dPatchZetaSparse_transmit(ka1,ka2,ord,'s',s); % zeta correction for SLP
ZCd = Helm3dPatchZetaSparse_transmit(ka1,ka2,ord,'d',s); % zeta correction for DLP
ZCsn = Helm3dPatchZetaSparse_transmit(ka1,ka2,ord,'sn',s); % zeta correction for SLP
ZCdn = Helm3dPatchZetaSparse_transmit(ka1,ka2,ord,'dn',s); % zeta correction for DLP
fprintf('%.2f sec\n',toc)

%----------------------------------------
% Fill system matrix
disp('filling matrices...')
% potential on surface (incl contribution from nearest neighbors)
tic
A = 0;
for lu = -1:1
    for lv = -1:1
        dl = [s.du*lu;s.dv*lv;0];
        [D1,Dn1] = Helm3dDLPmat(ka1,s,s,dl);
        [S1,Sn1] = Helm3dSLPmat(ka1,s,s,dl);
        [D2,Dn2] = Helm3dDLPmat(ka2,s,s,dl);
        [S2,Sn2] = Helm3dSLPmat(ka2,s,s,dl);
        A = A + alph_x^lu*alph_y^lv*[ D1-D2,    S1-S2;
                                     Dn1-Dn2, Sn1-Sn2];
    end
end
A = A + kron([1,0;0,-1],speye(N))+[ZCd,ZCs;ZCdn,ZCsn];
fprintf('\tA matrix %.2f sec\n',toc)
clear D1 D2 S1 S2 Dn1 Dn2 Sn1 Sn2
% proxy to surface
tic
[B11,Bn11] = Helm3dDSLPmat(ka1,s,px);
[B12,Bn12] = Helm3dDSLPmat(ka2,s,px2);
B = [ B11, -B12;
     Bn11, -Bn12];
fprintf('\tB matrix %.2f sec\n',toc)
clear B11 B12 Bn11 Bn12
% surface to left/right & front/back discrepancy
tic
C11rl = 0; C11bf = 0; C21rl = 0; C21bf = 0;
for l = -1:1
    % left/right discrepancy
    dl_r = [-s.du;l*s.dv;0]; dl_l = [s.du;l*s.dv;0];
        % layer 1
    [Crd, Crdn] = Helm3dDLPmat(ka1,Wr{1},s,dl_r);
    [Crs, Crsn] = Helm3dSLPmat(ka1,Wr{1},s,dl_r);
    [Cld, Cldn] = Helm3dDLPmat(ka1,Wl{1},s,dl_l);
    [Cls, Clsn] = Helm3dSLPmat(ka1,Wl{1},s,dl_l);
    C11rl = C11rl + alph_y^l*[alph_x^-1*Crd - alph_x^2*Cld,   alph_x^-1*Crs - alph_x^2*Cls;
                              alph_x^-1*Crdn - alph_x^2*Cldn, alph_x^-1*Crsn - alph_x^2*Clsn];
        % layer 2
    [Crd, Crdn] = Helm3dDLPmat(ka2,Wr{2},s,dl_r);
    [Crs, Crsn] = Helm3dSLPmat(ka2,Wr{2},s,dl_r);
    [Cld, Cldn] = Helm3dDLPmat(ka2,Wl{2},s,dl_l);
    [Cls, Clsn] = Helm3dSLPmat(ka2,Wl{2},s,dl_l);
    C21rl = C21rl + alph_y^l*[alph_x^-1*Crd - alph_x^2*Cld,   alph_x^-1*Crs - alph_x^2*Cls;
                              alph_x^-1*Crdn - alph_x^2*Cldn, alph_x^-1*Crsn - alph_x^2*Clsn];
    % front/back discrepancy
    dl_b = [l*s.du;-s.dv;0]; dl_f = [l*s.du;s.dv;0];
        % layer 1
    [Cbd, Cbdn] = Helm3dDLPmat(ka1,Wb{1},s,dl_b);
    [Cbs, Cbsn] = Helm3dSLPmat(ka1,Wb{1},s,dl_b);
    [Cfd, Cfdn] = Helm3dDLPmat(ka1,Wf{1},s,dl_f);
    [Cfs, Cfsn] = Helm3dSLPmat(ka1,Wf{1},s,dl_f);
    C11bf = C11bf + alph_x^l*[alph_y^-1*Cbd - alph_y^2*Cfd,   alph_y^-1*Cbs - alph_y^2*Cfs;
                              alph_y^-1*Cbdn - alph_y^2*Cfdn, alph_y^-1*Cbsn - alph_y^2*Cfsn];
        % layer 2
    [Cbd, Cbdn] = Helm3dDLPmat(ka2,Wb{2},s,dl_b);
    [Cbs, Cbsn] = Helm3dSLPmat(ka2,Wb{2},s,dl_b);
    [Cfd, Cfdn] = Helm3dDLPmat(ka2,Wf{2},s,dl_f);
    [Cfs, Cfsn] = Helm3dSLPmat(ka2,Wf{2},s,dl_f);
    C21bf = C21bf + alph_x^l*[alph_y^-1*Cbd - alph_y^2*Cfd,   alph_y^-1*Cbs - alph_y^2*Cfs;
                              alph_y^-1*Cbdn - alph_y^2*Cfdn, alph_y^-1*Cbsn - alph_y^2*Cfsn];
end
C = [C11rl;C11bf;C21rl;C21bf];
fprintf('\tC matrix %.2f sec\n',toc)
clear C11rl C11bf C21rl C21bf 
clear Crd Crdn Crs Crsn Cld Cldn Cls Clsn Cbd Cbdn Cbs Cbsn Cfd Cfdn Cfs Cfsn
% proxy to left/right & front/back discrepancies
tic
[Qrds, Qrdsn] = Helm3dDSLPmat(ka1,Wr{1},px);
[Qlds, Qldsn] = Helm3dDSLPmat(ka1,Wl{1},px);
[Qbds, Qbdsn] = Helm3dDSLPmat(ka1,Wb{1},px);
[Qfds, Qfdsn] = Helm3dDSLPmat(ka1,Wf{1},px);
Q1 = [Qrds - alph_x*Qlds;
     Qrdsn - alph_x*Qldsn;
     Qbds - alph_y*Qfds;
     Qbdsn - alph_y*Qfdsn];
[Qrds, Qrdsn] = Helm3dDSLPmat(ka2,Wr{2},px2);
[Qlds, Qldsn] = Helm3dDSLPmat(ka2,Wl{2},px2);
[Qbds, Qbdsn] = Helm3dDSLPmat(ka2,Wb{2},px2);
[Qfds, Qfdsn] = Helm3dDSLPmat(ka2,Wf{2},px2);
Q2 = [Qrds - alph_x*Qlds;
     Qrdsn - alph_x*Qldsn;
     Qbds - alph_y*Qfds;
     Qbdsn - alph_y*Qfdsn];
Q = blkdiag(Q1,Q2);
fprintf('\tQ matrix %.2f sec\n',toc)
clear Q1 Q2 Qrds Qrdsn Qlds Qldsn Qbds Qbdsn Qfds Qfdsn
% surface to upper & lower wall
tic
ZU = 0; ZD = 0;
for lu = -1:1
    for lv = -1:1
        dl = [s.du*lu;s.dv*lv;0];
        % upper wall
        [Zd, Zdn] = Helm3dDLPmat(ka1,Wu,s,dl);
        [Zs, Zsn] = Helm3dSLPmat(ka1,Wu,s,dl);
        ZU = ZU +alph_x^lu*alph_y^lv*[Zd, Zs; Zdn, Zsn];
        % lower wall
        [Zd, Zdn] = Helm3dDLPmat(ka2,Wd,s,dl);
        [Zs, Zsn] = Helm3dSLPmat(ka2,Wd,s,dl);
        ZD = ZD +alph_x^lu*alph_y^lv*[Zd, Zs; Zdn, Zsn];
    end
end
Z = [ZU; ZD];
clear ZU ZD Zd Zdn Zs Zsn
% proxy to upper wall
[Vds, Vdsn] = Helm3dDSLPmat(ka1,Wu,px);
VU = [Vds; Vdsn];
[Vds, Vdsn] = Helm3dDSLPmat(ka2,Wd,px2);
VD = [Vds; Vdsn];
V = blkdiag(VU,VD);
clear VU VD Vd Vdn Vs Vsn
% Rayleigh-Bloch expansion
[m1,m2] = meshgrid(-K:K); m1=m1(:)'; m2=m2(:)';
km = kv(1:2)+2*pi./[s.du;s.dv].*[m1;m2];    % horizontal wavenumbers
% upper layer
kmperpU = sqrt(ka1^2 - sum(km.^2,1));         % vertical wavenumbers
WU = -exp(1i*(km(1,:).*Wu.x(1,:)'+km(2,:).*Wu.x(2,:)'));
WU = [WU; 1i*kmperpU.*WU];
% lower layer
kmperpD = sqrt(ka2^2 - sum(km.^2,1));         % vertical wavenumbers
WD = -exp(1i*(km(1,:).*Wd.x(1,:)'+km(2,:).*Wd.x(2,:)'));
WD = [WD; -1i*kmperpD.*WD];
W = blkdiag(WU,WD);
clear WU WD
fprintf('\tZ,V,W matrices %.2f sec\n',toc)

%----------------------------------------
% Solve linear system
disp('solving linear system...')
tic
% Schur complement of Q
CZ = [C;Z]; QW = [Q,zeros(2*4*Mw^2,2*(2*K+1)^2);V,W];
QdagC = linsolve(QW,CZ,struct('RECT',true));
B0 = [B,zeros(2*N,2*(2*K+1)^2)];
Aper = A - B0*QdagC;
tau = Aper\f;
c_a = -QdagC*tau;
abs_res = norm([A*tau+B0*c_a-f;CZ*tau+QW*c_a]);
fprintf('\trel residual = %.2g\n',abs_res/norm(f))
c1 = c_a(1:px.n);
c2 = c_a(px.n+1:px.n+px2.n);
aU = c_a(px.n+px2.n+1:px.n+px2.n+(2*K+1)^2);
aD = c_a(px.n+px2.n+(2*K+1)^2+1:end);
fprintf('\tsolve time %.2f sec\n',toc)

% Conservation of flux via Bragg's coeffs (i.e. the "a" vector)
eng_bragg = (kmperpU>0)*(abs(aU).^2 .* kmperpU') + (kmperpD>0)*(abs(aD).^2 .* kmperpD');
eng_ka = -ka1*sin(th);
err_flux = abs((eng_bragg - eng_ka)/eng_ka);
fprintf('\trel flux err = %.2g\n',err_flux)

%----------------------------------------
% Accuracy at test point
if h == 0 % flat interface, exact solution known
    % set up exact solution
    th2 = -acos(ka1/ka2*cos(th));                          % refractive angle
    kv2 = ka2*[cos(th2)*cos(ph);cos(th2)*sin(ph);sin(th2)]; % refractive wave number
    rs=(ka1*sin(th)-ka2*sin(th2))/(ka1*sin(th)+ka2*sin(th2)); % reflection strength
    ts = 1+rs; % transmission strength
    
    % reflected wave
    t1.x = [-0.25;-0.25;0.25]; % test point
    u1 = Helm3dDSLPmat(ka1,t1,px)*c1; % proxy contribution
    for lu = -1:1    % add near surface potential contribution
        for lv = -1:1
            dl = [s.du*lu;s.dv*lv;0];
            u1=u1+alph_x^lu*alph_y^lv*(...
                Helm3dDLPmat(ka1,t1,s,dl)*tau(1:end/2) ...
                + Helm3dSLPmat(ka1,t1,s,dl)*tau(end/2+1:end));
        end
    end
    % exact solution at test point if flat surface
    uref1 = rs*exp(1i*(kv(1)*t1.x(1)+kv(2)*t1.x(2)-kv(3)*t1.x(3)));
    err_reflect = norm(u1-uref1)/norm(uref1); 
    fprintf('\trel err in u (reflected) = %.2g\n',err_reflect)
    
    % transmitted wave
    t2.x = [-0.25;-0.25;-0.25]; % test point
    u2 = Helm3dDSLPmat(ka2,t2,px2)*c2; % proxy contribution
    for lu = -1:1    % add near surface potential contribution
        for lv = -1:1
            dl = [s.du*lu;s.dv*lv;0];
            u2=u2+alph_x^lu*alph_y^lv*(...
                Helm3dDLPmat(ka2,t2,s,dl)*tau(1:end/2) ...
                + Helm3dSLPmat(ka2,t2,s,dl)*tau(end/2+1:end));
        end
    end
    % exact solution at test point if flat surface
    uref2 = ts*exp(1i*(kv2(1)*t2.x(1)+kv2(2)*t2.x(2)+kv2(3)*t2.x(3)));
    err_transmit = norm(u2-uref2)/norm(uref2); 
    fprintf('\trel err in u (transmitted) = %.2g\n',err_transmit)
end

%----------------------------------------
% Evaluate potential field
disp('evaluating potential field...')
% set up plotting mesh
Nt = 100;
x = (-Nt/2:Nt/2-1)/Nt*s.du; z = (1:Nt)/Nt*2-1;
[xx,zz] = meshgrid(x,z);
% find target pts in each layer
tmp = s.Z(xx(:)',0*xx(:)'); tmp = reshape(tmp(3,:),Nt,Nt);
ind1 = zz > tmp+0.02;
ind2 = zz < tmp-0.02;
% evaluate field in each layer
u = nan(size(xx))*(1+1i);
t1.x = [xx(ind1).*[1,0],zz(ind1)]';
t2.x = [xx(ind2).*[1,0],zz(ind2)]';
u(ind1) = Helm3dDSLPmat(ka1,t1,px)*c1;  % add proxy correction
u(ind2) = Helm3dDSLPmat(ka2,t2,px2)*c2;
for lu = -1:1    % add near surface potential contribution
    for lv = -1:1
        dl = [s.du*lu;s.dv*lv;0];
        u(ind1)=u(ind1)+alph_x^lu*alph_y^lv*( ...
                         Helm3dDLPmat(ka1,t1,s,dl)*tau(1:end/2) ...
                       + Helm3dSLPmat(ka1,t1,s,dl)*tau(end/2+1:end));
        u(ind2)=u(ind2)+alph_x^lu*alph_y^lv*( ...
                         Helm3dDLPmat(ka2,t2,s,dl)*tau(1:end/2) ...
                       + Helm3dSLPmat(ka2,t2,s,dl)*tau(end/2+1:end));
    end
end
u(ind1) = u(ind1) + exp(1i*sum(kv.*t1.x).'); % add incident wave

% Plot solution field (& neighboring fields)
subplot(1,2,2);
uu = [alph_x^-1*u,u,alph_x*u]; % include quasi-periodic neighbors
imagesc([x-s.du,x,x+s.du],z,real(uu)), axis xy equal tight
colorbar, colormap jet
title('quasi-periodic field (real part) in yz-plane')

%-------------------------------------------------------------------------

function [Wl,Wr,Wf,Wb,Wu,Wd,px,px2] = periWall2D(Mw,Mu,P,Rp,s,ifplot)
% setup proxy points & periodic wall
% Input:
%   Mw: such that Mw^2 points are used on each periodic wall
%   Mu: such that Mu^2 points are used on upper wall
%   P: order of proxy-sphere quadrature
%   Rp: proxy-sphere radius
%   s: scattering surface struct, see periSurf2D.m
%   ifplot: 0=no plots, 1=plot periodization setup
% Output:
%   Wl,Wr,Wf,Wb: left,right,front,back walls, on which quasi-periodic
%                condition is imposed
%   Wu: upper wall, on which radiation condition is imposed
%   px: proxy struct (points, normals, & weights)

if nargin < 6, ifplot = 1; end
%t = gauss(Mw); t= t(:)'; % t in [-1,1]
t = linspace(-1,1,Mw); t = t(:)'; % uniform grid
% build left/right walls (gauss collocation)
Wl = cell(2,1);
L0 = s.Z(-s.du/2+0*t,s.dv/2*t); L0 = L0(3,:);   % height of left bdry of scattering surf
[y,z]=meshgrid(t/2*s.dv,(t+1)/2); 
zz1=z.*(1-L0)+L0;                               % left wall, layer 1
Wl{1}.x = [-s.du/2*ones(numel(y),1),y(:),zz1(:)]'; 
Wl{1}.nx = [1;0;0]+0*Wl{1}.x;
zz2=z.*(L0+1)-1;                                % left wall, layer 2
Wl{2}.x = [-s.du/2*ones(numel(y),1),y(:),zz2(:)]'; 
Wl{2}.nx = [1;0;0]+0*Wl{2}.x;
Wr = Wl; 
Wr{1}.x(1,:) = Wl{1}.x(1,:)+s.du;         	% right wall, layer 1
Wr{2}.x(1,:) = Wl{2}.x(1,:)+s.du;         	% right wall, layer 2
% build front/back walls (gauss collocation)
Wf = cell(2,1);
F0 = s.Z(s.du/2*t,-s.dv/2+0*t); F0 = F0(3,:);   % height of front bdry of scattering surf
[x,z]=meshgrid(t/2*s.du,(t+1)/2); 
zz1=z.*(1-F0)+F0;                               % front wall, layer 1
Wf{1}.x = [x(:),-s.dv/2*ones(numel(x),1),zz1(:)]';   
Wf{1}.nx = [0;1;0]+0*Wf{1}.x;
zz2=z.*(F0+1)-1;                                % front wall, layer 2
Wf{2}.x = [x(:),-s.dv/2*ones(numel(x),1),zz2(:)]';   
Wf{2}.nx = [0;1;0]+0*Wf{2}.x;
Wb = Wf;
Wb{1}.x(2,:) = Wf{1}.x(2,:)+s.dv;         	% back wall, layer 1
Wb{2}.x(2,:) = Wf{2}.x(2,:)+s.dv;         	% back wall, layer 2
% build top wall (uniform collocation)
[x,y]=meshgrid((-Mu/2+0.5:Mu/2-0.5)'/Mu*s.du,(-Mu/2+0.5:Mu/2-0.5)'/Mu*s.dv);
Wu.x =[x(:),y(:),ones(numel(x),1)]';            % top wall
Wu.nx=[0;0;1]+0*Wu.x;
% build bottom wall
Wd = Wu;
Wd.x(3,:) = -1;
% build proxy pts (spherical grid)
if isempty(Rp), Rp = s.d*1.5; end
if ischar(P)        % Lebedev points
    [px.nx,w] = lebquad(P);             % proxy normal
    px.x = [0;0;0.4]+Rp*px.nx;      % proxy pts
    px.w = ones(size(px.x(1,:)));	% dummy wei
    %px.w = w;                       % lebedev wei
    px.n = numel(px.w);
elseif isfloat(P)  	% Gauss-Fourier mesh
    u = (0:2*P-1)/P*pi;
    %[v,w] = gauss(P); v = v(:);
    %X = cos(u).*sqrt(1-v.^2); Y = sin(u).*sqrt(1-v.^2); Z = ones(size(u)).*v;
    v = (1:P-1)/P*pi; v = v(:);
    X = cos(u).*sin(v); Y = sin(u).*sin(v); Z = ones(size(u)).*cos(v);
    px.nx = [X(:),Y(:),Z(:)]';              % proxy normals
    px.x = [0;0;0.4]+Rp*px.nx;              % proxy pts
    %px.w = repmat(Rp^2*2*pi/(2*P+1)*w(:)',1,2*P+1);	% proxy wei
    px.w = ones(size(px.x(1,:)));       % dummy wei
    px.n = numel(px.w);
end
px2 = px; px2.x(3,:) = px2.x(3,:) - 1; % proxy for the second layer
if ifplot % plot quadrature/collocation points
    % plot 2-peri surface
    N = [s.Nu,s.Nv];
    X = reshape(s.x(1,:),N); Y = reshape(s.x(2,:),N); Z = reshape(s.x(3,:),N);
    XX = [X-s.du,X-s.du,X-s.du;X,X,X;X+s.du,X+s.du,X+s.du];
    YY = [Y-s.dv,Y,Y+s.dv;Y-s.dv,Y,Y+s.dv;Y-s.dv,Y,Y+s.dv];
    ZZ = kron(ones(3),Z);
    mesh(XX,YY,ZZ,'edgecolor','k'), hold on
    % plot proxies
    X = px.x(1,:); Y = px.x(2,:); Z = px.x(3,:);
    plot3(X,Y,Z,'.')
    X = px2.x(1,:); Y = px2.x(2,:); Z = px2.x(3,:);
    plot3(X,Y,Z,'.')
    % plot L,R,F,B,U,D walls
    for i = 1:2
    X = [Wl{i}.x(1,:),Wr{i}.x(1,:)]; Y = [Wl{i}.x(2,:),Wr{i}.x(2,:)]; Z = [Wl{i}.x(3,:),Wr{i}.x(3,:)];
    plot3(X,Y,Z,'.')
    X = [Wf{i}.x(1,:),Wb{i}.x(1,:)]; Y = [Wf{i}.x(2,:),Wb{i}.x(2,:)]; Z = [Wf{i}.x(3,:),Wb{i}.x(3,:)];
    plot3(X,Y,Z,'.')
    end
    X = Wu.x(1,:); Y = Wu.x(2,:); Z = Wu.x(3,:);
    plot3(X,Y,Z,'.')
    X = Wd.x(1,:); Y = Wd.x(2,:); Z = Wd.x(3,:);
    plot3(X,Y,Z,'.'), hold off
    axis equal
end
end