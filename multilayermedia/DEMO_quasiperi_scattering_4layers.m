% Quasi-periodic scattering in 3D
% This code solves a multilayered media transmission problem
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
h = 0.2; % h=0 is flat interface, h>0 is corrugate inteface
ns = 3; % num of interfaces (number of layers = ns + 1)
Nu = 40; Nv = Nu;   % # of pts on periodic surface
N = Nu*Nv;          % total # of pts
Mw = 20;            % # of gauss pts on L, R, F, B walls
Mu = 20;            % # of pts on U wall is Mu^2
P = 30;     % mesh order on proxy sphere ( (P-1)x(2P) mesh )
Rp = 1.5;     % proxy sphere radius
K = 10;     % use 2K+1 terms in Rayleigh-Bloch expansion

% Set up interfaces, separated by distance = 1 apart 
s = cell(ns,1);
s{1} = periSurf2D(Nu,Nv); % scattering surface
for i = 2:ns
    s{i} = s{i-1}; s{i}.x(3,:) = s{i}.x(3,:) - 1; s{i}.X{3} = s{i}.X{3} - 1;
end

% Set up side walls, top & bottom walls, proxy spheres ...
% and plot them
clf, subplot(1,2,1)
[Wl,Wr,Wf,Wb,Wu,Wd,px] = periWall2D(Mw,Mu,P,Rp,s);

%----------------------------------------
% Set scattering params
ka = [5,10,8,12];   % wavenumbers for 4 layers
% ka = [5,10,8,repmat([12,8],1,ceil(ns/2))]; % wave numbers for ns+1 layers

% Incident wave vector = ka(1)*(cos(th)*cos(ph), cos(th)*sin(ph), sin(th))
% where incident angle -pi < th < 0, phase angle -pi/2 <= ph <= pi/2
%th = -acos(1-2*pi/ka(1));	% incident angle (Wood anomaly: th = -acos(1-2*pi*n/s.d/ka) for any integer n)
th = -pi/3;
ph = pi/4;              % phase angle
kv = ka(1)*[cos(th)*cos(ph);cos(th)*sin(ph);sin(th)]; % wave vector
alph_x = exp(1i*kv(1)*s{1}.du);  % Bloch phases (quasi-periodicity parameter)
alph_y = exp(1i*kv(2)*s{1}.dv);

% Transmission BC
f = -exp(1i*sum(kv.*(s{1}.x+[0.5;0.5;0])).'); % f=[ -u^inc; -d/dn{u^inc} ]
f = [f;-exp(1i*sum(kv.*(s{1}.x+[0.5;0.5;0])).').*(1i*sum(kv.*s{1}.nx).')];

%----------------------------------------
% Precompute zeta weights
fprintf('precomputing zeta weights... ')
tic
ord = 5; % zeta quadr order
ZCs = cell(ns,1); ZCd = cell(ns,1);
ZCsn = cell(ns,1); ZCdn = cell(ns,1);
for i = 1:ns
    s{i}.alph_x=alph_x; s{i}.alph_y=alph_y; % allow zeta correction account for Bloch phases
    ZCs{i} = Helm3dPatchZetaSparse_transmit(ka(i),ka(i+1),ord,'s',s{i}); % zeta correction for SLP
    ZCd{i} = Helm3dPatchZetaSparse_transmit(ka(i),ka(i+1),ord,'d',s{i}); % zeta correction for DLP
    ZCsn{i} = Helm3dPatchZetaSparse_transmit(ka(i),ka(i+1),ord,'sn',s{i}); % zeta correction for SLPn
    ZCdn{i} = Helm3dPatchZetaSparse_transmit(ka(i),ka(i+1),ord,'dn',s{i}); % zeta correction for DLPn
end
fprintf('%.2f sec (%d sparse matrices)\n',toc,ns*4)

%----------------------------------------
% Fill system matrix
disp('filling matrices...')
% A: potential on surface (incl contribution from nearest neighbors)
tic
A = num2cell(zeros(ns));
for lu = -1:1
    for lv = -1:1
        dl = [s{1}.du*lu;s{1}.dv*lv;0]; % periodic shift
        for i = 1:ns
            [D1,Dn1] = Helm3dDLPmat(ka(i),  s{i},s{i},dl);
            [S1,Sn1] = Helm3dSLPmat(ka(i),  s{i},s{i},dl);
            [D2,Dn2] = Helm3dDLPmat(ka(i+1),s{i},s{i},dl);
            [S2,Sn2] = Helm3dSLPmat(ka(i+1),s{i},s{i},dl);
            A{i,i} = A{i,i} + alph_x^lu*alph_y^lv*[ D1-D2,   S1-S2;
                                                    Dn1-Dn2, Sn1-Sn2];
        end
        for i = 2:ns
            [D1,Dn1] = Helm3dDLPmat(ka(i),s{i},s{i-1},dl);
            [S1,Sn1] = Helm3dSLPmat(ka(i),s{i},s{i-1},dl);
            A{i,i-1} = A{i,i-1} + alph_x^lu*alph_y^lv*[ D1,  S1;
                                                        Dn1, Sn1];
            [D2,Dn2] = Helm3dDLPmat(ka(i),s{i-1},s{i},dl);
            [S2,Sn2] = Helm3dSLPmat(ka(i),s{i-1},s{i},dl);
            A{i-1,i} = A{i-1,i} + alph_x^lu*alph_y^lv*[ -D2,  -S2;
                                                        -Dn2, -Sn2];
        end
    end
end
for i = 1:ns
    A{i,i} = A{i,i} + kron([1,0;0,-1],speye(N))+[ZCd{i},ZCs{i};ZCdn{i},ZCsn{i}];
end
fprintf('\tA matrix %.2f sec\n',toc)
clear D1 D2 S1 S2 Dn1 Dn2 Sn1 Sn2
clear ZCs ZCsn ZCd ZCdn

% B: proxy to surface
tic
B = num2cell(zeros(ns,ns+1));
for i = 1:ns
    [DS1,DSn1] = Helm3dDSLPmat(ka(i),  s{i},px{i});
    [DS2,DSn2] = Helm3dDSLPmat(ka(i+1),s{i},px{i+1});
    B{i,i}   =  [DS1; DSn1];    
    B{i,i+1} = -[DS2; DSn2];
end
fprintf('\tB matrix %.2f sec\n',toc)
clear DS1 DSn1 DS2 DSn2

% C: surface to left/right & front/back discrepancy
tic
C = num2cell(zeros(ns+1,ns));
for l = -1:1
    for i = 1:ns
        dl_r = [-s{1}.du;l*s{1}.dv;0]; dl_l = [s{1}.du;l*s{1}.dv;0]; % left/right periodic shift
        dl_b = [l*s{1}.du;-s{1}.dv;0]; dl_f = [l*s{1}.du;s{1}.dv;0]; % front/back periodic shift
        % interface i to layer i
        [Dr, Dnr] = Helm3dDLPmat(ka(i),Wr{i},s{i},dl_r);    % left/right descrepancy
        [Sr, Snr] = Helm3dSLPmat(ka(i),Wr{i},s{i},dl_r);
        [Dl, Dnl] = Helm3dDLPmat(ka(i),Wl{i},s{i},dl_l);
        [Sl, Snl] = Helm3dSLPmat(ka(i),Wl{i},s{i},dl_l);
        [Db, Dnb] = Helm3dDLPmat(ka(i),Wb{i},s{i},dl_b);    % front/back descrepancy
        [Sb, Snb] = Helm3dSLPmat(ka(i),Wb{i},s{i},dl_b);
        [Df, Dnf] = Helm3dDLPmat(ka(i),Wf{i},s{i},dl_f);
        [Sf, Snf] = Helm3dSLPmat(ka(i),Wf{i},s{i},dl_f);
        C{i,i} = C{i,i} + [alph_y^l*[alph_x^-1*Dr  - alph_x^2*Dl,  alph_x^-1*Sr  - alph_x^2*Sl;
                                     alph_x^-1*Dnr - alph_x^2*Dnl, alph_x^-1*Snr - alph_x^2*Snl];
                           alph_x^l*[alph_y^-1*Db  - alph_y^2*Df,  alph_y^-1*Sb  - alph_y^2*Sf;
                                     alph_y^-1*Dnb - alph_y^2*Dnf, alph_y^-1*Snb - alph_y^2*Snf]];
        % interface i to layer i+1
        [Dr, Dnr] = Helm3dDLPmat(ka(i+1),Wr{i+1},s{i},dl_r);% left/right descrepancy
        [Sr, Snr] = Helm3dSLPmat(ka(i+1),Wr{i+1},s{i},dl_r);
        [Dl, Dnl] = Helm3dDLPmat(ka(i+1),Wl{i+1},s{i},dl_l);
        [Sl, Snl] = Helm3dSLPmat(ka(i+1),Wl{i+1},s{i},dl_l);
        [Db, Dnb] = Helm3dDLPmat(ka(i+1),Wb{i+1},s{i},dl_b);% front/back descrepancy
        [Sb, Snb] = Helm3dSLPmat(ka(i+1),Wb{i+1},s{i},dl_b);
        [Df, Dnf] = Helm3dDLPmat(ka(i+1),Wf{i+1},s{i},dl_f);
        [Sf, Snf] = Helm3dSLPmat(ka(i+1),Wf{i+1},s{i},dl_f);
        C{i+1,i} = C{i+1,i} + [alph_y^l*[alph_x^-1*Dr  - alph_x^2*Dl,  alph_x^-1*Sr  - alph_x^2*Sl;
                                         alph_x^-1*Dnr - alph_x^2*Dnl, alph_x^-1*Snr - alph_x^2*Snl];
                               alph_x^l*[alph_y^-1*Db  - alph_y^2*Df,  alph_y^-1*Sb  - alph_y^2*Sf;
                                         alph_y^-1*Dnb - alph_y^2*Dnf, alph_y^-1*Snb - alph_y^2*Snf]];
    end
end
fprintf('\tC matrix %.2f sec\n',toc)
clear Dr Dnr Sr Snr Dl Dnl Sl Snl Db Dnb Sb Snb Df Dnf Sf Snf

% Q: proxy to left/right & front/back discrepancies
tic
Q = cell(ns+1,1);
for i = 1:ns+1
    [DSr, DSnr] = Helm3dDSLPmat(ka(i),Wr{i},px{i});
    [DSl, DSnl] = Helm3dDSLPmat(ka(i),Wl{i},px{i});
    [DSb, DSnb] = Helm3dDSLPmat(ka(i),Wb{i},px{i});
    [DSf, DSnf] = Helm3dDSLPmat(ka(i),Wf{i},px{i});
    Q{i} = [DSr  - alph_x*DSl;
            DSnr - alph_x*DSnl;
            DSb  - alph_y*DSf;
            DSnb - alph_y*DSnf];
end
fprintf('\tQ matrix %.2f sec\n',toc)
clear DSr DSnr DSl DSnl DSb DSnb DSf DSnf

% Z: surface to upper & lower wall
tic
ZU = 0; ZD = 0;
for lu = -1:1
    for lv = -1:1
        dl = [s{1}.du*lu;s{1}.dv*lv;0];
        % upper wall
        [Zd, Zdn] = Helm3dDLPmat(ka(1),Wu,s{1},dl);
        [Zs, Zsn] = Helm3dSLPmat(ka(1),Wu,s{1},dl);
        ZU = ZU +alph_x^lu*alph_y^lv*[Zd, Zs; Zdn, Zsn];
        % lower wall
        [Zd, Zdn] = Helm3dDLPmat(ka(ns+1),Wd,s{ns},dl);
        [Zs, Zsn] = Helm3dSLPmat(ka(ns+1),Wd,s{ns},dl);
        ZD = ZD +alph_x^lu*alph_y^lv*[Zd, Zs; Zdn, Zsn];
    end
end
clear Zd Zdn Zs Zsn

% V: proxy to upper wall
[Vds, Vdsn] = Helm3dDSLPmat(ka(1),Wu,px{1});
VU = [Vds; Vdsn];
[Vds, Vdsn] = Helm3dDSLPmat(ka(ns+1),Wd,px{ns+1});
VD = [Vds; Vdsn];
clear Vds Vdsn

% W: Rayleigh-Bloch expansion
[m1,m2] = meshgrid(-K:K); m1=m1(:)'; m2=m2(:)';
km = kv(1:2)+2*pi./[s{1}.du;s{1}.dv].*[m1;m2];  % horizontal wavenumbers
% upper layer
kmperpU = sqrt(ka(1)^2 - sum(km.^2,1));         % vertical wavenumbers
WU = -exp(1i*(km(1,:).*Wu.x(1,:)'+km(2,:).*Wu.x(2,:)'));
WU = [WU; 1i*kmperpU.*WU];
% lower layer
kmperpD = sqrt(ka(ns+1)^2 - sum(km.^2,1));      % vertical wavenumbers
WD = -exp(1i*(km(1,:).*Wd.x(1,:)'+km(2,:).*Wd.x(2,:)'));
WD = [WD; -1i*kmperpD.*WD];
fprintf('\tZ,V,W matrices %.2f sec\n',toc)

%----------------------------------------
% Solve the linear system
disp('solving the system...')
tic
% step 1: rearrange system
B{1,1}     = [B{1,1},    zeros(2*N,(2*K+1)^2)];
B{ns,ns+1} = [B{ns,ns+1},zeros(2*N,(2*K+1)^2)];
C{1,1}     = [C{1,1};    ZU];
C{ns+1,ns} = [C{ns+1,ns};ZD];
Q{1}       = [Q{1},    zeros(4*Mw^2,(2*K+1)^2);
              VU,      WU];
Q{ns+1}    = [Q{ns+1}, zeros(4*Mw^2,(2*K+1)^2);
              VD,      WD];
clear ZU ZD VU VD WU WD
% step 2: reduction to tridiag system
% NOTE: overwrite C by QdagC to save memory
warning('off','MATLAB:rankDeficientMatrix')
% layer 1
C{1,1} = linsolve(Q{1},C{1,1},struct('RECT',true));
A{1,1} = A{1,1} - B{1,1}*C{1,1}; % NOTE: C overwritten by QdagC
% layer 2 to ns
for i = 2:ns
    C{i,i-1} = linsolve(Q{i},C{i,i-1},struct('RECT',true));
    A{i-1,i-1} = A{i-1,i-1} - B{i-1,i}*C{i,i-1};
    A{i,  i-1} = A{i,  i-1} - B{i,  i}*C{i,i-1};
    C{i,i} = linsolve(Q{i},C{i,i},struct('RECT',true));
    A{i-1,i  } = A{i-1,i  } - B{i-1,i}*C{i,i};
    A{i,  i  } = A{i,  i  } - B{i,  i}*C{i,i};
end
% layer ns+1
C{ns+1,ns} = linsolve(Q{ns+1},C{ns+1,ns},struct('RECT',true));
A{ns,ns} = A{ns,ns} - B{ns,ns+1}*C{ns+1,ns};
fprintf('\treduction to tridiag system, %.2f sec\n',toc)
% step 3: solve block-tridiag system (c.f. Golub-Van Loan, 4.5.1)
tic
% blk-LU factorization of A + forward sweep of blk-tridiag solve.
% NOTE1: A_{i,i-1} overwritten by L_{i-1} factor from blk-LU
% NOTE2: A_{i,i} overwritten by U_i factor = A_{i,i}-L_{i-1}*A_{i-1,i-1}
% NOTE3: use same cell arrary for soln (tau) and for RHS to save memory
tau = cell(ns,1); % cell array for RHS (to be overwritten by the soln)
tau{1} = f;
for i = 2:ns
    A{i,i-1} = A{i,i-1}/A{i-1,i-1};     % L_{i-1} factor
    A{i,i} = A{i,i}-A{i,i-1}*A{i-1,i};  % U_i factor
    tau{i} = 0-A{i,i-1}*tau{i-1};       % forward sweep of blk-tridiag solve
end
% backward sweep of blk-tridiag solve. NOTE: RHS cell array overwritten by soln array
tau{ns} = A{ns,ns}\tau{ns}; 
for i = ns-1:-1:1
    tau{i} = A{i,i}\(tau{i}-A{i,i+1}*tau{i+1});
end
% post processing
c = num2cell(zeros(ns+1,1));
for i = 1:ns
    c{i} = c{i}-C{i,i}*tau{i}; % NOTE: C is overwritten by QdagC
    c{i+1} = c{i+1}-C{i+1,i}*tau{i};
end
aU =    c{1}(   px{1}.n+1:px{1}.n+(2*K+1)^2);    c{1}    = c{1}(1:px{1}.n);
aD = c{ns+1}(px{ns+1}.n+1:px{ns+1}.n+(2*K+1)^2); c{ns+1} = c{ns+1}(1:px{ns+1}.n);
fprintf('\tsolve the tridiag system, %.2f sec\n',toc)

% Conservation of flux via Bragg's coeffs (i.e. the "a" vector)
eng_bragg = (kmperpU>0)*(abs(aU).^2 .* kmperpU') + (kmperpD>0)*(abs(aD).^2 .* kmperpD');
eng_ka = -ka(1)*sin(th);
err_flux = abs((eng_bragg - eng_ka)/eng_ka);
fprintf('\trel flux err = %.2e\n',err_flux)

%----------------------------------------
% Evaluate potential field
fprintf('evaluating potential field... ')
% set up plotting mesh
Nt1 = 100; Nt2 = 50*(ns+1);
x = (-Nt1/2:Nt1/2-1)/Nt1*s{1}.du; z = (1:Nt2)/Nt2*(ns+1)-ns;
[xx,zz] = meshgrid(x,z);
% find target pts in each layer
ztmp = cell(ns,1);
ztmp{1} = s{1}.Z(xx(:)',0*xx(:)'); ztmp{1} = reshape(ztmp{1}(3,:),Nt2,Nt1);
for i = 2:ns
    ztmp{i} = ztmp{i-1} - 1;
end
gap = 0.02;
ind = cell(ns+1,1);
ind{1} = zz > ztmp{1}+gap;
for i = 2:ns
    ind{i} = zz < ztmp{i-1}-gap & zz > ztmp{i}+gap; 
end
ind{ns+1} = zz < ztmp{ns}-gap;
t = cell(ns+1,1);
for i = 1:ns+1
    t{i}.x = [xx(ind{i}).*[1,0],zz(ind{i})]';
end
% evaluate field in each layer
tic
u = nan(size(xx))*(1+1i);
for i = 1:ns+1
    u(ind{i}) = Helm3dDSLPmat(ka(i),t{i},px{i})*c{i};  % add proxy correction
end
for lu = -1:1    % add near surface potential contribution
    for lv = -1:1
        dl = [s{1}.du*lu;s{1}.dv*lv;0];
        u(ind{1})=u(ind{1})+alph_x^lu*alph_y^lv*( ...
                  Helm3dDLPmat(ka(1),t{1},s{1},dl)*tau{1}(1:end/2) ...
                + Helm3dSLPmat(ka(1),t{1},s{1},dl)*tau{1}(end/2+1:end));
        for i = 2:ns
            u(ind{i})=u(ind{i})+alph_x^lu*alph_y^lv*( ...
                Helm3dDLPmat(ka(i),t{i},s{i-1},dl)*tau{i-1}(1:end/2) ...
                + Helm3dSLPmat(ka(i),t{i},s{i-1},dl)*tau{i-1}(end/2+1:end) ...
                + Helm3dDLPmat(ka(i),t{i},s{i},dl)*tau{i}(1:end/2) ...
                + Helm3dSLPmat(ka(i),t{i},s{i},dl)*tau{i}(end/2+1:end));
        end
        u(ind{ns+1})=u(ind{ns+1})+alph_x^lu*alph_y^lv*( ...
                  Helm3dDLPmat(ka(ns+1),t{ns+1},s{ns},dl)*tau{ns}(1:end/2) ...
                + Helm3dSLPmat(ka(ns+1),t{ns+1},s{ns},dl)*tau{ns}(end/2+1:end));
    end
end
u(ind{1}) = u(ind{1}) + exp(1i*sum(kv.*(t{1}.x+[0.5;0.5;0])).'); % add incident wave
fprintf('%.2f sec\n',toc)

% Plot solution field (& neighboring fields)
subplot(1,2,2);
np = 1;
uu = kron(alph_x.^(-np:np),u);  % include quasi-periodic neighbors
xp = cell(1,2*np+1);
for i = -np:np, xp{i+np+1} = x+i*s{1}.du; end % periodic x coords
xp = cell2mat(xp);
imagesc(xp,z,real(uu)), axis xy equal tight
colorbar, colormap jet, %caxis([-2,2])
title('quasi-periodic field (real part) in yz-plane')

%-------------------------------------------------------------------------

function [Wl,Wr,Wf,Wb,Wu,Wd,px] = periWall2D(Mw,Mu,P,Rp,s,ifplot)
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
ns = numel(s); % numel of interfaces
%t = gauss(Mw); t= t(:)'; % t in [-1,1]
t = linspace(-1,1,Mw); t = t(:)'; % uniform grid
% build left/right walls (gauss collocation)
Wl = cell(ns+1,1);
L0 = cell(ns+2,1);
L0{1} = 1;                                          % height of top wall
L0{ns+2} = -ns;                                      % height of bottom wall
% height of left bdry of interfaces
L0{2} = s{1}.Z(-s{1}.du/2+0*t,s{1}.dv/2*t); 
L0{2} = L0{2}(3,:);                                 
for i = 3:ns+1
    L0{i} = L0{i-1}-1;
end
[y,z]=meshgrid(t/2*s{1}.dv,(t+1)/2);
for i = 1:ns+1
    zz = z.*(L0{i}-L0{i+1})+L0{i+1};
    Wl{i}.x = [-s{1}.du/2*ones(numel(y),1),y(:),zz(:)]';  % left wall, layer i
    Wl{i}.nx = [1;0;0]+0*Wl{i}.x;
end
Wr = Wl;
for i = 1:ns+1
    Wr{i}.x(1,:) = Wl{i}.x(1,:)+s{1}.du;         	% right wall, layer i
end
% build front/back walls (gauss collocation)
Wf = cell(ns+1,1);
F0 = cell(ns+2,1);
F0{1} = 1;                                          % height of top wall
F0{ns+2} = -ns;                                        % height of bottom wall
% height of front bdry of interfaces
F0{2} = s{1}.Z(s{1}.du/2*t,-s{1}.dv/2+0*t); 
F0{2} = F0{2}(3,:);
for i = 3:ns+1
    F0{i} = F0{i-1}-1;
end
[x,z]=meshgrid(t/2*s{1}.du,(t+1)/2); 
for i = 1:ns+1
    zz = z.*(F0{i}-F0{i+1})+F0{i+1};
    Wf{i}.x = [x(:),-s{1}.dv/2*ones(numel(x),1),zz(:)]';  % front wall, layer i
    Wf{i}.nx = [0;1;0]+0*Wf{i}.x;
end
Wb = Wf;
for i = 1:ns+1
    Wb{i}.x(2,:) = Wf{i}.x(2,:)+s{1}.dv;         	% back wall, layer i
end
% build top wall (uniform collocation)
[x,y]=meshgrid((-Mu/2+0.5:Mu/2-0.5)'/Mu*s{1}.du,(-Mu/2+0.5:Mu/2-0.5)'/Mu*s{1}.dv);
Wu.x =[x(:),y(:),ones(numel(x),1)]';            % top wall
Wu.nx=[0;0;1]+0*Wu.x;
% build bottom wall
Wd = Wu;
Wd.x(3,:) = -ns;
% build proxy pts (spherical grid)
px = cell(ns+1,1);
if isempty(Rp), Rp = s{1}.d*1.5; end
if ischar(P)        % Lebedev points
    [px{1}.nx,w] = lebquad(P);             % proxy normal
    px{1}.x = [0;0;0.5]+Rp*px{1}.nx;      % proxy pts
    px{1}.w = ones(size(px{1}.x(1,:)));	% dummy wei
    %px{1}.w = w;                       % lebedev wei
    px{1}.n = numel(px{1}.w);
elseif isfloat(P)  	% Gauss-Fourier mesh
    u = (0:2*P-1)/P*pi;
    %[v,w] = gauss(P); v = v(:);
    %X = cos(u).*sqrt(1-v.^2); Y = sin(u).*sqrt(1-v.^2); Z = ones(size(u)).*v;
    v = (1:P-1)/P*pi; v = v(:);
    X = cos(u).*sin(v); Y = sin(u).*sin(v); Z = ones(size(u)).*cos(v);
    px{1}.nx = [X(:),Y(:),Z(:)]';              % proxy normals
    px{1}.x = [0;0;0.4]+Rp*px{1}.nx;              % proxy pts
    %px.w = repmat(Rp^2*2*pi/(2*P+1)*w(:)',1,2*P+1);	% proxy wei
    px{1}.w = ones(size(px{1}.x(1,:)));       % dummy wei
    px{1}.n = numel(px{1}.w);
end
for i = 2:ns+1
    px{i} = px{i-1}; px{i}.x(3,:) = px{i}.x(3,:) - 1; % proxy for layer i
end
if ifplot % plot quadrature/collocation points
    % plot 2-peri surfaces
    for i = 1:ns
        N = [s{i}.Nu,s{i}.Nv];
        X = reshape(s{i}.x(1,:),N); Y = reshape(s{i}.x(2,:),N); Z = reshape(s{i}.x(3,:),N);
        XX = [X-s{i}.du,X-s{i}.du,X-s{i}.du;X,X,X;X+s{i}.du,X+s{i}.du,X+s{i}.du];
        YY = [Y-s{i}.dv,Y,Y+s{i}.dv;Y-s{i}.dv,Y,Y+s{i}.dv;Y-s{i}.dv,Y,Y+s{i}.dv];
        ZZ = kron(ones(3),Z);
        mesh(XX,YY,ZZ,'edgecolor','k'), hold on
    end
    % plot proxies
    for i = 1:ns+1
        X = px{i}.x(1,:); Y = px{i}.x(2,:); Z = px{i}.x(3,:);
        plot3(X,Y,Z,'.'), hold on
    end
    % plot L,R,F,B,U,D walls
    for i = 1:ns+1
        X = [Wl{i}.x(1,:),Wr{i}.x(1,:)]; Y = [Wl{i}.x(2,:),Wr{i}.x(2,:)]; Z = [Wl{i}.x(3,:),Wr{i}.x(3,:)];
        plot3(X,Y,Z,'.'), hold on
        X = [Wf{i}.x(1,:),Wb{i}.x(1,:)]; Y = [Wf{i}.x(2,:),Wb{i}.x(2,:)]; Z = [Wf{i}.x(3,:),Wb{i}.x(3,:)];
        plot3(X,Y,Z,'.'), hold on
    end
    X = Wu.x(1,:); Y = Wu.x(2,:); Z = Wu.x(3,:);
    plot3(X,Y,Z,'.'), hold on
    X = Wd.x(1,:); Y = Wd.x(2,:); Z = Wd.x(3,:);
    plot3(X,Y,Z,'.'), hold off
    axis equal
end
end