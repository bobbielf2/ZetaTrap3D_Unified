% This example runs a convergence test on evaluating the Helmholtz
% potentials on a surface patch in 3D using zeta corrected trapezoidal
% rules. ONLY the first fundamental form is needed (instead of requiring
% higher-order parametric derivatives as in the first version). 
% See the manuscript [1], Figure 6, for more details.
%
% [1] Wu, B., & Martinsson, P.G. (2022, arXiv:xxxx). A Unified Trapezoidal
%     Quadrature Method for Singular and Hypersingular Boundary Integral
%     Operators on Curved Surfaces.

% Bowei Wu 11/29/2020; updated 9th order quadrature 7/8/2021;
%          updated hypersingular 11/10/21

addpath utils
addpath zetafunc

ord = 5; % desired order of convergence, can pick ord=1,3,5,7,9 
ka = 5+1i; % wavenumber

% generate random quartic patch
rng(2023);
[ru,rv,patch] = randpatch;
% define density function (with compact support)
rpat = 0.6; % radius of the patch
a = randn(); b = randn();
f = @(u,v) (a*cos(a+u)+b*sin(b+v)).* exp(-40*((u.^2+v.^2)/rpat^2).^6); 


% precompute zeta vals
E=dot(ru,ru); F=dot(ru,rv); G=dot(rv,rv); % compute first fundamental form
pre = zeta_precomp(ord,E,F,G);

% convergence test
hh = 0.5.^([4:.3:8,10]);
val_d = []; val_s = []; val_sn = []; val_dn = [];
for h = hh
    % parametric grid & density func
    n = ceil(rpat/h);
    [u,v] = ndgrid(h*(-n:n)); u = u(:).'; v = v(:).';
    sig = f(u,v).'; % density
    
    % define geometry (quartic surface patch centered at 0)
    [rvec,r_u,r_v] = quarticSurf(u,v,patch);	
	% plot surface
   	if h == hh(2), figure(1); subplot(1,2,1), PlotPatch(n,rvec,sig,ru,rv); end 
    
    % define kernels via punctured trapezoidal rule
    r2 = dot(rvec,rvec);  	% r^2
    r = sqrt(r2);          	% r
    ker_r1 = 1./r*h^2;      % 1/r kernel (for SLP)
    ker_r3 = ker_r1./r2;    % 1/r^3 kernel (for DLP & SLPn & DLPn)
    ker_r5 = ker_r3./r2;    % 1/r^5 kernel (for DLPn)
    ind = find(r2==0);     	% locate singular point
    ker_r1(ind) = 0;       	% punctured trapezoidal rule
    ker_r3(ind) = 0;
    ker_r5(ind) = 0;
    
    % smooth components for each kernel
    nJ = cross(r_u,r_v);    % normal vec r_u x r_v (not normalized)
    J = vecnorm(nJ);        % jacobian
    n0 = cross(ru,rv)/norm(cross(ru,rv));   % unit normal at 0
    ikr = 1i*ka*r;          % i*k*r
    eikr = exp(ikr);        % e^(i*k*r)
    sn_fac = sum(rvec.*n0,1).*J;    % SLPn smooth factor: (r.n0)*jacobian
    d_fac = -dot(rvec,nJ);          % DLP smooth factor: -r.(r_u x r_v)
    dn_fac1 = sum(nJ.*n0,1);        % DLPn smooth factor for n0.(r_u x r_v)/r^3
    dn_fac2 = sum(rvec.*n0,1).*dot(rvec,nJ); % DLPn smooth factor for -3(r.n0)*[r.(r_u x r_v)]/r^5 and ka^2*(r.n0)*[r.(r_u x r_v)]/r
    
    % kernels
    ker_s = eikr.*J.*ker_r1;    % SLP
    ker_sn = eikr.*(1-ikr).*sn_fac.*ker_r3;    % SLPn
    ker_d = eikr.*(1-ikr).*d_fac.*ker_r3;      % DLP
    ker_dn = eikr.*(((1-ikr).*dn_fac1 + ka^2*dn_fac2).*ker_r3 - 3*(1-ikr) .*dn_fac2.*ker_r5); % DLPn kernel
    ker_s(ind) = 1i*ka*h^2*J(ind);             % SLP smooth part diag limit
    ker_dn(ind) = 1i*ka^3/3*h^2*J(ind);        % DLPn smooth part diag limit
    
    % zeta correction
    Q = E*u.^2+2*F*u.*v+G*v.^2; % 1st fundamental form
    r2mQ = r2-Q;                % smooth factors for product integration
    % SLP correction: cos(ka*r)*J
    p = -1; offset=0; c = cos(ka*r).*J; 
    ker_s = zeta_correct_rp(ker_s,p,ind,h,r2mQ,c,ord,offset,pre);
    % SLPn correction: [cos(ka*r)+ka*r*sin(ka*r)]*(r.n0)*J
    p = -3; offset=1;  c = (cos(ka*r)+ka*r.*sin(ka*r)).*sn_fac;
    ker_sn = zeta_correct_rp(ker_sn,p,ind,h,r2mQ,c,ord,offset,pre);
    % DLP correction: [cos(ka*r)+ka*r*sin(ka*r)]*[-r.(r_u x r_v)]
    c = (cos(ka*r)+ka*r.*sin(ka*r)).*d_fac;
    ker_d = zeta_correct_rp(ker_d,p,ind,h,r2mQ,c,ord,offset,pre);
    if ord <= 7
        % DLPn correction
        p = [-3,-3,-5]; offset = [0,2,2]; 
        c = [(cos(ka*r)+ka*r.*sin(ka*r)).*dn_fac1;
             ka^2*cos(ka*r).*dn_fac2;
             -3*(cos(ka*r)+ka*r.*sin(ka*r)).*dn_fac2];
        ker_dn = zeta_correct_rp(ker_dn,p,ind,h,r2mQ,c,ord,offset,pre);	% DLPn kernel correction
    end
    
    val_s = [val_s; ker_s*sig];     % SLP vals
    val_d = [val_d; ker_d*sig];     % DLP vals
    val_sn = [val_sn; ker_sn*sig];  % SLP normal gradient 
    val_dn = [val_dn; ker_dn*sig];  % DLP normal gradient 
end

% plot err
err_s = abs(val_s(1:end-1)-val_s(end));
err_d = abs(val_d(1:end-1)-val_d(end));
err_sn= abs(val_sn(1:end-1)-val_sn(end));
err_dn= abs(val_dn(1:end-1)-val_dn(end))*(ord<=7);
hh = hh(1:end-1).';
subplot(1,2,2)
loglog(hh,(hh/hh(1)).^ord*err_s(1)*2,'k--','linewidth',1.5), hold on
loglog(hh,err_s,'o-',hh,err_d,'v-',hh,err_sn,'*-',hh,err_dn,'d-','linewidth',1.5);
title(sprintf('Helmholtz $\\kappa=%.2f%+.2fi$',real(ka),imag(ka)),'interpreter','latex')
xlabel('h'), ylabel('error','rotation',90)
legend({['$O(h^',num2str(ord),')$'],'SLP','DLP','SLPn','DLPn'},'interpreter','latex','location','nw')
hold off

function [ru,rv,patch] = randpatch
% generate random quartic patch by sampling a point on a random surface
m = 3; % petal number of generating curve
n = 2; % twist number along toroidal direction
a = 0.2;
s = wobblytorus(m,n,a); % toroidal surface
u = (rand()*2-1)*pi;
v = (rand()*2-1)*pi;
% derivatives at the sampled point (for surface patch construction only)
[r,ru,rv,ruu,ruv,rvv,ruuu,ruuv,ruvv,rvvv,...
    ruuuu,ruuuv,ruuvv,ruvvv,rvvvv] = s.get_x_xp_xpp_xppp(u,v);
patch.r = r;
patch.ru = ru;
patch.rv = rv;
patch.ruu = ruu;
patch.ruv = ruv;
patch.rvv = rvv;
patch.ruuu = ruuu;
patch.ruuv = ruuv;
patch.ruvv = ruvv;
patch.rvvv = rvvv;
patch.ruuuu = ruuuu;
patch.ruuuv = ruuuv;
patch.ruuvv = ruuvv;
patch.ruvvv = ruvvv;
patch.rvvvv = rvvvv;
end

function [rvec,r_u,r_v] = quarticSurf(u,v,patch)
ru = patch.ru;
rv = patch.rv;
ruu = patch.ruu;
ruv = patch.ruv;
rvv = patch.rvv;
ruuu = patch.ruuu;
ruuv = patch.ruuv;
ruvv = patch.ruvv;
rvvv = patch.rvvv;
ruuuu = patch.ruuuu;
ruuuv = patch.ruuuv;
ruuvv = patch.ruuvv;
ruvvv = patch.ruvvv;
rvvvv = patch.rvvvv;
% construct a quartic surface centered at (0,0) from the given derivatives
% Define quartic surface, r
d1 = ru.*u+rv.*v;
d2 = (ruu.*u.^2+2*ruv.*u.*v+rvv.*v.^2)/2;
d3 = (ruuu.*u.^3+3*ruuv.*(u.^2.*v)+3*ruvv.*(u.*v.^2)+rvvv.*v.^3)/6;
d4 = (ruuuu.*u.^4+4*ruuuv.*(u.^3.*v)+6*ruuvv.*(u.^2.*v.^2)+4*ruvvv.*(u.*v.^3)+rvvvv.*v.^4)/24;
rvec = d1 + d2 + d3 + d4;
% Define derivatives of the surface, dr/du, dr/dv
dr2u = ruu.*u+ruv.*v;
dr2v = ruv.*u+rvv.*v;
dr3u = (ruuu.*u.^2+2*ruuv.*(u.*v)+ruvv.*v.^2)/2;
dr3v = (ruuv.*u.^2+2*ruvv.*(u.*v)+rvvv.*v.^2)/2;
dr4u = (ruuuu.*u.^3+3*ruuuv.*(u.^2.*v)+3*ruuvv.*(u.*v.^2)+ruvvv.*v.^3)/6;
dr4v = (ruuuv.*u.^3+3*ruuvv.*(u.^2.*v)+3*ruvvv.*(u.*v.^2)+rvvvv.*v.^3)/6;
r_u = ru+dr2u+dr3u+dr4u;
r_v = rv+dr2v+dr3v+dr4v;
end

function pre = zeta_precomp(ord,E,F,G)
if ord <= 1
    Zs = epstein_zeta(3,E,F,G);
    pre.s = 3;
    pre.zeta={Zs};
elseif ord <= 3
    [Zs,Zsd1,Zsd2,Zsd3,Zsd4,Zs2,Zs2d1]=epstein_zeta_d4_s2d1(1,E,F,G);
    pre.s = [1;3];
    pre.zeta={Zs,Zsd1,Zsd2,Zsd3, Zsd4;
              Zs2,Zs2d1,  [],  [], []};
elseif ord <= 5
    [Zs,Zsd1,Zsd2,Zsd3,Zsd4,Zsd5,Zsd6,Zsd7]=epstein_zeta_d7(-1,E,F,G);
    [Zs2,Zs2d1,Zs2d2,Zs2d3,Zs2d4,Zs4,Zs4d1]=epstein_zeta_d4_s2d1(1,E,F,G);
    pre.s = [-1;1;3];
    pre.zeta={  Zs,  Zsd1, Zsd2, Zsd3, Zsd4,Zsd5,Zsd6,Zsd7;
                Zs2, Zs2d1,Zs2d2,Zs2d3,Zs2d4, [],  [],[];
                Zs4, Zs4d1,  [],  [], [],  [], [],[]};
elseif ord <= 7
    [Zs,Zsd1,Zsd2,Zsd3,Zsd4,Zsd5,Zsd6,Zsd7,Zsd8,Zsd9,Zsd10]=epstein_zeta_d10(-3,E,F,G);
    [Zs2,Zs2d1,Zs2d2,Zs2d3,Zs2d4,Zs2d5,Zs2d6,Zs2d7]=epstein_zeta_d7(-1,E,F,G);
    [Zs4,Zs4d1,Zs4d2,Zs4d3,Zs4d4,Zs6,Zs6d1] = epstein_zeta_d4_s2d1(1,E,F,G);
    pre.s = [-3;-1;1;3];
    pre.zeta={  Zs,  Zsd1, Zsd2, Zsd3, Zsd4,Zsd5,Zsd6,Zsd7,Zsd8,Zsd9,Zsd10;
                Zs2, Zs2d1,Zs2d2,Zs2d3,Zs2d4,Zs2d5,Zs2d6,Zs2d7,[],[],[];
                Zs4, Zs4d1,Zs4d2,Zs4d3,Zs4d4,[],[],[],[],[],[];
                Zs6,Zs6d1,[],[],[],[],[],[], [],[],[]};
elseif ord <= 9
    [Zs,Zsd1,Zsd2,Zsd3,Zsd4,Zsd5,Zsd6,Zsd7,Zsd8,Zsd9,Zsd10]=epstein_zeta_d10(-5,E,F,G);
    [Zs2,Zs2d1,Zs2d2,Zs2d3,Zs2d4,Zs2d5,Zs2d6,Zs2d7]=epstein_zeta_d7(-3,E,F,G);
    [Zs4,Zs4d1,Zs4d2,Zs4d3,Zs4d4,Zs6,Zs6d1] = epstein_zeta_d4_s2d1(-1,E,F,G);
    pre.s = [-5;-3;-1;1];
    pre.zeta={  Zs,  Zsd1, Zsd2, Zsd3, Zsd4,Zsd5,Zsd6,Zsd7,Zsd8,Zsd9,Zsd10;
                Zs2, Zs2d1,Zs2d2,Zs2d3,Zs2d4,Zs2d5,Zs2d6,Zs2d7,[],[],[];
                Zs4, Zs4d1,Zs4d2,Zs4d3,Zs4d4,[],[],[],[],[],[];
                Zs6,Zs6d1,[],[],[],[],[],[], [],[],[]};
else
    error('not implemented.')
end
end

function ker = zeta_correct_rp(ker,ps,ind,h,r2mQ,cs,ord,qs,pre)
% zeta correction for r^p kernel on surface
% Input:
%   ker = kernel matrix via punctured trapezoidal rule
%   ps  = [p1,p2,... ] powers in r^p for p = p1, p2, ...
%   ind = index location of the singular point
%   h = mesh spacing
%   r2mQ = r^2-Q
%   c = [c1; c2; ...] smooth components for the r^(p_i) kernel correction
%   ord = desired order of convergence
%   E,F,G = first fundamental form coeffs
%   qs = [q1,q2,...] stencil offset parameters, assume an O(h^(2*q_i)) extra
%            smoothness for the non-singular part of the kernel r^(p_i)
if all(mod(ps,2)==0 & ps >= 0), warning('p=%d,r^p is smooth'); return; end
siz = sqrt(numel(ker))*[1,1];   % mesh size
[sub1,sub2] = ind2sub(siz,ind); % subscripts of singular point
np = numel(ps);
for i = 1:np
    p = ps(i); q = qs(i);
    % calculate M, such that 2M+1 terms in the kernel expansion needs correction
    M = ceil((ord-p)/2)-2-q;
    if M < 0, continue; end
    for m = 0:2*M               % correct 2M+1 terms in the kernel expansion
        l1=ceil(3*m/2)+q;	% stencil inner layer = l1
        l2=M+m+q;        	% stencil outer layer = l2+1
        Qpow = m-p/2;           % power of the quadratic form in the m-th term
        fac = binom(p/2,m)*h^(2+p)*cs(i,:).*(r2mQ/h^2).^m; 	% compute smooth factors
        tau = zeta_weights(l1,l2,Qpow,pre);       % compute correction weights
        ker = zeta_correct(ker,sub1,sub2,l1,l2,fac,tau); % apply correction
    end
end
end

function ker = zeta_correct(ker,sub1,sub2,l1,l2,fac,tau)
% Generic zeta correction for a kernel matrix

nt = numel(tau);              % correction weights DoF
siz = sqrt(numel(ker))*[1,1]; % mesh size
% 1. inner layer (symmetric about axes)
u=l1:-1:0; v=0:l1;
taus=tau(1:l1+1);
if l1>1, u=[u,-l1+1:-1]; v=[v,1:l1-1]; taus=[taus,tau(2:l1)]; end % off-axis copies
if l1>0, u=[u,-u]; v=[v,-v]; taus = kron([1,1],taus); end % central symmetry copies
inds = sub2ind(siz,sub1+u,sub2+v);      % kernel indices
ker(inds) = ker(inds)+taus.*fac(inds); 
% 2. intermediate layers (symmetric about center)
if l2 > l1
    u=[]; v=[]; % assemble stencil layers l, l1<l<=l2
    for l=l1+1:l2, u=[u,(l:-1:-l+1)]; v=[v,0:l-1,l:-1:1]; end
    inds=sub2ind(siz,sub1+[u,-u],sub2+[v,-v]);	% kernel indices
    taus = kron([1,1],tau(l1+2:nt-l2));         % tau copies
    ker(inds) = ker(inds) + taus.*fac(inds);
end
% 3. outer layer (anti-symmetric about axes)
if l2 > 0
    u=[l2:-1:1,-l2:-1]'; v=[1:l2,1:l2]';
    inds=sub2ind(siz,sub1+[u;-u],sub2+[v;-v]);	% kernel indices
    taus=kron([1,-1,1,-1],tau(nt-l2+1:nt));     % tau copies
    ker(inds) = ker(inds) + taus.*fac(inds);
end
end

function tau = zeta_weights(l1,l2,Qpow,pre)
% Compute zeta correction weights

C = stencil_inv(l1,l2); % use precomputed inverse of A
if isempty(C)
    % form monomial powers (a,b) for each row u^a * v^b
    % such that a+b=2*l, for all l1<=l<=l2
    a = []; b = [];
    for p=2*(l1:l2), a=[a;(p:-1:0)']; b=[b;(0:p)']; end
    % Construct fitting matrix A on the (u,v) reduced stencil
    n = numel(a); A = zeros(n);
    % 1. inner layer (symmetric about axes)
    u=l1:-1:0; v=0:l1;
    A(:,1:l1+1)=u.^a.*v.^b;
    if l1>0, A(:,1:l1+1)=2*A(:,1:l1+1); end % on-axis pts
    if l1>1, A(:,2:l1)=A(:,2:l1)+2*(-u(2:l1)).^a.*v(2:l1).^b; end % off-axis pts
    % 2. intermediate layers (symmetric about center)
    if l2 > l1
        u=[]; v=[]; % assemble stencil layers l, l1<l<=l2
        for l=l1+1:l2, u=[u,l:-1:-l+1]; v=[v,0:l-1,l:-1:1]; end
        A(:,l1+2:n-l2)=2*u.^a.*v.^b;
    end
    % 3. outer layer (anti-symmetric about axes)
    if l2>0, u=l2:-1:1; v=1:l2; A(:,n-l2+1:n)=2*(u.^a.*v.^b-(-u).^a.*v.^b); end
end
% Compute Wigner limits
W = [];
for l = l1:l2
    s = Qpow-l;
    fallpow = 1/prod(-s-l+1:-s); % Knuth falling factorial
    if 1 % switch off to collect the pairs [2*s,l] for zeta precomputation
        Z = pre.zeta{2*s == pre.s,l+1}; % use precomputed zeta vals
        W = [W;-fallpow*Z];
    else
        disp([2*s,l]) % collect input combinations to accelerate epstein zeta
        W = zeros(size(A,2),1);
    end
end
if ~isempty(C)
    tau = (C*W).';
else
    tau = double(sym(A)\W).'; % solve for zeta weights (row vector)
end
end

function c = binom(n,k)
% generalized binomial coeff, n = any real number, k = positive integer
c = prod(n-k+1:n)/gamma(k+1); 
end

function PlotPatch(n,rvec,sig,ru,rv)
scatter3(0,0,0,100,'k','filled');
hold on;
%quiver3(0,0,0,ru(1),ru(2),ru(3),0.5,'linewidth',2)  % r_u vector
%quiver3(0,0,0,rv(1),rv(2),rv(3),0.5,'linewidth',2)  % r_v vector
surf(reshape(rvec(1,:),2*n+1,2*n+1),reshape(rvec(2,:),2*n+1,2*n+1),reshape(rvec(3,:),2*n+1,2*n+1),reshape(sig,2*n+1,2*n+1))
hold off
axis equal
colorbar
title('random quartic patch','interpreter','latex')
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
zlabel('$z$','interpreter','latex')
%legend({'target','$\mathbf{r}_u$','$\mathbf{r}_v$'},'interpreter','latex')
end
