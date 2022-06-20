% This example runs a convergence test on solving Helmholtz BVPs exterior 
% to a toroidal surface using zeta corrected trapezoidal rules.
% See the manuscript [1], Figure 8, for more details.
%
% [1] Wu, B., & Martinsson, P.G. (2022, arXiv:xxxx). A Unified Trapezoidal
%     Quadrature Method for Singular and Hypersingular Boundary Integral
%     Operators on Curved Surfaces.

addpath potentials
addpath utils
addpath zetafunc

ord = 5;    % desired order of convergence, can pick ord = 1,3,5,7
ka = 2+1i;	% wavenumber
hyper = 0;  % 1=hypersingular CFIE for Neumann, 0=regularized CFIE for Neumann

% define toroidal surface
rng(2021);
m = 3; % petal number of generating curve
n = 2; % twist number along toroidal direction
a = 0.2;
k = 2; % aspect ratio of (u,v) space
s = wobblytorus2(m,n,a,k);
% plot
s = quadr_doubleptr_patch2(s, [20*k,20]);
figure(1); subplot(1,3,1); showsurf(s); hold off

% Def exact solution
zs = [0.3,-0.9,0; 0.5,0.85,0; -0.99,-0.1,0].';    % source & strength inside torus
strength = ones(3,1)*2.*randn(3,1);
th = rand(1,20)*2*pi; ph = (rand(1,20)-0.5)*pi/4; % test points outside torus
zt = 2.5*[cos(th).*cos(ph);sin(th).*cos(ph);sin(ph)]; 
hold on; plot3(zs(1,:),zs(2,:),zs(3,:),'k.','markersize',20); 
plot3(zt(1,:),zt(2,:),zt(3,:),'r.','markersize',20); hold off;
uexac = Helm3dSLPmat(ka,struct('x',zt),struct('x',zs,'w',1))*strength;

% convergence test (WARNING: form full matrix, slow for large N!)
NN = 16*(1:4); 
errN = []; errD = [];
j = 1;
for Nv = NN     % num of nodes in the v-direction
    % set up quadr & precompute local correction weights
    s = quadr_doubleptr_patch2(s, [k*Nv,Nv]);
    if hyper
        lptypes = {'s','d','sn','dn'};
        ZZ = Helm3dPatchZetaSparse_multi(ka,ord,lptypes,s);
        Zs = ZZ{1}; Zd = ZZ{2}; Zsn = ZZ{3}; Zdn = ZZ{4};
    else
        lptypes = {'s','d','sn'};
        ZZ = Helm3dPatchZetaSparse_multi(ka,ord,lptypes,s);
        Zs = ZZ{1}; Zd = ZZ{2}; Zsn = ZZ{3};
    end
    
    % Diri & Neu data on grid
    [f, g] = Helm3dSLPmat(ka,s,struct('x',zs,'w',1));
    f = f*strength; g = g*strength;
    
    ieta = 1i*abs(ka);
    % solve ext Diri BVP. Ansatz: u = (D - i*eta*S)*tau
    AD = Helm3dDLPmat(ka,s,s) - ieta * Helm3dSLPmat(ka,s,s);
    AD = AD + (Zd - ieta*Zs) + 0.5*eye(size(AD));
    tauD = AD\f;
    uD = (Helm3dDLPmat(ka,struct('x',zt),s) - ieta * Helm3dSLPmat(ka,struct('x',zt),s))*tauD;
    err = max(abs(uD - uexac));
    fprintf('Ns=[%d,%d], N=%d:  \tDiri err = %.3g\n',s.Nu,s.Nv,s.N,err);
    errD = [errD;err];
    
    % solve ext Neu BVP.
    if hyper    % Hypersingular ansatz: u = (S + 1i*eta*D)*tau
        [~,ANs] = Helm3dSLPmat(ka,s,s); [~,ANd] = Helm3dDLPmat(ka,s,s);
        AN = ANs + ieta*ANd; AN = AN + (Zsn + ieta*Zdn) - 0.5*eye(size(AN));
        tauN = AN\g;
        uN = (Helm3dSLPmat(ka,struct('x',zt),s) + ieta * Helm3dDLPmat(ka,struct('x',zt),s))*tauN;
    else        % Regularized ansatz: u = (D*S - i*eta*S)*tau
        [As,Asn] = Helm3dSLPmat(ka,s,s); As = As + Zs; Asn = Asn + Zsn;
        AN = Asn*Asn-ieta*Asn+(0.5*ieta-0.25)*eye(size(Asn));     % Combined Field BIO: A=(0.5*i*eta-0.25)+D'^2-i*eta*D'
        tauN = AN\g;
        uN = Helm3dDLPmat(ka,struct('x',zt),s)*(As*tauN)-ieta*(Helm3dSLPmat(ka,struct('x',zt),s)*tauN);
    end
    err = max(abs(uN - uexac));
    fprintf('Ns=[%d,%d], N=%d:  \tNeu err = %.3g\n',s.Nu,s.Nv,s.N,err);
    errN = [errN;err];
    
    j = j+1;
end

% plot error
Ntot = k*NN.*NN;
err_ref = (NN/NN(end)).^-ord;

subplot(1,3,2)
errDrel = errD/max(abs(uexac));
loglog(Ntot,errDrel,'o-'); hold on
loglog(Ntot,err_ref*errDrel(end),'k--'); hold off
legend({'Dirichlet',['O(h^{',num2str(ord),'})']})
xlabel('N pts'), ylabel('rel. err','rotation',90)
    
subplot(1,3,3)
errNrel = errN/max(abs(uexac));
loglog(Ntot,errNrel,'o-'); hold on
loglog(Ntot,err_ref*errNrel(end),'k--'); hold off
legend({'Neumann',['O(h^{',num2str(ord),'})']})
xlabel('N pts'), ylabel('rel. err','rotation',90)
