function s = patch_cmpl(L,n,a,b,L0)
t = linspace(-L,L,n+1)';
Nu = n+1; Nv = n+1;
s.Nu = Nu; s.Nv = Nv;
s.h = 2*L/n;
[c,cp,cpp] = mypath(t,a,b,L0);
[U,V] = ndgrid(t);      
s.u = U; s.v = V; U = U(:).'; V = V(:).';
[X,Y] = ndgrid(t+1i*c); X = X(:).'; Y = Y(:).'; 
[Xu, ~] = ndgrid(1+1i*cp,0*t); Xu = Xu(:).';
[ ~,Yv] = ndgrid(0*t,1+1i*cp); Yv = Yv(:).';
[Xuu,~] = ndgrid(1i*cpp,0*t); Xuu = Xuu(:).';
[~,Yvv] = ndgrid(0*t,1i*cpp); Yvv = Yvv(:).';
if 0 % flat
    s.x   = [ X  ;   Y; 0*U];
    s.xu  = [ Xu ; 0*V; 0*U];
    s.xv  = [ 0*U;  Yv; 0*U];
    s.xuu = [ Xuu; 0*V; 0*U];
    s.xuv = [ 0*U; 0*V; 0*U];
    s.xvv = [ 0*U; Yvv; 0*U];
else % bump
    aa = 0.05; H = -L/5;
    s.x   = [ X  ;   Y; H*exp(-aa*(U.^2+V.^2))];
    s.xu  = [ Xu ; 0*V; -2*H*aa*U.*exp(-aa*(U.^2+V.^2))];
    s.xv  = [ 0*U;  Yv; -2*H*aa*V.*exp(-aa*(U.^2+V.^2))];
    s.xuu = [ Xuu; 0*V; H*(-2*aa+4*aa^2*U.^2).*exp(-aa*(U.^2+V.^2))];
    s.xuv = [ 0*U; 0*V; 4*H*aa^2*U.*V.*exp(-aa*(U.^2+V.^2))];
    s.xvv = [ 0*U; Yvv; H*(-2*aa+4*aa^2*V.^2).*exp(-aa*(U.^2+V.^2))];
end
s.nxj = cross(s.xu,s.xv); % unit normal times jacobian
s.sp = sqrt(sum(s.nxj.^2));  % length of normal (jacobian), or "speeds"
s.nx = s.nxj ./ s.sp;   % unit normal
s.w  = s.h^2*s.sp;  % quadr weights
s.E = sum(s.xu.*s.xu); 	% first fundamental form E,F,G
s.F = sum(s.xu.*s.xv);
s.G = sum(s.xv.*s.xv);
s.e = sum(s.xuu.*s.nx); % second fundamental form e,f,g
s.f = sum(s.xuv.*s.nx);
s.g = sum(s.xvv.*s.nx);
% Nu*Nv mesh version
s.X = {reshape(s.x(1,:),[Nu,Nv]),reshape(s.x(2,:),[Nu,Nv]),reshape(s.x(3,:),[Nu,Nv])};
s.NX = {reshape(s.nx(1,:),[Nu,Nv]),reshape(s.nx(2,:),[Nu,Nv]),reshape(s.nx(3,:),[Nu,Nv])};
s.J = reshape(s.sp,[Nu,Nv]);

% generate imaginary part of the complex integration path
function [z,zp,zpp] = mypath(x,a,b,L0)
phi = @(x) x.*erfc(x)-exp(-x.^2)/sqrt(pi);
phid = @(x) erfc(x);
phidd=@(x) -exp(-x.^2)*2/sqrt(pi);

z = (phi(b*(x+L0)) - phi(-b*(x-L0)))*a;
zp = (phid(b*(x+L0)) + phid(-b*(x-L0)))*b*a;
zpp = (phidd(b*(x+L0)) - phidd(-b*(x-L0)))*b^2*a;