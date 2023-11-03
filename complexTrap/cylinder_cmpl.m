function s = cylinder_cmpl(cylax,r,L,ncirc,a,b,L0)
% cylinder
% cylax = cylinder axis (e.g. [0.5,0.5,1])
% r = cylinder radius
% ncirc = #pts on cross section circle


u = (0:(ncirc-1))*2*pi/ncirc;
s.h = u(2)-u(1);
nv = round(L/pi*ncirc);
v = (0:nv)*s.h-L;
s.Nu = numel(u); s.Nv = numel(v);
[U,V] = ndgrid(u,v); 
s.u = U; s.v = V; U = U(:).'; V = V(:).';
[Vc,Vcp,Vcpp] = mypath(V,a,b,L0);

s.x = [r*cos(U); r*sin(U); 0*U]+cylax(:).*V;
s.xu = [-r*sin(U);  r*cos(U);  0*U];
s.xv = cylax(:).*ones(size(V));
s.xuu = [-r*cos(U); -r*sin(U);  0*U];

s.xuv = zeros(size(s.x));
s.xvv = zeros(size(s.x));

% complexify
for ii=1:3
    s.x  = s.x  + 1i*cylax(:).*Vc;
    s.xv = s.xv + 1i*cylax(:).*Vcp;
    s.xvv= s.xvv+ 1i*cylax(:).*Vcpp;
end
% geometric quantities
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

% E = real(E); F = real(F); G = real(G);
% L = real(L); M = real(M); N = real(N);


% generate imaginary part of the complex integration path
function [z,zp,zpp] = mypath(x,a,b,L0)
phi = @(x) x.*erfc(x)-exp(-x.^2)/sqrt(pi);
phid = @(x) erfc(x);
phidd=@(x) -exp(-x.^2)*2/sqrt(pi);

z = (phi(b*(x+L0)) - phi(-b*(x-L0)))*a;
zp = (phid(b*(x+L0)) + phid(-b*(x-L0)))*b*a;
zpp = (phidd(b*(x+L0)) - phidd(-b*(x-L0)))*b^2*a;