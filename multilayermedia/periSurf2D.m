function s = periSurf2D(Nu,Nv,h)
% set up a periodic scattering surface on [-.5,.5]^2
% discretized on a (Nu*Nv) grid
s.Nu = Nu; s.Nv = Nv;
a = -0.5; b = -0.5; % shift param
if nargin < 3, h = 0.2; end % height
s.Z = @(u,v) [u;v;h*sin(2*pi*(u-a)).*cos(2*pi*(v-b))];  % periodic rough surface
s.Zu= @(u,v) [ones(size(u));zeros(size(v));h*2*pi*cos(2*pi*(u-a)).*cos(2*pi*(v-b))];
s.Zv= @(u,v) [zeros(size(u));ones(size(v));-h*2*pi*sin(2*pi*(u-a)).*sin(2*pi*(v-b))];
s.du = 1; s.dv = 1;    	% periods in each direction
s.type = 'periodic'; s.L = 1; % param required to compute zeta weights
s.hu = s.du/Nu; s.hv = s.dv/Nv; s.h = s.hu; % assume spacing h=hu=hv
s.u = (-Nu/2:Nu/2-1)'/Nu*s.du;
s.v = (-Nv/2:Nv/2-1)'/Nv*s.dv;
[U,V] = ndgrid(s.u,s.v); U = U(:)'; V = V(:)';
s.x = s.Z(U,V); s.xu = s.Zu(U,V); s.xv = s.Zv(U,V);
nxj = cross(s.xu,s.xv);
s.sp = vecnorm(nxj);    % speed/jacobian
s.nx = nxj./s.sp;       % normal
s.w = s.du*s.dv/Nu/Nv*s.sp;
s.E = dot(s.xu,s.xu);   % first fundamental form
s.F = dot(s.xu,s.xv);
s.G = dot(s.xv,s.xv);
% setup variables for zeta quad
N = [Nu,Nv];
s.supp = true(N);
s.X = {reshape(s.x(1,:),N),reshape(s.x(2,:),N),reshape(s.x(3,:),N)};
s.NX = {reshape(s.nx(1,:),N),reshape(s.nx(2,:),N),reshape(s.nx(3,:),N)};
s.J = reshape(s.sp,N);
end