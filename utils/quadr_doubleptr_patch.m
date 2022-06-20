function s = quadr_doubleptr_patch(s, Ns)
% construct doubly periodic trapezoidal rule on a torus-like surface
% This function is from Alex Barnett's BIE3D package

if nargin == 0, test_quadr_doubleptr; return; end

Nu = Ns(1); Nv = Ns(2);
s.Nu=Nu; s.Nv=Nv; s.N = Nu*Nv;
s.hu=2*pi/Nu; s.hv=2*pi/Nv;
s.h=s.hu; % assume hu=hv for now
% u = (0:Nu-1)'*s.hu;
% v = (0:Nv-1)'*s.hv;
u = (-Nu/2+1:Nu/2)'*s.hu;
v = (-Nv/2+1:Nv/2)'*s.hv;
[u,v] = ndgrid(u,v); % use 'ndgrid' instead of 'meshgrid' to avoid annoying dims swapping
u = u(:)'; v = v(:)'; % turn into row vectors
s.u = u; s.v = v;

% 3-by-n arrays of points and partials and normals
s.x = s.Z(u,v);
s.xu = s.Zu(u,v);
s.xv = s.Zv(u,v);
s.nxj = cross(s.xu,s.xv); % unit normal times jacobian
s.sp = vecnorm(s.nxj);  % length of normal (jacobian), or "speeds"
s.nx = s.nxj ./ s.sp;   % unit normal
s.w  = s.hu*s.hv*s.sp;  % quadr weights
s.E = dot(s.xu,s.xu); 	% first fundamental form E,F,G
s.F = dot(s.xu,s.xv);
s.G = dot(s.xv,s.xv);
% Nu*Nv mesh version
s.X = {reshape(s.x(1,:),[Nu,Nv]),reshape(s.x(2,:),[Nu,Nv]),reshape(s.x(3,:),[Nu,Nv])};
s.NX = {reshape(s.nx(1,:),[Nu,Nv]),reshape(s.nx(2,:),[Nu,Nv]),reshape(s.nx(3,:),[Nu,Nv])};
s.J = reshape(s.sp,[Nu,Nv]);
%s.supp = true(Nu,Nv);   % paremeterization support
%s.type = 'periodic';
if ~isfield(s,'L'), s.L = 0; end

function test_quadr_doubleptr

% Define & show torus-like surface
disp('torus-like surface')
m = 3; % petal number of generating curve
n = 5; % twist number along toroidal direction
a = 0.3;
s = wobblytorus(m,n,a);
Ns = [2,1]*50;
s = quadr_doubleptr_patch(s, Ns);
showsurf(s)

disp('Gauss'' law flux convergence...')
ns = 10;
s_th = rand(1,ns)*2*pi; s_z = a/2 * (rand(1,ns)*2-1); % random src pts inside torus about unit circle (cylindrical coord)
s_str = rand(1,ns)*2 - 1; % random src strengths
err = [];
for nn = 10:5:80
    Nu = 2*nn; Nv = nn;
    s = quadr_doubleptr_patch(s, [Nu,Nv]);
    flux = 0;
    for j = 1:ns
        z = [cos(s_th(j)); sin(s_th(j)); s_z(j)]; % src location
        hold on; plot3(z(1),z(2),z(3),'k.','markersize',20); 
        d = bsxfun(@minus,s.x,z); r = sqrt(sum(d.^2,1));
        ddotn = sum(d.*s.nx,1);
        flux = flux + s_str(j)*sum(s.w.*ddotn./r.^3)/(4*pi);      % surf flux of monopole at z
    end
    fprintf('N=[%d,%d], N=%d:  \terr = %.3g\n',Nu,Nv,s.N,flux - sum(s_str));
    err = [err; abs(flux - sum(s_str))];
end
hold off;
keyboard
