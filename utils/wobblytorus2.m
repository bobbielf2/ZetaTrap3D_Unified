function s = wobblytorus2(m,n,a,k)
% Geometry generation following Barnett's idea in BIE3D wobblytorus.
% Generate a wobbly torus by revolving a generating curve about z-axis.
% Steps:
%   1. define a planar curve rho = f(th) in polar form, where:
%       f(th) = 1 + a*cos(m*th), a=wobbly ampl, m=wobbly wave number
%   2. define generating curve [x,z] = [R,0] + r*f(th)*[cos(th),sin(th)], rotate about z-axis
%   3. add flavor by also rotating f(th) based on the major (toroidal) angle p=phi
%       f(th,phi) = 1 + a*cos(m*th + n*phi), n=wave number along toroidal direction
%
% Input:
%   m,n = wavenumbers in the poloidal and toroidal directions, resp.
%   a = wobbly amplitude
%   k = aspect ratio of parametric space (phi,theta), [0,2*k*pi)x[0,2*pi)
% Output:
%   s = struct containing geom info
%   s.Z =  [x(u,v); y(u,v); z(u,v)] func handle for surface parameterization
%   s.Zu = [xu(u,v); yu(u,v); zu(u,v)] func handle for the partials w.r.t. u
%   s.Zv = [xv(u,v); yv(u,v); zv(u,v)] func handle for the partials w.r.t. v
%       All of above assume u and v are row vectors.
% 
% Bowei Wu 4/28/20
% Added function "get_x_xp_xpp_xppp" computing sample points and
% derivatives more efficiently

if nargin == 0 && nargout == 0, test_wobblytorus; return; end
if nargin < 4, k = 1; end
if nargin < 3, a = 0.2; end
if nargin < 2, n = 2; end
if nargin < 1, m = 3; end

% planar curve
ps = randn()*3+pi*1.341*0; % phase shift
ts = randn()*3+pi*0.234*0; % t shift
% ps = 0; % phase shift
% ts = 0; % t shift
f = @(p,t) 1 + a * cos(m*t+n/k*p+ps); %generating curve
fp= @(p,t) -n/k*a*sin(m*t+n/k*p+ps); %partials
ft= @(p,t) -m*a*sin(m*t+n/k*p+ps);
fpp= @(p,t) -(n/k)^2*a*cos(m*t+n/k*p+ps); %2nd partials
ftt= @(p,t) -m^2*a*cos(m*t+n/k*p+ps);
fpt= @(p,t) -n/k*m*a*cos(m*t+n/k*p+ps);

%t = linspace(0,2*pi,100); plot(f(0,t).*cos(t),f(0,t).*sin(t)); axis equal %plot it

% toroidal surface
R = 1; r = 0.5;
x = @(p,t) (R + r * f(p,t) .* cos(t+ts)) .* cos(p/k);
y = @(p,t) (R + r * f(p,t) .* cos(t+ts)) .* sin(p/k);
z = @(p,t) r * f(p,t) .* sin(t+ts);

% derivatives
xp = @(p,t) -R/k * sin(p/k) + r * (fp(p,t) .* cos(p/k) - 1/k*f(p,t) .* sin(p/k)) .* cos(t+ts);
yp = @(p,t)  R/k * cos(p/k) + r * (fp(p,t) .* sin(p/k) + 1/k*f(p,t) .* cos(p/k)) .* cos(t+ts);
zp = @(p,t) r * fp(p,t) .* sin(t+ts);

xt = @(p,t) r * (ft(p,t) .* cos(t+ts) - f(p,t) .* sin(t+ts)) .* cos(p/k);
yt = @(p,t) r * (ft(p,t) .* cos(t+ts) - f(p,t) .* sin(t+ts)) .* sin(p/k);
zt = @(p,t) r * (ft(p,t) .* sin(t+ts) + f(p,t) .* cos(t+ts));

% 2nd derivatives
xpp = @(p,t) -R/k^2 * cos(p/k) + r * (fpp(p,t) .* cos(p/k) - 2/k * fp(p,t) .* sin(p/k) - 1/k^2*f(p,t) .* cos(p/k)) .* cos(t+ts);
ypp = @(p,t) -R/k^2 * sin(p/k) + r * (fpp(p,t) .* sin(p/k) + 2/k * fp(p,t) .* cos(p/k) - 1/k^2*f(p,t) .* sin(p/k)) .* cos(t+ts);
zpp = @(p,t) r * fpp(p,t) .* sin(t+ts);

xtt = @(p,t) r * (ftt(p,t) .* cos(t+ts) - 2 * ft(p,t) .* sin(t+ts) - f(p,t) .* cos(t+ts)) .* cos(p/k);
ytt = @(p,t) r * (ftt(p,t) .* cos(t+ts) - 2 * ft(p,t) .* sin(t+ts) - f(p,t) .* cos(t+ts)) .* sin(p/k);
ztt = @(p,t) r * (ftt(p,t) .* sin(t+ts) + 2 * ft(p,t) .* cos(t+ts) - f(p,t) .* sin(t+ts));

xpt = @(p,t) r * ( (fpt(p,t) .* cos(p/k) - 1/k*ft(p,t) .* sin(p/k)) .* cos(t+ts) - (fp(p,t) .* cos(p/k) - 1/k*f(p,t) .* sin(p/k)) .* sin(t+ts) );
ypt = @(p,t) r * ( (fpt(p,t) .* sin(p/k) + 1/k*ft(p,t) .* cos(p/k)) .* cos(t+ts) - (fp(p,t) .* sin(p/k) + 1/k*f(p,t) .* cos(p/k)) .* sin(t+ts) );
zpt = @(p,t) r * ( fpt(p,t) .* sin(t+ts) + fp(p,t) .* cos(t+ts) );

% output
s.Z  = @(u,v) [x(u,v); y(u,v); z(u,v)];
s.Zu = @(u,v) [xp(u,v); yp(u,v); zp(u,v)];
s.Zv = @(u,v) [xt(u,v); yt(u,v); zt(u,v)];
s.Zuu= @(u,v) [xpp(u,v);ypp(u,v);zpp(u,v)];
s.Zuv= @(u,v) [xpt(u,v);ypt(u,v);zpt(u,v)];
s.Zvv= @(u,v) [xtt(u,v);ytt(u,v);ztt(u,v)];
s.topo = 't';
s.type = 'periodic'; s.L = 0;



function test_wobblytorus

k = 3;
s = wobblytorus2(3,2,0.2,k);

n = 30;
p = n*k+1; q = n+1;
[u,v] = meshgrid(linspace(0,2*k*pi,p),linspace(0,2*pi,q));
u = u(:)'; v = v(:)'; % turn into row vectors

% 3-by-n arrays of points and partials and normals
X = s.Z(u,v);
Xu= s.Zu(u,v);
Xv= s.Zv(u,v);

NX = cross(Xu,Xv);
J = vecnorm(NX); %length of normal (jacobian)
NX = NX./J; % unit normal

% plot it
x = reshape(X(1,:),q,p);
y = reshape(X(2,:),q,p);
z = reshape(X(3,:),q,p);
nx = reshape(NX(1,:),q,p);
ny = reshape(NX(2,:),q,p);
nz = reshape(NX(3,:),q,p);
mesh(x,y,z,'FaceAlpha',0); axis equal; hold on 
quiver3(x,y,z,nx,ny,nz); axis equal; hold off 

