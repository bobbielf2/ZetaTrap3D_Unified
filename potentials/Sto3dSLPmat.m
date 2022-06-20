function [A, P, T] = Sto3dSLPmat(mu,t,s)
% Stokes single-layer surface potential, pressure, and traction matrices,
% assumming source density f=[fx;fy;fz] in node-fast,component-slow format
% Input:
%   mu : viscosity
%   t : target struct with fields:
%       x (3*M) nodes, and, if traction T requested, nx (3*M) normals
%   s : source surface struct with fields:
%       x (3*N) nodes, nx (3*N) unit normals, w (1*N) quadr weights
% Outputs:
%  A : (3M*3N) single-layer potential matrix, [ux;uy;uz]=A*[fx;fy;fz]
%  P : (M*3N) single-layer pressure matrix, p=P*[fx;fy;fz]
%  T : (3M*3N) single-layer traction matrix, [tx;ty;tz]=T*[fx;fy;fz]

if nargin == 0, testcode; return; end
if isempty(mu), mu = 1; end % default viscosity
if size(s.x,2)==size(t.x,2) && max(vecnorm(s.x - t.x))<1e-14, self=1; else, self=0; end

d1 = bsxfun(@minus,t.x(1,:)',s.x(1,:)); % 3 coords of displacement matrix (M*N)
d2 = bsxfun(@minus,t.x(2,:)',s.x(2,:));
d3 = bsxfun(@minus,t.x(3,:)',s.x(3,:));
r2 = d1.^2+d2.^2+d3.^2;   % dist^2 mat
ir = s.w ./ sqrt(r2);   % 1/r, incl wei
ir3 = ir ./ r2;        	% 1/r^3, incl wei
if self, ir(diagind(ir)) = 0; ir3(diagind(ir3)) = 0; end
d1d2 = d1.*d2.*ir3;
d1d3 = d1.*d3.*ir3;
d2d3 = d2.*d3.*ir3;
% tensor product (d x d)/r^3, incl wei
dxd_ir3 = [d1.^2 .*ir3,        d1d2,        d1d3; 
                  d1d2, d2.^2 .*ir3,        d2d3;
                  d1d3,        d2d3, d3.^2 .*ir3];
A = kron((1/8/pi/mu)*eye(3),ir) + (1/8/pi/mu)*dxd_ir3; % SLP, incl prefac
if nargout > 1   % pressure
    ir3 = (1/4/pi) * ir3;     % 1/(4*pi)/r^3, incl wei
    P = [d1.*ir3,d2.*ir3,d3.*ir3];  % SLP pressure
    if nargout > 2 % traction
        ddnx_ir2 = (d1.*t.nx(1,:)'+d2.*t.nx(2,:)'+d3.*t.nx(3,:)')./r2; % (d.nx)/r^2
        if self, ddnx_ir2(diagind(ddnx_ir2)) = 0; end
        T = kron((-3/4/pi)*ones(3),ddnx_ir2).*dxd_ir3; % SLP traction, incl prefac
    end
end

function testcode   % unittest
mu = 0.7;   % viscosity
t = wobblytorus; nu = 50;nv = nu; 	% define target geometry
t = quadr_doubleptr_patch(t, [nu,nv]);  	% discretize targ surf
xin = [0.9;-0.2;0.1]; xout = [1.9;0.7;1.0];	% pts both "far" from surf
so.x = xout; so.nx = randn(3,1); so.w = 1;  % outside src, dummy normal & wei
%si.x = xin; si.nx = randn(3,1); si.w = 1;  % inside src, dummy normal & wei
[A, P, T] = Sto3dSLPmat(mu,t,so);
tau = randn(3,1);                   % rand density vector
uu = A*tau;
pp = P*tau;
tt = T*tau;
uu = reshape(uu,[nu*nv,3])';
flux = dot(dot(uu,t.nx),t.w); % net outflux
disp(['out flux =', num2str(flux)]) % should be zero w/ xout as src (since incompressible)
%showsurffunc(t, pp);
keyboard

function i = diagind(A)
% DIAGIND  Return indices of diagonal of square matrix
%
% Example usage:   A = randn(3,3); A(diagind(A)) = 0;
N = size(A,1); i = sub2ind([N,N], 1:N, 1:N); i = i(:);