function [A, P, T] = Sto3dDLPmat(mu,t,s)
% Stokes double-layer surface potential, pressure, and traction matrices,
% assumming source density f=[fx;fy;fz] in node-fast,component-slow format
% Input:
%   mu : viscosity
%   t : target struct with fields:
%       x (3*M) nodes, and, if traction T requested, nx (3*M) normals
%   s : source surface struct with fields:
%       x (3*N) nodes, nx (3*N) unit normals, w (1*N) quadr weights
% Outputs:
%  A : (3M*3N) double-layer potential matrix, [ux;uy;uz]=A*[fx;fy;fz]
%  P : (M*3N) double-layer pressure matrix, p=P*[fx;fy;fz]
%  T : (3M*3N) double-layer traction matrix, [tx;ty;tz]=T*[fx;fy;fz]

if nargin == 0, testcode; return; end
if isempty(mu), mu = 1; end % default viscosity
if size(s.x,2)==size(t.x,2) && max(vecnorm(s.x - t.x))<1e-14, self=1; else, self=0; end

d1 = bsxfun(@minus,t.x(1,:)',s.x(1,:)); % 3 coords of displacement matrix (M*N)
d2 = bsxfun(@minus,t.x(2,:)',s.x(2,:));
d3 = bsxfun(@minus,t.x(3,:)',s.x(3,:));
d1d2 = d1.*d2;
d1d3 = d1.*d3;
d2d3 = d2.*d3;
% tensor product d x d
dxd = [d1.^2,  d1d2,  d1d3; 
        d1d2, d2.^2,  d2d3;
        d1d3,  d2d3, d3.^2];
ir2 = 1./(d1.^2+d2.^2+d3.^2);     % 1/dist^2 mat
if self, ir2(diagind(ir2))=0; end
ir3 = s.w.*sqrt(ir2).*ir2;	 % 1/r^3 * w, weights included
ir5 = ir3.*ir2;            % 1/r^5 * w, weights included
dnyir5 = (d1.*s.nx(1,:)+d2.*s.nx(2,:)+d3.*s.nx(3,:)).*ir5;  % (d.ny)/r^5 * w
A = kron((3/4/pi)*ones(3),dnyir5).*dxd; % DLP matrix
if nargout > 1   % pressure
    nyir3 = {s.nx(1,:).*ir3;s.nx(2,:).*ir3;s.nx(3,:).*ir3};
    P = kron((3*mu/2/pi)*ones(1,3),dnyir5).*[d1,d2,d3] ...
            - (mu/2/pi)*[nyir3{1}, nyir3{2}, nyir3{3}];  % DLP pressure
    if nargout > 2 % traction
        dnx = d1.*t.nx(1,:)'+d2.*t.nx(2,:)'+d3.*t.nx(3,:)'; % d.nx
        dnxir5 = dnx.*ir5;
        ddnxir5 = {d1.*dnxir5;d2.*dnxir5;d3.*dnxir5}; % (d.nx)/r^5 * w
        ddnyir5 = {d1.*dnyir5;d2.*dnyir5;d3.*dnyir5}; % (d.ny)/r^5 * w
        nxnyir5 = (t.nx(1,:)'.*s.nx(1,:)+t.nx(2,:)'.*s.nx(2,:)+t.nx(3,:)'.*s.nx(3,:)).*ir5; % (nx.ny)/r^5 * w
        dnxdnyir5 = dnx.*dnyir5; % (d.nx)(d.ny)/r^5 * w
        T = kron((3*mu/4/pi)*ones(3),nxnyir5-10*dnxdnyir5.*ir2).*dxd; % (3(nx.ny)/r^5-30(d.nx)(d.ny)/r^7)*(d x d), incl prefac
        T = T + kron((3*mu/4/pi)*eye(3),dnxdnyir5); % + 3(d.nx)(d.ny)/r^5 * w * I
        % tensor prod terms
        TP = cell(3);
        for i = 1:3
          for j = 1:3
            TP{i,j} = 2*t.nx(i,:)'.*nyir3{j} + ...	 % 2*(nx)x(ny)/r^3 * w
                      3*s.nx(i,:).*ddnxir5{j} + ...  % 3*(d.nx)*[(ny)x(d)]/r^5 * w
                      3*ddnyir5{i}.*t.nx(j,:)';      % 3*(d.ny)*[(d)x(nx)]/r^5 * w
          end
        end
        T = T + (mu/4/pi)*cell2mat(TP); % DLP traction, incl prefac
    end
end

function testcode   % unittest
mu = 0.7;   % viscosity
t = wobblytorus; nu = 50;nv = nu; 	% define target geometry
t = quadr_doubleptr_patch(t, [nu,nv]);  	% discretize targ surf
xin = [0.9;-0.2;0.1]; xout = [1.9;0.7;1.0];	% pts both "far" from surf
so.x = xout; so.nx = randn(3,1); so.w = 1;  % outside src, dummy normal & wei
%si.x = xin; si.nx = randn(3,1); si.w = 1;  % inside src, dummy normal & wei
[A, P, T] = Sto3dDLPmat(mu,t,so);
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