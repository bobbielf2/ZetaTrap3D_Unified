function [A, An] = Helm3dDSLPmat(k,t,s,a,eta)
% HELM3DDSLPMAT.  Combined-field Helmholtz Nystrom eval matrix DLP-i*eta*SLP
% from sources to targets.
% 
% [A, An] = Helm3dDSLPmat(k,t,s,a)
% 
% Inputs:
%  k - wavenumber, satisfies imag(k)>0
%  s - source surface struct with fields:
%      x (3*N) nodes, w (1*N) quadr weights
%  t - target struct with fields:
%      x (3*M) nodes, and, if An requested, nx (3*M) normals
%  a - (optional) (3*1) vector that shifts s.x -> s.x + a
%  eta - coupling constant, default eta=abs(k)
% Outputs:
%  A - (M*N) matrix getting potentials from density values
%  An - (M*N) matrix getting t.nx directional derivs from density values
%
% Bowei Wu 2/5/21

if nargin==0, testcode; return; end % unittest
if nargin>=4 && ~isempty(a), s.x = s.x + a; end % shift source pts
if nargin<5 || isempty(eta), eta=abs(k); end

d1 = bsxfun(@minus,t.x(1,:)',s.x(1,:)); % 3 coords of displacement matrix (M*N)
d2 = bsxfun(@minus,t.x(2,:)',s.x(2,:));
d3 = bsxfun(@minus,t.x(3,:)',s.x(3,:));
rr = d1.^2+d2.^2+d3.^2;   % dist^2 mat
r  = sqrt(rr);
ddotsn = bsxfun(@times,d1,s.nx(1,:))+bsxfun(@times,d2,s.nx(2,:))+bsxfun(@times,d3,s.nx(3,:));
one_ikr_ir2 = (1-1i*k*r)./rr;
% A = DLP-i|k|*SLP, including src quadr wei
A = bsxfun(@times, exp(1i*k*r)./r.*(one_ikr_ir2.*ddotsn-1i*eta), s.w*(1/4/pi));
if size(s.x,2) == size(t.x,2) && norm(s.x - t.x,Inf) < 1e-14
    A(diagind(A)) = 0;
end
if nargout>1                  % targ deriv wanted
  ddottn = bsxfun(@times,d1,t.nx(1,:)')+bsxfun(@times,d2,t.nx(2,:)')+bsxfun(@times,d3,t.nx(3,:)');
  tndotsn = bsxfun(@times,s.nx(1,:),t.nx(1,:)')+bsxfun(@times,s.nx(2,:),t.nx(2,:)')+bsxfun(@times,s.nx(3,:),t.nx(3,:)');
  % An = d/dn ( DLP-i|k|*SLP )
  An = exp(1i*k*r)./r.* ((k^2-3*one_ikr_ir2).*ddotsn.*ddottn./rr+one_ikr_ir2.*tndotsn ...
                               +1i*eta*one_ikr_ir2.*ddottn);
  An = bsxfun(@times, An, s.w*(1/4/pi));  % incl src quad wei
  if size(s.x,2) == size(t.x,2) && norm(s.x - t.x) < 1e-14
      An(diagind(An)) = 0;
  end
end

function i = diagind(A)
% DIAGIND  Return indices of diagonal of square matrix
%
% Example usage:   A = randn(3,3); A(diagind(A)) = 0;
N = size(A,1); i = sub2ind([N,N], 1:N, 1:N); i = i(:);

function testcode % unittest
% check consistency
s = wobblytorus; nu = 80;nv = nu;       % define src geometry
s = quadr_doubleptr_patch(s, [nu,nv]); 	% discretize src surf
s.L = 0;
t.x = randn(3,20); % 20 rand targ pts
t.nx = randn(3,20); t.nx = t.nx./vecnorm(t.nx); % rand normal
t.x = t.x./vecnorm(t.x)*1.8; % away from src
%showsurf(s); hold on; plot3(t.x(1,:),t.x(2,:),t.x(3,:),'.','markersize',15); hold off
ka = pi+1i*exp(1);
[ADS, ADSn] = Helm3dDSLPmat(ka,t,s);
[AD, ADn] = Helm3dDLPmat(ka,t,s);
[AS, ASn] = Helm3dSLPmat(ka,t,s);
err_DS = max(abs(ADS - (AD-1i*abs(ka)*AS)),[],'all');
err_DSn = max(abs(ADSn - (ADn-1i*abs(ka)*ASn)),[],'all');
fprintf('eval err of D-i|k|S: %.3g\neval err of Dn-i|k|Sn: %.3g\n',err_DS,err_DSn)
