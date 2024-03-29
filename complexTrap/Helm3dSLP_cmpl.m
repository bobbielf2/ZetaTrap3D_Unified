function [A, An] = Helm3dSLP_cmpl(k,t,s,a)
% HELM3DSLP_CMPL.  dense Helmholtz SLP Nystrom eval matrix from sources to targets
% 
% [A, An] = Helm3dSLP_cmpl(k,t,s,a)
% 
% Inputs:
%  k - wavenumber, satisfies imag(k)>0
%  s - source surface struct with fields:
%      x (3*N) nodes, w (1*N) quadr weights
%  t - target struct with fields:
%      x (3*M) nodes, and, if An requested, nx (3*M) normals
%  a - (optional) (3*1) vector that shifts s.x -> s.x + a
%  if_jac - whether or not to include the jacobian (speed) s.sp
% Outputs:
%  A - (M*N) matrix getting potentials from density values
%  An - (M*N) matrix getting t.nx directional derivs from density values
%
% Bowei Wu 3/13/20; updated 2/5/21

if nargin>=4 && ~isempty(a), s.x = s.x + a; end % shift source pts
if size(s.x,2)==size(t.x,2) && max(abs(s.x-t.x),[],'all')<1e-14, self=1; else, self=0; end

d1 = bsxfun(@minus,t.x(1,:).',s.x(1,:)); % 3 coords of displacement matrix (M*N)
d2 = bsxfun(@minus,t.x(2,:).',s.x(2,:));
d3 = bsxfun(@minus,t.x(3,:).',s.x(3,:));
rr = d1.^2+d2.^2+d3.^2;   % dist^2 mat
r  = sqrt(rr);
A = bsxfun(@times, exp(1i*k*r)./r, s.w*(1/4/pi));	% including src quadr wei
if self
    A(diagind(A)) = 0;
end
if nargout>1                  % targ deriv wanted
  ddottn = bsxfun(@times,d1,t.nx(1,:)')+bsxfun(@times,d2,t.nx(2,:)')+bsxfun(@times,d3,t.nx(3,:)');
  An = bsxfun(@times, (1-1i*k*r).*exp(1i*k*r).*ddottn./(r.*rr), s.w*(-1/4/pi));	% monopole deriv, incl src quad wei
  if self
      An(diagind(An)) = 0;
  end
end

function i = diagind(A)
% DIAGIND  Return indices of diagonal of square matrix
%
% Example usage:   A = randn(3,3); A(diagind(A)) = 0;
N = size(A,1); i = sub2ind([N,N], 1:N, 1:N); i = i(:);
