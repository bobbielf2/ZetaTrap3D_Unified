function [Z, in_supp] = patchZetaSparseAssembly(zwei,sten_indsten,s_supp)
% Assemble sparse matrix for zeta correction on a compactly supported patch
%
% Inputs:
%   zwei = zeta correction weights
%   sten_indsten = "sten" or "indsten"
%       sten = struct of local stencil info, to be used with s.ui & s.vi (see below)
%       indsten = global indices (rectglr mesh) of stencil pts assoc. with zwei
%   s_supp = "s" or "supp"
%       s = struct with node coordinates given by s.ui & s.vi
%       supp = logicals (characteristic func) indicating the support of the
%          parametric func of the geometry on the rectangular mesh; if the
%          patch is doubly-periodic, supp = all ones.
% Output:
%   Z = N-by-N sparse matrix for zeta correction, where N is the number of
%       points in the support (N = numel(find(supp)

if nargin == 0, testcode; return; end % unit test

N = size(zwei,2);	% total num of pts in support
if isstruct(sten_indsten)   % assembly via local stencil search
    % assume the columns of zwei, s.ui, s.vi are permuted the same way
    sten = sten_indsten;
    s = s_supp;
    [VU, p] = sortrows([s.vi(:),s.ui(:)]);
    i_corr = (1:N)'.*ones(1,size(zwei,1));
    j_corr = sub2ind([s.Nu,s.Nv], mod(VU(:,2) + sten.i-1,s.Nu)+1, ...
                                  mod(VU(:,1) + sten.j-1,s.Nv)+1);
    zwei = zwei(:,p).';
    in_supp = ~isnan(zwei);
    Z = sparse(i_corr(in_supp),j_corr(in_supp),zwei(in_supp),N,N);
    Z(p,p) = Z;
else                        % assembly based on pre-defined node-ordering
    indsten = sten_indsten;
    if nargin < 3       % supp not given, assume supp = entire patch
        j_corr = indsten;
        i_corr = repmat(1:N,size(zwei,1),1);
        in_supp = ~isnan(zwei); % non-zero indices (inside rectangular domain)
        Z = sparse(i_corr(in_supp),j_corr(in_supp),zwei(in_supp),N,N);
    else                % supp is given
        supp = s_supp;
        j_supp = zeros(size(supp)); % points outside support has index 0
        j_supp(supp) = 1:N;         % assign linear indexing to on-support points
        j_corr = j_supp(indsten);   % global col indices for correction
        i_corr = repmat(1:N,size(j_corr,1),1); % global row indices for correction
        % ind_corr = (j_supp(indsten)-1)*N+(1:N); % global linear indices for correction
        in_supp = ~(~j_corr) & ~isnan(zwei); % non-zero indices (inside support)
        Z = sparse(i_corr(in_supp),j_corr(in_supp),zwei(in_supp),N,N);
    end
end

function testcode
Nu = 50;
ord = 5;
s = wobblytorus;
s = quadr_doubleptr_patch(s, [Nu,Nu]);

% precompute weights in original ordering
[Dd,indsten,sten] = Lap3dPatchZeta(ord,'d',s);
% ka = 2*pi; [Dd,indsten,sten] = Helm3dPatchZeta(ka,ord,'d',s);

% assemble oringinal matrix
Z0 = patchZetaSparseAssembly(Dd,indsten);

% permute the nodes
ip = randperm(Nu^2);        % a random permutation of 1:Nu^2
Dd = Dd(:,ip); s.ui = s.ui(ip); s.vi = s.vi(ip);    % permute global coordinates & zeta weights accordingly

% assemble permuted matrix
Z = patchZetaSparseAssembly(Dd,sten,s);

% error
v0 = randn(Nu^2,1); % random vector for testing
v = v0(ip); % permuted vector
err1 = norm(Z0(ip,ip) - Z,Inf);
err2 = norm(Z0(ip,:)*v0 - Z*v,Inf);
fprintf('matrix err = %.2g,\tmat-vec err = %.2g\n',err1,err2)