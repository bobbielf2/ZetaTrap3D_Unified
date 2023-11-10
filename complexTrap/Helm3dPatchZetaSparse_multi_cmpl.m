function Zs = Helm3dPatchZetaSparse_multi_cmpl(ka,ord,lptypes,P)
% Compute sparse zeta correction matrix for the Helmholtz potentials on a
% curved-rectangular patch
%
% Input:
%   ka: wavenumber
%   ord: order of correction, can be 3, 5
%   lptypes: cell array specifying the types of layer-potential to be computed
%            e.g. lptypes = {'s','d'} computes the SLP and DLP
%            possible types: 's', 'd', 'sn', 'dn'
%   P: a struct containing info of a rectangular surface patch,
%      fields under the struct P include:
%       h  = mesh size
%       Nu = number of points in u-direction
%       Nv = number of points in u-direction
%       x  = 3-by-(Nu*Nv) array, locations of surface points x(u,v)
%       nx = 3-by-(Nu*Nv) array, normal vectors at x(u,v)
%       sp = 1-by-(Nu*Nv) array, Jacobian or "speed weight"
%       E, F, G: 1-by-(Nu*Nv) arrays, coefficients of first fundamental form
%
% Output:
%   Zs: cell array containing (Nu*Nv)-by-(Nu*Nv) sparse matrices,
%       Zs{i} gives zeta correction for the potential specified by lptypes{i}

np = numel(lptypes);
Zs = cell(1,np);
% compute extrinsic zeta weight
if ~isfield(P,'X')
    Ns = [P.Nu,P.Nv];
    P.X  = {reshape(P.x(1,:),Ns),reshape(P.x(2,:),Ns),reshape(P.x(3,:),Ns)};
    P.NX = {reshape(P.nx(1,:),Ns),reshape(P.nx(2,:),Ns),reshape(P.nx(3,:),Ns)};
    P.J  = reshape(P.sp,[P.Nu,P.Nv]);
end
[zweis, indstens] = Helm3dPointZeta_multi_cmpl(ka,ord,lptypes,P);
% assemble sparse matrix
for i = 1:np
    Zs{i} = patchZetaSparseAssembly_cmpl(zweis{i},indstens{i});
end