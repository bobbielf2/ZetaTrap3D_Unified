function Zs = Helm3dPatchZetaSparse_multi(ka,ord,lptypes,P)
% Compute sparse zeta correction matrix for the Laplace potentials on a
% curved-rectangular patch

np = numel(lptypes);
Zs = cell(1,np);
% compute extrinsic zeta weight
[zweis, indstens] = Helm3dPointZeta_multi(ka,ord,lptypes,P);
% assemble sparse matrix
for i = 1:np
    Zs{i} = patchZetaSparseAssembly(zweis{i},indstens{i});
end