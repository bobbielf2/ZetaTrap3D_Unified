function [Z, zwei, in_supp] = Helm3dPatchZetaSparse_transmit(ka1,ka2,ord,lptype,P)
% Compute sparse zeta correction matrix for the Helmholtz potentials on a
% curved-rectangular patch

if ord <= 1
    N = size(P.x,2); Z = sparse(N,N); return
else
    % compute extrinsic zeta weight
    [zwei, indsten] = Helm3dPointZeta_transmit(ka1,ka2,ord,lptype,P);
    % assemble sparse matrix
    if isfield(P,'supp')
        [Z, in_supp] = patchZetaSparseAssembly(zwei,indsten,P.supp);
    else
        Z = patchZetaSparseAssembly(zwei,indsten);
    end
end