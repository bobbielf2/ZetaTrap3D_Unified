function [zweis, indstens, stens] = Helm3dPointZeta_multi_cmpl(ka,ord,lptypes,P,inds)
% Compute zeta correction weights & their associated global column index in
% the matrix for the Helmholtz potentials on a curved-rectangular patch.
% Correction order <= 9 for weakly singular kernels & order <= 7 for
% hypersingular kernel

N = size(P.x,2);
if nargin < 5, inds = 1:N; end

np = numel(lptypes);
zweis = cell(1,np);
r2mQh2s = cell(1,np);
indstens = cell(1,np);
stens = cell(1,np);
ps = cell(1,np);
qs = cell(1,np);
cs = cell(1,np);
E = P.E(inds); F = P.F(inds); G = P.G(inds);
for i = 1:np
    lptype = lptypes{i};
    if (ord <= 1 && ~strcmp(lptype,'dn')) || ord <= -1
        zweis{i} = zeros(1,numel(inds)); indstens{i} = inds; stens{i} = []; continue;
    end

    % set up kernel param
    if strcmp(lptype,'s')       % SLP kernel 1/(4*pi)*e^(i*k*r)*J/r
        p = -1; q = 0;
    elseif strcmp(lptype,'d')   % DLP kernel 1/(4*pi)*e^(i*k*r)*(1-i*k*r)*dot(ny,x-y)*J/r^3
        p = -3; q = 1;
    elseif strcmp(lptype,'sn')  % SLPn kernel 1/(4*pi)*e^(i*k*r)*(1-i*k*r)*dot(nx,y-x)*J/r^3
        p = -3; q = 1;
    elseif strcmp(lptype,'dn')  % DLPn kernel e^(i*k*r)/(4*pi)*((1-i*k*r)*dot(nx,ny)/r^3+k^2*dot(nx,x-y)*dot(ny,x-y)/r^3-3*(1-i*k*r)*dot(nx,x-y)*dot(ny,x-y)/r^5)*J
        p = [-3,-3,-5];  q = [0,2,2];
    end
    ps{i} = p; qs{i} = q;

    % setup stencil & geometric quantities: r & |r|^2-Q
    [d,r2mQh2,sten,indsten] = zeta_setup(ord,P,p,q,inds); % d = x - y = targ - src
    kr = ka*sqrt(d{1}.^2+d{2}.^2+d{3}.^2); % ka*|x-y|
    muy = (P.NX{1}(indsten).*d{1}+P.NX{2}(indsten).*d{2}+P.NX{3}(indsten).*d{3}); % dot(ny,x-y)
    mux = (P.nx(1,inds).*d{1}+P.nx(2,inds).*d{2}+P.nx(3,inds).*d{3}); % dot(nx,x-y)
    J = P.J(indsten); % jacobian
    ckr = cos(kr);
    cskr = ckr+kr.*sin(kr);
    if strcmp(lptype,'s')       % SLP kernel 1/(4*pi)*r^(-1)*J
        cs{i} = {1/4/pi * ckr .* J};
    elseif strcmp(lptype,'d')   % DLP kernel 1/(4*pi)*r^(-3)*dot(ny,x-y)*J
        cs{i} = {1/4/pi * cskr .* J .* muy};
    elseif strcmp(lptype,'sn')  % SLPn kernel 1/(4*pi)*r^(-3)*dot(nx,y-x)*J
        cs{i} = {-1/4/pi * cskr .* J .* mux};
    elseif strcmp(lptype,'dn')  % DLPn kernel 1/(4*pi)*(r^(-3)*dot(nx,ny)-3*r^(-5)*dot(nx,x-y)*dot(ny,x-y))*J
        if ord <=1
            cs{i} = {1/4/pi * J};
        elseif ord<=7
            nynx = P.nx(1,inds).*P.NX{1}(indsten)+P.nx(2,inds).*P.NX{2}(indsten)+P.nx(3,inds).*P.NX{3}(indsten);
        cs{i} = {1/4/pi * cskr .* J .* nynx,... % 1/r^3 hypersingular part
                ka^2/4/pi * ckr .* J .* muy .* mux,...   % 1/r^3 weakly singular part
                -3/4/pi * cskr .* J .* muy .* mux};  % 1/r^5 weakly singular part
        end
    end
    r2mQh2s{i} = r2mQh2;
    indstens{i} = indsten;
    stens{i} = sten;
end

% precompute singular zeta weights
if ord <=1
    if any(strcmp(lptypes,'dn'))
        Zs = epstein_zeta_cmpl(3,E,F,G);
        pre.s = 3;
        pre.zeta={Zs};
    end
elseif ord <= 3  % O(h^3) correction. DLPn need bigger stencil, others need 1-pt.
    if any(strcmp(lptypes,'dn'))
        [Zs,Zsd1,Zsd2,Zsd3,Zsd4,Zs2,Zs2d1]=epstein_zeta_d4_s2d1_cmpl(1,E,F,G);
        pre.s = [1;3];
        pre.zeta={ Zs,Zsd1,Zsd2,Zsd3,Zsd4;
                  Zs2,Zs2d1,[], [],  []}; % pre.zeta{i,j}=D^(j+1)Z(pre.s(i))
    else
        [Zs,Zsd1]=epstein_zeta_d4_s2d1_cmpl(1,E,F,G);
        pre.s = 1;
        pre.zeta={Zs,Zsd1};
    end
elseif ord <= 5 % O(h^5) correction. Need 4th zeta deriv. (DLPn need higher, not implemented!)
    if any(strcmp(lptypes,'dn'))
        [Zs,Zsd1,Zsd2,Zsd3,Zsd4,Zsd5,Zsd6,Zsd7]=epstein_zeta_d7(-1,E,F,G);
        [Zs2,Zs2d1,Zs2d2,Zs2d3,Zs2d4,Zs4,Zs4d1]=epstein_zeta_d4_s2d1_cmpl(1,E,F,G);
        pre.s = [-1;1;3];
        pre.zeta={  Zs,  Zsd1, Zsd2, Zsd3, Zsd4,Zsd5,Zsd6,Zsd7;
                   Zs2, Zs2d1,Zs2d2,Zs2d3,Zs2d4, [],  [],[];
                   Zs4, Zs4d1,  [],  [], [],  [], [],[]};
    else
        [Zs,Zsd1,Zsd2,Zsd3,Zsd4,Zs2,Zs2d1]=epstein_zeta_d4_s2d1_cmpl(-1,E,F,G);
        pre.s = [-1;1];
        pre.zeta={ Zs,Zsd1,Zsd2,Zsd3,Zsd4;
                  Zs2,Zs2d1,  [],  [], []}; % pre.zeta{i,j}=D^(j+1)Z(pre.s(i))
    end
elseif ord <= 7 % O(h^7) correction. Need 7th zeta deriv
    if any(strcmp(lptypes,'dn'))
        [Zs,Zsd1,Zsd2,Zsd3,Zsd4,Zsd5,Zsd6,Zsd7,Zsd8,Zsd9,Zsd10]=epstein_zeta_d10(-3,E,F,G);
        [Zs2,Zs2d1,Zs2d2,Zs2d3,Zs2d4,Zs2d5,Zs2d6,Zs2d7]=epstein_zeta_d7(-1,E,F,G);
        [Zs4,Zs4d1,Zs4d2,Zs4d3,Zs4d4,Zs6,Zs6d1] = epstein_zeta_d4_s2d1(1,E,F,G);
        pre.s = [-3;-1;1;3];
        pre.zeta={  Zs,  Zsd1, Zsd2, Zsd3, Zsd4,Zsd5,Zsd6,Zsd7,Zsd8,Zsd9,Zsd10;
                    Zs2, Zs2d1,Zs2d2,Zs2d3,Zs2d4,Zs2d5,Zs2d6,Zs2d7,[],[],[];
                    Zs4, Zs4d1,Zs4d2,Zs4d3,Zs4d4,[],[],[],[],[],[];
                    Zs6,Zs6d1,[],[],[],[],[],[], [],[],[]};
    else
        [Zs,Zsd1,Zsd2,Zsd3,Zsd4,Zsd5,Zsd6,Zsd7]=epstein_zeta_d7(-3,E,F,G);
        [Zs2,Zs2d1,Zs2d2,Zs2d3,Zs2d4,Zs4,Zs4d1]=epstein_zeta_d4_s2d1(-1,E,F,G);
        pre.s = [-3;-1;1];
        pre.zeta={  Zs,  Zsd1, Zsd2, Zsd3, Zsd4,Zsd5,Zsd6,Zsd7;
                   Zs2, Zs2d1,Zs2d2,Zs2d3,Zs2d4,[], [], [];
                   Zs4, Zs4d1,[],[],[],[],[],[]};
    end
elseif ord <= 9
    if any(strcmp(lptypes,'dn'))
        error('DLPn: correction order >= %d not yet implemented!',ord)
    else
        [Zs,Zsd1,Zsd2,Zsd3,Zsd4,Zsd5,Zsd6,Zsd7,Zsd8,Zsd9,Zsd10]=epstein_zeta_d10(-5,E,F,G);
        [Zs2,Zs2d1,Zs2d2,Zs2d3,Zs2d4,Zs2d5,Zs2d6,Zs2d7]=epstein_zeta_d7(-3,E,F,G);
        [Zs4,Zs4d1,Zs4d2,Zs4d3,Zs4d4,Zs6,Zs6d1] = epstein_zeta_d4_s2d1(-1,E,F,G);
        pre.s = [-5;-3;-1;1];
        pre.zeta={  Zs,  Zsd1, Zsd2, Zsd3, Zsd4,Zsd5,Zsd6,Zsd7,Zsd8,Zsd9,Zsd10;
                   Zs2, Zs2d1,Zs2d2,Zs2d3,Zs2d4,Zs2d5,Zs2d6,Zs2d7,[],[],[];
                   Zs4, Zs4d1,Zs4d2,Zs4d3,Zs4d4,[],[],[],[],[],[];
                   Zs6,Zs6d1,[],[],[],[],[],[], [],[],[]};
    end
else
    error('correction order >= %d not yet implemented!',ord)
end

% compute geometric zeta weights
for i = 1:np
    if isempty(zweis{i})
        zweis{i}=zeta_wei(ord,stens{i},ps{i},qs{i},cs{i},r2mQh2s{i},P.h,pre);
        % diagonal limit of the smooth component for SLP and DLPn
        if strcmp(lptypes{i},'s')
            zweis{i}(1,:) = zweis{i}(1,:) + 1i*ka/4/pi*P.h^2 * P.sp(inds);
        elseif strcmp(lptypes{i},'dn')
            zweis{i}(1,:) = zweis{i}(1,:) + 1i*ka^3/12/pi*P.h^2 * P.sp(inds);
        end
    end
end
end

function [dloc,r2mQh2,sten,indsten] = zeta_setup(ord,P,p,q,inds)
% TODO: handle cancellation errors in r^2-Q?

% setup for zeta correction
% Inputs:
%   ord = desired order of convergence
%   P = struct containing info of a geometric patch
%   p = power in r^p
%   q = smoothness offset, assume kernel be multiplied by a O(h^(2*q)) function
% Outputs:
%   dloc = x-y, displacement b/w targ x & src y that live on x's local stencil
%   r2mQ = r^2-Q on the local stencil assoc. with each point on the patch
%          where r = |x-y| = |dloc| local displacement
%                Q = first fundamental form on local parameter space
%   sten = struct containing stencil info, assumed ndgrid (not meshgrid)
%           ordering
%   indsten = global indices of stencil points assoc. with dloc & r2mQ,
%           assumed ndgrid (not meshgrid) ordering

% determine how many layers l in the diamond-shaped stencil l1 <= l <= l2+1
lmax = max(3*ceil((ord-p)/2)-2*q-6); % 3*M+q, where M = ceil((ord-p)/2)-2-q
lmin = min(q);
%if lmax > 0, lmax = lmax+1; end % extra layer if 3M+q>0
sten = zeta_stencil(lmin,lmax);	% build full stencil
% compute local 1st fund. form
uloc = P.h*sten.i';
vloc = P.h*sten.j';
Qloc = P.E(inds).*uloc.^2+2*P.F(inds).*uloc.*vloc+P.G(inds).*vloc.^2;
% compute local r^2
if ~isfield(P,'supp'), P.supp=true(P.Nu,P.Nv); end
[indu,indv]=find(P.supp); indu = indu(inds); indv = indv(inds);
indu = (indu+sten.i)'; 	% global u indices of local stencils (#stcl pts * #surf pts)
indv = (indv+sten.j)';  % global v indices of local stencils (#stcl pts * #surf pts)
% fix indices outside patch
if isfield(P,'type') && strcmp(P.type,'periodic')   % periodic boundary
    indu_p = indu>P.Nu; % u indices > Nu ("plus" side)
    indu(indu_p) = indu(indu_p) - P.Nu; 
    indu_m = indu<1;    % u indices < 1 ("minus" side)
    indu(indu_m) = indu(indu_m) + P.Nu;
    indv_p = indv>P.Nv; % v indices > Nv ("plus" side)
    indv(indv_p) = indv(indv_p) - P.Nv; 
    indv_m = indv<1;    % v indices < 1 ("minus" side)
    indv(indv_m) = indv(indv_m) + P.Nv;
else                                % vanishing boundary (default)
    indu_out = indu>P.Nu | indu < 1; % u indices outside param domain
    indu(indu_out) = 1;
    indv_out = indv>P.Nv | indv < 1; % v indices outside param domain
    indv(indv_out) = 1;
end
indsten=sub2ind([P.Nu,P.Nv],indu,indv); % global indices of stencils
% compute local r vector and |r|^2-Q
dloc = {P.x(1,inds)-P.X{1}(indsten), ...
        P.x(2,inds)-P.X{2}(indsten), ...
        P.x(3,inds)-P.X{3}(indsten)};  % dloc = x - y = targ - src
if isfield(P,'type') && strcmp(P.type,'periodic') % correct periodic shifts
    dloc{1}(indu_p)=dloc{1}(indu_p)-P.L;    % assume param is (u,v,Z(u,v))
    dloc{1}(indu_m)=dloc{1}(indu_m)+P.L;
    dloc{2}(indv_p)=dloc{2}(indv_p)-P.L;
    dloc{2}(indv_m)=dloc{2}(indv_m)+P.L;
    r2loc = dloc{1}.^2+dloc{2}.^2+dloc{3}.^2;
    r2mQ = r2loc-Qloc;      % |r|^2-Q
else
    r2loc = dloc{1}.^2+dloc{2}.^2+dloc{3}.^2;
    r2mQ = r2loc-Qloc;      % |r|^2-Q
    r2mQ(indu_out)=nan;     % set out-of-domain entries 0
    r2mQ(indv_out)=nan;
end
r2mQh2 = r2mQ/P.h^2;        % (|r|^2-Q)/h^2
end

function sten = zeta_stencil(l1,l2)
% Form stencil layers l1 to l2
% Input:
%   l1 = inner layer
%   l2 = outer layer 
%       (Note: also include off-axis nodes of layer l2+1)
% Output:
%   sten = struct containing stencil info
%   sten.i = stencil i coordinates
%   sten.j = stencil j coordinates
%   sten.ind(k+1) = starting index of layer-k in sten.i and sten.j
%   sten.ind(end)-1 = length of sten.i and sten.j
%   q = smoothness offset, assume kernel be multiplied by a O(h^(2*q)) function
assert(l1>=0,'layer l1 must be non-negative!')
i = []; j = [];
ind = []; % starting indices of each layer
if l2 >= l1
    ind = 1;
    if l1 == 0 
        i = 0; j = 0;
        ind = [ind,2];
        for l = l1+1:l2
            il=l:-1:-l+1; jl=[0:l-1,l:-1:1];
            i = [i,il,-il]; j = [j,jl,-jl];
            ind = [ind,ind(end)+2*numel(il)];
        end
    else
        for l = l1:l2
            il=l:-1:-l+1; jl=[0:l-1,l:-1:1];
            i = [i,il,-il]; j = [j,jl,-jl];
            ind = [ind,ind(end)+2*numel(il)];
        end 
    end
    if l2 > 0
        il=[l2:-1:1,-1:-1:-l2]; jl=[1:l2,l2:-1:1];
        i = [i,il,-il]; j = [j,jl,-jl];
        ind = [ind,ind(end)+2*numel(il)];
    end
end
sten.i = i;
sten.j = j;
sten.ind = ind;
sten.l0 = l1;
% scatter(i,j,'k'), axis equal	% plot stencil
end