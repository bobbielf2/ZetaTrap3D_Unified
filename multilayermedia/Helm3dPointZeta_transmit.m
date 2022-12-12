function [zwei, indsten, sten] = Helm3dPointZeta_transmit(ka1,ka2,ord,lptype,P,inds)
% Compute zeta correction weights & their associated global column index in
% the matrix for the Helmholtz potentials on a curved-rectangular patch.
% Correction order <= 9 for weakly singular kernels & order <= 7 for
% hypersingular kernel

N = size(P.x,2);
if nargin < 6, inds = 1:N; end

if ord <= 3 && ~strcmp(lptype,'dn') || ord <= 1
    zwei = zeros(1,numel(inds)); indsten = inds; sten = []; return
end

% set up kernel param
if strcmp(lptype,'s')       % SLP kernel 1/(4*pi)*e^(i*k*r)*J/r
    p = 1; q = 0;
elseif strcmp(lptype,'d')   % DLP kernel 1/(4*pi)*e^(i*k*r)*(1-i*k*r)*dot(ny,x-y)*J/r^3
    p = -1; q = 1;
elseif strcmp(lptype,'sn')  % SLPn kernel 1/(4*pi)*e^(i*k*r)*(1-i*k*r)*dot(nx,y-x)*J/r^3
    p = -1; q = 1;
elseif strcmp(lptype,'dn')  % DLPn kernel e^(i*k*r)/(4*pi)*((1-i*k*r)*dot(nx,ny)/r^3+k^2*dot(nx,x-y)*dot(ny,x-y)/r^3-3*(1-i*k*r)*dot(nx,x-y)*dot(ny,x-y)/r^5)*J
    p = [-1,-3];  q = [0,2];
end

% setup stencil & geometric quantities: r & |r|^2-Q
E = P.E(inds); F = P.F(inds); G = P.G(inds);
[d,r2mQh2,sten,indsten] = zeta_setup(ord,P,p,q,inds); % d = x - y = targ - src
r = sqrt(d{1}.^2+d{2}.^2+d{3}.^2);
irr = (1./r).^2;
k1r = ka1*r; % ka_1*|x-y|
k2r = ka2*r; % ka_2*|x-y|
muy = (P.NX{1}(indsten).*d{1}+P.NX{2}(indsten).*d{2}+P.NX{3}(indsten).*d{3}); % dot(ny,x-y)
mux = (P.nx(1,inds).*d{1}+P.nx(2,inds).*d{2}+P.nx(3,inds).*d{3}); % dot(nx,x-y)
J = P.J(indsten); % jacobian
ck1r = cos(k1r);
ck2r = cos(k2r);
csk1r = ck1r+k1r.*sin(k1r);
csk2r = ck2r+k2r.*sin(k2r);
if isfield(P,'type') && strcmp(P.type,'periodic')
    % Quasi-periodic scattering, account for Bloch phase
    J(sten.indu_p) = J(sten.indu_p)*P.alph_x;
    J(sten.indu_m) = J(sten.indu_m)/P.alph_x;
    J(sten.indv_p) = J(sten.indv_p)*P.alph_y;
    J(sten.indv_m) = J(sten.indv_m)/P.alph_y;
end
if strcmp(lptype,'s')       % SLP smooth part 1/(4*pi)*Re[e^(i*k*r)]*J
    ckr_irr = (ck1r-ck2r).*irr;
    ckr_irr(1,:) = (ka2^2-ka1^2)/2;
    c = {1/4/pi * ckr_irr .* J};
elseif strcmp(lptype,'d')   % DLP smooth part 1/(4*pi)*Re[e^(i*k*r)*(1-i*k*r)]*dot(ny,x-y)*J
    cskr_irr = (csk1r - csk2r).*irr;
    c = {1/4/pi * cskr_irr .* J .* muy};
elseif strcmp(lptype,'sn')  % SLPn kernel 1/(4*pi)*r^(-3)*dot(nx,y-x)*J
    cskr_irr = (csk1r - csk2r).*irr;
    c = {-1/4/pi * cskr_irr .* J .* mux};
elseif strcmp(lptype,'dn')  % DLPn kernel 1/(4*pi)*(r^(-3)*dot(nx,ny)-3*r^(-5)*dot(nx,x-y)*dot(ny,x-y))*J
    if ord <=3
        c = {1/4/pi * (ka1^2-ka2^2)/2 * J};
    elseif ord<=9
        cskr_irr = (csk1r - csk2r).*irr;
        cskr_irr(1,:) = (ka1^2-ka2^2)/2;
        kc3cskr_irr = ka1^2*ck1r - ka2^2*ck2r - 3*cskr_irr;
        nynx = P.nx(1,inds).*P.NX{1}(indsten)+P.nx(2,inds).*P.NX{2}(indsten)+P.nx(3,inds).*P.NX{3}(indsten);
    c = {1/4/pi * cskr_irr .* J .* nynx,... % 1/r^3 hypersingular part
        1/4/pi * kc3cskr_irr .* J .* muy .* mux};   % 1/r^3 weakly singular part
    end
end

% precompute singular zeta weights
if ord <= 5  % O(h^5) correction. DLPn need bigger stencil, others need 1-pt.
    if strcmp(lptype,'dn')
        [~,Zsd1,Zsd2,Zsd3,~,Zs2]=epstein_zeta_d4_s2d1(-1,E,F,G);
        pre.s = [-1;1];
        pre.zeta={ [],Zsd1,Zsd2,Zsd3;
                  Zs2,  [],  [],  []}; % pre.zeta{i,j}=D^(j+1)Z(pre.s(i))
    end
elseif ord <= 7 % O(h^7) correction. Need 4th zeta deriv. (DLPn need higher, not implemented!)
    if strcmp(lptype,'dn')
        [Zs,Zsd1,Zsd2,Zsd3,Zsd4,Zsd5,Zsd6]=epstein_zeta_d7(-3,E,F,G);
        [Zs2,Zs2d1,Zs2d2,Zs2d3,~,Zs4]=epstein_zeta_d4_s2d1(-1,E,F,G);
        pre.s = [-3;-1;1];
        pre.zeta={  Zs,  Zsd1, Zsd2, Zsd3, Zsd4,Zsd5,Zsd6;
                   Zs2, Zs2d1,Zs2d2,Zs2d3, [], [],  [];
                   Zs4, [],  [],  [], [],  [], []}; % pre.zeta{i,j}=D^(j+1)Z(pre.s(i))
    else
        [Zs,Zsd1,Zsd2,Zsd3,Zsd4,Zs2,Zs2d1]=epstein_zeta_d4_s2d1(-3,E,F,G);
        pre.s = [-3;-1];
        pre.zeta={ Zs,Zsd1,Zsd2,Zsd3,Zsd4;
                  Zs2,Zs2d1,  [],  [], []}; % pre.zeta{i,j}=D^(j+1)Z(pre.s(i))
    end
elseif ord <= 9 % O(h^9) correction. Need 7th zeta deriv
    if strcmp(lptype,'dn')
        [Zs,Zsd1,Zsd2,Zsd3,Zsd4,Zsd5,Zsd6,Zsd7,Zsd8,Zsd9]=epstein_zeta_d10(-5,E,F,G);
        [Zs2,Zs2d1,Zs2d2,Zs2d3,Zs2d4,Zs2d5,Zs2d6]=epstein_zeta_d7(-3,E,F,G);
        [Zs4,Zs4d1,Zs4d2,Zs4d3,~,Zs6] = epstein_zeta_d4_s2d1(-1,E,F,G);
        pre.s = [-5;-3;-1;1];
        pre.zeta={  Zs,  Zsd1, Zsd2, Zsd3, Zsd4,Zsd5,Zsd6,Zsd7,Zsd8,Zsd9;
                   Zs2, Zs2d1,Zs2d2,Zs2d3,Zs2d4,Zs2d5,Zs2d6,[],[], [];
                   Zs4,Zs4d1,Zs4d2,Zs4d3,[],[], [], [],[], [];
                   Zs6, [],[],[],[],[],[],[],[], []}; % pre.zeta{i,j}=D^(j+1)Z(pre.s(i))
    else
        [Zs,Zsd1,Zsd2,Zsd3,Zsd4,Zsd5,Zsd6,Zsd7]=epstein_zeta_d7(-5,E,F,G);
        [Zs2,Zs2d1,Zs2d2,Zs2d3,Zs2d4,Zs4,Zs4d1]=epstein_zeta_d4_s2d1(-3,E,F,G);
        pre.s = [-5;-3;-1];
        pre.zeta={  Zs,  Zsd1, Zsd2, Zsd3, Zsd4,Zsd5,Zsd6,Zsd7;
                   Zs2, Zs2d1,Zs2d2,Zs2d3,Zs2d4,[], [], [];
                   Zs4, Zs4d1,[],[],[],[],[],[]}; % pre.zeta{i,j}=D^(j+1)Z(pre.s(i))
    end
elseif ord <= 11
    if strcmp(lptype,'dn')
        error('DLPn: correction order >= %d not yet implemented!',ord)
    else
        [Zs,Zsd1,Zsd2,Zsd3,Zsd4,Zsd5,Zsd6,Zsd7,Zsd8,Zsd9,Zsd10]=epstein_zeta_d10(-7,E,F,G);
        [Zs2,Zs2d1,Zs2d2,Zs2d3,Zs2d4,Zs2d5,Zs2d6,Zs2d7]=epstein_zeta_d7(-5,E,F,G);
        [Zs4,Zs4d1,Zs4d2,Zs4d3,Zs4d4,Zs6,Zs6d1] = epstein_zeta_d4_s2d1(-3,E,F,G);
        pre.s = [-7;-5;-3;-1];
        pre.zeta={  Zs,  Zsd1, Zsd2, Zsd3, Zsd4,Zsd5,Zsd6,Zsd7,Zsd8,Zsd9,Zsd10;
                   Zs2, Zs2d1,Zs2d2,Zs2d3,Zs2d4,Zs2d5,Zs2d6,Zs2d7,[],[],[];
                   Zs4, Zs4d1,Zs4d2,Zs4d3,Zs4d4,[],[],[],[],[],[];
                   Zs6,Zs6d1,[],[],[],[],[],[], [],[],[]};
    end
else
    error('correction order >= %d not yet implemented!',ord)
end

% compute geometric zeta weights
if ~exist('pre','var')
    zwei=zeta_wei(ord,sten,p,q,c,r2mQh2,P.h,E,F,G);
else
    zwei=zeta_wei(ord,sten,p,q,c,r2mQh2,P.h,pre);
end

% diagonal limit of the smooth component for SLP and DLPn
if strcmp(lptype,'s')
    zwei(1,:) = zwei(1,:) + 1i*(ka1-ka2)/4/pi * P.w(inds);
elseif strcmp(lptype,'dn')
    zwei(1,:) = zwei(1,:) + 1i*(ka1^3-ka2^3)/12/pi * P.w(inds);
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
    sten.indu_p = indu_p;
    sten.indu_m = indu_m;
    sten.indv_p = indv_p;
    sten.indv_m = indv_m;
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