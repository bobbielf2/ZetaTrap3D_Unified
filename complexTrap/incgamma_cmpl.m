function [gs,gsp1,gsp2,gsp3,gsp4] = incgamma_cmpl(s,x)
% My customed fast implementation of the (scaled) upper incomplet gamma 
% function for all real s, including s<=0. (Still slow for complex s.)
%
% Input:
%   (s,x) are the two arguments of the incomplete gamma function defined as
%       incgamma(s,x) = \int_1^\infty t^{s-1}*exp(-x*t) dt
%                     = x^{-s} \int_x^\infty t^{s-1}*exp(-t) dt
%   x>0 and at least one of s and x must be a scalar.
% Output:
%   gs   == igamma(s,x)*x^(-s)
%   gsp1 == igamma(s+1,x)*x^(-(s+1))
%   gsp2 == igamma(s+2,x)*x^(-(s+2))
%   gsp3 == igamma(s+3,x)*x^(-(s+3))
%   gsp4 == igamma(s+3,x)*x^(-(s+4))
%
% Note: this implementation is Fast in the following sense:
% When isreal(s), we use gammainc(x,s,'upper') or expint(x) (FAST) 
% instead of igamma(s,x) (SLOW); in particular:
%   - When s<0, use recursion, since gammainc(x,s) only take s>0
%   - When s=0, use expint(x) == igamma(0,x)
% When s is complex, we fall back to igamma(s,x) (Symbolic Toolbox
% required), so still slow.

% Bowei Wu 5/2020. Updated 11/2020.

if nargin == 0, test_incgamma; return; end % unit test


gs = igam_mod_fast(s,x); % customized incomplete gamma
if nargout > 1
    ex = exp(-x);
    gsp1 = (s.*gs+ex)./x;
    gsp2 = ((s+1).*gsp1+ex)./x;
    gsp3 = ((s+2).*gsp2+ex)./x;
    gsp4 = ((s+3).*gsp3+ex)./x;
end

end

function test_incgamma
% unittest
x = 2.^(-5:5)+1i*randn(1,11);
for s = [-4.5:0.5,1:0.5:5]
    g = cell(1,5);
    [g{1},g{2},g{3},g{4},g{5}] = incgamma_cmpl(s,x);
    
    % verify results
    err_abs = zeros(1,5);
    err_rel = zeros(1,5);
    g_exact = cell(1,5);
    for i = 1:5
        ss = s+i-1;
        g_exact{i} = igamma(ss,x).*x.^(-ss);
        err_abs(i) = norm(g_exact{i}-g{i});
        err_rel(i) = norm((g_exact{i}-g{i})./g_exact{i});
    end
    fprintf('s = %.1f, x = [',s)
    fprintf('%.2f ',x)
    fprintf([']:\n\t',...
              'g(s,x):\t\tabs.err = %.1e \trel.err = %.1e\n\t',...
              'g(s+1,x):\tabs.err = %.1e \trel.err = %.1e\n\t',...
              'g(s+2,x):\tabs.err = %.1e \trel.err = %.1e\n\t',...
              'g(s+3,x):\tabs.err = %.1e \trel.err = %.1e\n\t',...
              'g(s+4,x):\tabs.err = %.1e \trel.err = %.1e\n'],...
              [err_abs;err_rel]);
end
end


function [Fout] = igam_mod_fast(S,X)
% This function returns Gamma(s,x)*x^(-s), where Gamma(s,x) is the upper
% incomplete gamma function.
% REQUIRE: x is not in (-infinity,0), s is not negative integers or 0.
% The algorithm evaluates the continued fraction representations of
% Gamma(s,x) using the "modified Lentz algorithm".
% 
% implemented by J. Hoskins (B. Wu modified 6/26/23)

itmax = 1000;
eps   = 10^(-14);

Fout = zeros(size(X));
inds = find(abs(X)<1);

if (numel(inds) ~= 0)
xs = X(inds);

F = 10^(-30);
D = 0;
C = 10^(-30);

b = S;
a = 1;
D = 1./(b+a.*D);
C = b + a./C;
F = F.*C.*D;

for ii=1:itmax
   
    b = S+2*ii-1;
    a = -xs.*(S+ii-1);
    D = 1./(b+a.*D);
    C = b + a./C;
    F = F.*C.*D;
    if (max(abs(C.*D-1))<eps)
        break
    end
    
    b = S+2*ii;
    a = xs*ii;
    D = 1./(b+a.*D);
    C = b + a./C;
    F = F.*C.*D;
    if (max(abs(C.*D-1))<eps)
        break
    end
    
end   

F = F.*exp(-xs);

if (ii >itmax -1)
    disp("bombing in fast gamma");    
end


F = -F+gamma(S).*(xs.^(-S));
Fout(inds) = F;

end

inds = find(abs(X)>=1);
if (numel(inds)~=0)
xs = X(inds);

F = 10^(-120);
D = 0;
C = 10^(-120);

b = xs;
a = 1;
D = 1./(b+a.*D);
C = b + a./C;
F = F.*C.*D;

for ii=1:itmax
   
    b = 1;
    a = ii-S;
    D = 1./(b+a.*D);
    C = b + a./C;
    F = F.*C.*D;
    if (max(abs(C.*D-1))<eps)
        break
    end
    
    b = xs;
    a = ii;
    D = 1./(b+a.*D);
    C = b + a./C;
    F = F.*C.*D;
    if (max(abs(C.*D-1))<eps)
        break
    end
    
end   

F = F.*exp(-xs);
Fout(inds) = F;

if (ii >itmax -1)
    disp("bombing in fast gamma");    
end
end

end
