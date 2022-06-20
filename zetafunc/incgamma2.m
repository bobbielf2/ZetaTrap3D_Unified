function gs = incgamma2(s,x,K)
% My customed fast implementation of the (scaled) upper incomplet gamma 
% function for all real s, including s<=0. (Still slow for complex s.)
%
% Input:
%   (s,x) are the two arguments of the incomplete gamma function defined as
%       incgamma(s,x) = \int_1^\infty t^{s-1}*exp(-x*t) dt
%                     = x^{-s} \int_x^\infty t^{s-1}*exp(-t) dt
%                     = gammainc(x,s,'upper')*gamma(s)*x^-s
%                     = gammainc(x,s,'scaledupper')*exp(-x)/s
%   x>0 and at least one of s and x must be a scalar.
% Output:
%   gs{0+1} == igamma(s,x)*x^(-s)
%   gs{1+1} == igamma(s+1,x)*x^(-(s+1))
%   gs{2+1} == igamma(s+2,x)*x^(-(s+2))
%   gs{3+1} == igamma(s+3,x)*x^(-(s+3))
%   gs{4+1} == igamma(s+3,x)*x^(-(s+4))
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
if nargin < 2, K = 0; end

if isreal(s)
    if numel(s) == 1 && s <= 0
        K0 = -floor(s);
        gs = cell(max(K,K0)+1,1);
        ss = s + K0;
        ex = exp(-x);
        if ss == 0  % s is integer
            % expint(x) is equivalent to (but MUCH faster than) igamma(0,x) from the symbolic toolbox
            gs{K0+1} = expint(x).*ones(size(x));
        else
            % gammainc(x,s,'upper').*gamma(s) is equivalent to (but MUCH faster than) igamma(s,x) from the symbolic toolbox
            gs{K0+1} = gammainc(x,ss,'scaledupper').*ex./ss;
        end
        for k = K0-1:-1:0
            ss = ss - 1;
            gs{k+1} = (gs{k+2}.*x-ex)./ss;
        end
        if K > K0
            ss = s+K0-1;
            for k = K0+1:K
                ss = ss + 1;
                gs{k+1} = (ss.*gs{k}+ex)./x;
            end
        end
    else    % handle s>0 case or numel(s)>1 and numel(x)==1 case
        % TODO: this currently only works if all(s>=0)
        gs = cell(K+1,1);
        ex = exp(-x);
        gs{0+1} = gammainc(x,s,'scaledupper').*ex./s;
        if numel(s)>1, gs{1}(s==0) = expint(x); end % correct where s=0
        if K > 0
            ss = s - 1;
            for k = 1:K
                ss = ss + 1;
                gs{k+1} = (ss.*gs{k}+ex)./x;
            end
        end
    end
else
    gs = cell(K+1,1);
    gs{0+1} = igamma(s,x).*x.^(-s); % slow
    if K > 0
        ex = exp(-x);
        ss = s - 1;
        for k = 1:K
            ss = ss + 1;
            gs{k+1} = (ss.*gs{k}+ex)./x;
        end
    end
end
end

function test_incgamma
% unittest
x = 2.^(-5:5);
for s = -5:0.5:5
    g = incgamma2(s,x,4);
    
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