function [gs,gsp] = incgamma7(s,x)
% My customed fast implementation of the (scaled) upper incomplet gamma 
% function for all real s, including s<=0. (Still slow for complex s.)
%
% Input:
%   (s,x) are the two arguments of the incomplete gamma function defined as
%       incgamma(s,x) = \int_1^\infty t^{s-1}*exp(-x*t) dt
%                     = x^{-s} \int_x^\infty t^{s-1}*exp(-t) dt
%   x>0 and at least one of s and x must be a scalar.
% Output:
%   gs       = igamma(s,x)*x^(-s)
%   gsp(k,:) = igamma(s+k,x)*x^(-(s+k)), k = 1,...,7
%
% Note: this implementation is Fast in the following sense:
% When isreal(s), we use gammainc(x,s,'upper') or expint(x) (FAST) 
% instead of igamma(s,x) (SLOW); in particular:
%   - When s<0, use recursion, since gammainc(x,s) only take s>0
%   - When s=0, use expint(x) == igamma(0,x)
% When s is complex, we fall back to igamma(s,x) (Symbolic Toolbox
% required), so still slow.

% Bowei Wu 5/2020. Updated 7/2021.

if nargin == 0, test_incgamma; return; end % unit test

if numel(s) == 1 && s <= 0
    k = -floor(s);
    ss = s + k;
    gg = zeros(k+1,length(x));
    if ss == 0  % s is integer
        % expint(x) is equivalent to (but MUCH faster than) igamma(0,x) from the symbolic toolbox
        gg(k+1,:) = expint(x).*x.^(-ss);
    else
        % gammainc(x,s,'upper').*gamma(s) is equivalent to (but MUCH faster than) igamma(s,x) from the symbolic toolbox
        gg(k+1,:) = gammainc(x,ss,'upper').*gamma(ss).*x.^(-ss);
    end
    ex = exp(-x);
    for i = k:-1:1
        ss = ss - 1;
        gg(i,:) = (gg(i+1,:).*x-ex)./ss;
    end
    gs = gg(1,:);
    if nargout > 1
        if k >= 7
            gsp = gg(2:8,:);
        elseif k==0
            gsp = zeros(7,size(gg,2));
            gsp(1,:) = ex./x;
            for i = 1:6
                gsp(i+1,:) = ((s+i).*gsp(i,:)+ex)./x;
            end
        else
            gsp = [gg(2:k+1,:);zeros(7-k,size(gg,2))];
            for i = k:6
                gsp(i+1,:) = ((s+i).*gsp(i,:)+ex)./x;
            end
        end
    end
else    % handle s>0 case or numel(s)>1 and numel(x)==1 case
    % TODO: this currently only works if all(s>=0)
    gs = gammainc(x,s,'upper').*gamma(s).*x.^(-s);
    if numel(s)>1, gs(s==0) = expint(x); end % correct where s=0
    if nargout > 1
        ex = exp(-x);
        gsp(1,:) = (s.*gs+ex)./x;
        for k = 1:6
            gsp(k+1,:) = ((s+k).*gsp(k,:)+ex)./x;
        end
    end
end
end

function test_incgamma
% unittest
x = 2.^(-5:5);
for s = -5:0.5:5
    [g,gp] = incgamma7(s,x);
    
    % verify results
    err_abs = zeros(1,8);
    err_rel = zeros(1,8);
    g_exact = igamma(s,x).*x.^(-s);
    err_abs(1) = norm(g_exact-g);
    err_rel(1) = norm((g_exact-g)./g_exact);
    for i = 1:7
        ss = s+i;
        gp_exact = igamma(ss,x).*x.^(-ss);
        err_abs(i+1) = norm(gp_exact-gp(i,:));
        err_rel(i+1) = norm((gp_exact-gp(i,:))./gp_exact);
    end
    fprintf('s = %.1f, x = [',s)
    fprintf('%.2f ',x)
    fprintf([']:\n\t',...
              'g(s,x):\t\tabs.err = %.1e \trel.err = %.1e\n\t',...
              'g(s+1,x):\tabs.err = %.1e \trel.err = %.1e\n\t',...
              'g(s+2,x):\tabs.err = %.1e \trel.err = %.1e\n\t',...
              'g(s+3,x):\tabs.err = %.1e \trel.err = %.1e\n\t',...
              'g(s+4,x):\tabs.err = %.1e \trel.err = %.1e\n\t',...
              'g(s+5,x):\tabs.err = %.1e \trel.err = %.1e\n\t',...
              'g(s+6,x):\tabs.err = %.1e \trel.err = %.1e\n\t',...
              'g(s+7,x):\tabs.err = %.1e \trel.err = %.1e\n'],...
              [err_abs;err_rel]);
end
end