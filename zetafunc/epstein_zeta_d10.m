function [Zs,Zsd1,Zsd2,Zsd3,Zsd4,Zsd5,Zsd6,Zsd7,Zsd8,Zsd9,Zsd10] = epstein_zeta_d10(s,E,F,G)
% Compute all the partial derivatives of the Epstein zeta function
% 
%   (d/dE)^(n-l)*(0.5*d/dF)^l*Z(s;E,F,G), l=0,...,n; n = 1,...,8
%
% by combining mixed derivatives
%
%   (a * d/dE + b * 0.5*d/dF)^l*Z(s;E,F,G),     l=0,...,n
%   (a * 0.5*d/dF + b * d/dG)^l*Z(s;E,F,G),     l=0,...,n
%
% for different values of a's and b's.

if nargin == 0, testcode; return; end % unit test

% compute mixed derivatives of Z(s;E,F,G)
N = 10;
ZdEF_mix = cell(N+1,1);
ZdFG_mix = cell(N+1,1);
[Zs,ZdEF_mix{1}] = epstein_zeta10(s,E,F,G,1,0,0);
for k = 1:N     % NOTE: only N+1 func calls
    ak = cos(k*pi/2/N); bk = sin(k*pi/2/N);
    [~,ZdEF_mix{k+1}] = epstein_zeta10(s,E,F,G,ak,bk/2,0);
    [~,ZdFG_mix{k+1}] = epstein_zeta10(s,E,F,G,0,ak/2,bk);
end
ZdFG_mix{1} = ZdEF_mix{N+1};
% solve for partial derivatives
C = unmixing_mat; % get matrix for unmixing the derivatives
for n=1:N
    K = round(linspace(0,N,n+1))+1;
    ZdEF_n = zeros(n+1,size(Zs,2));
    ZdFG_n = ZdEF_n;
    for i=1:n+1
        ZdEF_n(i,:) = ZdEF_mix{K(i)}(n,:);
        ZdFG_n(i,:) = ZdFG_mix{K(i)}(n,:);
    end
    ZdEF = C{n}*ZdEF_n;
    ZdFG = C{n}(2:end,:)*ZdFG_n;
    eval(['Zsd',num2str(n),' = [ZdEF;ZdFG];']);
end


% Algorithm description: 
% Use Only N+1 function calls to get mixed derivatives, from which compute
% all partial derivatives of orders 1 through N.
%
% Steps:
% 1. Pick N+1 values a_k = cos(k*pi/(2*N)), b_k = sin(k*pi/(2*N))
% 2. N+1 function calls to compute
%       (a_k*d/dx+b_k*d/dy)^n*f(x,y) = c_{n,k},   k=0,1,...,N
%    In the k-th function call, the mixed derivatives are computed for all
%    n = 1,2,...,N
% 3. There are n+1 unknown n-th derivatives
%
%       (d/dx)^(n-l)*(d/dy)^l*f,  l=0,1,...,n 
%
%    but N+1 mixed n-th derivatives are available: c_{n,k}, k=0,1,...,N. 
%    We pick n+1 of the c_{n,k} as follows:
%
%      (a_k*d/dx+b_k*d/dy)^n*f(x,y) = c_{n,k},  k=round(N/n*i), i=0,1,...,n
%
%    This picks out n+1 equations for our n+1 unknowns.
%
% Example (N=7): suppose we want (d/dx)^(n-l)*(d/dy)^l*f(x,y), n=1,...,7
% and l = 0,1,...,n. 
% 1. Define a_k = cos(k*pi/14), b_k = sin(k*pi/14), k=0,1,...,7
% 2. For each k evaluate
%       c_{n,k} = (a_k*d/dx+b_k*d/dy)^n*f(x,y),  n=1,...,7
% 3. For each n=1,...,7, we solve
%       A{n}(i+1,:)*x = c_{n,k(i)},  i=0,...,n
%    where the matrix A{n} is
%       A{n}(i,l) := nchoosek(n,l) * a_k(i)^(n-l) * b_k(i)^l,  
%           l=0,...,n and k(i)=round(7/n*i), i=0,1,...,n.
%    Then the solution x satisfies 
%           x(l+1)=(d/dx)^(n-l)*(d/dy)^l*f(x,y), l=0,...,n
% Specifically, we solve the following 7 linear systems:
%   (n=1):  A{1}(i,:)*x = c_{n,k(i)}, k(i)=0,7
%   (n=2):  A{2}(i,:)*x = c_{n,k(i)}, k(i)=0,4,7
%   (n=3):  A{3}(i,:)*x = c_{n,k(i)}, k(i)=0,2,5,7
%   (n=4):  A{4}(i,:)*x = c_{n,k(i)}, k(i)=0,2,4,5,7
%   (n=5):  A{5}(i,:)*x = c_{n,k(i)}, k(i)=0,1,3,4,6,7
%   (n=6):  A{6}(i,:)*x = c_{n,k(i)}, k(i)=0,1,2,4,5,6,7
%   (n=7):  A{7}(i,:)*x = c_{n,k(i)}, k(i)=0,1,2,3,4,5,6,7
%
% Note: the set of k(i) is given by round(linspace(0,N,n+1))


function C = unmixing_mat
% The matrices C{n} = A{n}^-1 are precomputed symbolically in MATHEMATICA
% and then converted to double. The code for precomputation is as follows.
% N = 9;
% C = cell(N,1);
% for n = 1:N
%     A = sym(zeros(n+1));
%     i = 0;
%     for k = sym(round(linspace(0,N,n+1)))
%         a = cos(k*pi/2/N); b = sin(k*pi/2/N);
%         for l = 0:n
%             A(i+1,l+1) = nchoosek(n,l)*a^(n-l)*b^l;
%         end
%         i = i+1;
%     end
%     C{n} = double(A^-1);
% end
% This code is VERY SLOW for n=8, so need to do it in MATHEMATICA:
% (WARNING: the rounding can be different in different software! e.g. in
% MATLAB: round(5/2)==3, but in MATHEMATICA: Round[5/2]==2)
%
% Unprotect[Power]; Power[0, 0] = 1; Protect[Power];(*set 0^0=1*)
% NN = 10;
% n = 4;
% K = Round[Subdivide[0, NN, n]]; (* K[[2]]++ when NN=10 and n=4 due to different roundings in MATLAB*)
% A = Table[0, {i, n+1}, {j, n+1}];
% For[i = 1, i <= n + 1, i++,
%   k = K[[i]];
%   a = Cos[k Pi/2/NN];
%   b = Sin[k Pi/2/NN];
%   For[l = 0, l <= n, l++,
%   	A[[i, l+1]] = Binomial[n, l] a^(n-l) b^l;
%   ]
% ]
% Inverse[N[A, 50]] // MatrixForm (*show numerical values*)
%
% Then copy the output to MATLAB via Right-click->"Copy As"->"Plain Text"

C{1} = eye(2);

C{2} = [...
     1.000000000000000e+00                         0                         0
    -5.000000000000000e-01     1.000000000000000e+00    -5.000000000000000e-01
                         0                         0     1.000000000000000e+00];

C{3} = [...
     1.000000000000000e+00                         0                         0                         0
    -8.240453183331932e-01     1.249146268978545e+00    -6.364718141855817e-01     3.333333333333333e-01
     3.333333333333333e-01    -6.364718141855817e-01     1.249146268978545e+00    -8.240453183331932e-01
                         0                         0                         0     1.000000000000000e+00];
C{4} = [...
     1.000000000000000e+00                         0                         0                         0                         0
    -8.218825504345142e-01     1.902113032590307e+00    -1.538841768587627e+00     8.506508083520399e-01    -3.920395219202061e-01
     4.875368351683835e-01    -1.680098014226801e+00     2.346764680893468e+00    -1.680098014226801e+00     8.592278457250845e-01
    -1.594227023180611e-01     6.180339887498949e-01    -9.813052527525753e-01     1.669496212988162e+00    -1.146802246667421e+00
                         0                         0                         0                         0     1.000000000000000e+00];

C{5} = [...
     1.000000000000000e+00                         0                         0                         0                         0                         0
    -1.101105536376939e+00     2.094427190999916e+00    -1.781626183058571e+00     1.294427190999916e+00    -6.805206466816319e-01     2.000000000000000e-01
     9.155417527999327e-01    -2.542471396094693e+00     3.678297101099823e+00    -3.093024164283163e+00     1.762755348299891e+00    -5.505527681884694e-01
    -5.505527681884694e-01     1.762755348299891e+00    -3.093024164283163e+00     3.678297101099823e+00    -2.542471396094693e+00     9.155417527999327e-01
     2.000000000000000e-01    -6.805206466816319e-01     1.294427190999916e+00    -1.781626183058571e+00     2.094427190999916e+00    -1.101105536376939e+00
                         0                         0                         0                         0                         0     1.000000000000000e+00];

C{6} = [...
     1.000000000000000e+00                         0                         0                         0                         0                         0                         0
    -1.145789864734623e+00     3.797145017283432e+00    -3.797145017283432e+00     2.013431347560135e+00    -1.934742021726871e+00     1.233767205568027e+00    -1.666666666666667e-01
     1.085762465487592e+00    -5.767188378570133e+00     7.460825981015437e+00    -4.731353616272051e+00     4.926018599179518e+00    -3.232380996734214e+00     4.583159458938493e-01
    -8.143218491156938e-01     5.240403364622822e+00    -7.570567323443125e+00     6.288971615871994e+00    -7.570567323443125e+00     5.240403364622822e+00    -8.143218491156938e-01
     4.583159458938493e-01    -3.232380996734214e+00     4.926018599179518e+00    -4.731353616272051e+00     7.460825981015437e+00    -5.767188378570133e+00     1.085762465487592e+00
    -1.666666666666667e-01     1.233767205568027e+00    -1.934742021726871e+00     2.013431347560135e+00    -3.797145017283432e+00     3.797145017283432e+00    -1.145789864734623e+00
                         0                         0                         0                         0                         0                         0     1.000000000000000e+00];

C{7} = [...
     1.000000000000000e+00                         0                         0                         0                         0                         0                         0                         0
    -1.578170908353670e+00     2.301417420745564e+00    -2.215368586762247e+00     2.156105238364911e+00    -1.566502150527244e+00     1.128786674965872e+00    -3.645087101379241e-01     1.428571428571428e-01
     1.800434670159959e+00    -3.631210807887028e+00     6.708568707471886e+00    -6.950431223012186e+00     5.389108806759670e+00    -3.964927502411185e+00     1.323022262406578e+00    -5.260569694512235e-01
    -1.592605254073598e+00     3.646979559205116e+00    -8.852446472005632e+00     1.056422309449650e+01    -9.496366045599940e+00     7.323549111529842e+00    -2.630623616751042e+00     1.080260802095975e+00
     1.080260802095975e+00    -2.630623616751042e+00     7.323549111529842e+00    -9.496366045599940e+00     1.056422309449650e+01    -8.852446472005632e+00     3.646979559205116e+00    -1.592605254073598e+00
    -5.260569694512235e-01     1.323022262406578e+00    -3.964927502411185e+00     5.389108806759670e+00    -6.950431223012186e+00     6.708568707471886e+00    -3.631210807887028e+00     1.800434670159959e+00
     1.428571428571428e-01    -3.645087101379241e-01     1.128786674965872e+00    -1.566502150527244e+00     2.156105238364911e+00    -2.215368586762247e+00     2.301417420745564e+00    -1.578170908353670e+00
                         0                         0                         0                         0                         0                         0                         0     1.000000000000000e+00];

C{8} = [...
     1.000000000000000e+00                         0                         0                         0                         0                         0                         0                         0                         0
    -1.622207954610534e+00     3.156875757337521e+00    -2.260073510670101e+00     2.422533358786554e+00    -3.097858762962574e+00     1.760073510670101e+00    -7.343423985509674e-01     5.000000000000000e-01    -1.250000000000000e-01
     2.063057315216251e+00    -6.010612179322845e+00     6.392766508432429e+00    -8.029854846483220e+00     1.060146007352575e+01    -6.160821068499740e+00     2.654698952465470e+00    -1.831325599508534e+00     4.634879870315812e-01
    -2.072394705695320e+00     7.076506621954125e+00    -8.813188636453791e+00     1.446522715144631e+01    -2.026351069298394e+01     1.228648086917695e+01    -5.628679840203635e+00     3.981087890367424e+00    -1.031528657608125e+00
     1.657915764556256e+00    -6.127229156306119e+00     8.276727558020619e+00    -1.620310826415762e+01     2.487710248148801e+01    -1.620310826415762e+01     8.276727558020619e+00    -6.127229156306119e+00     1.657915764556256e+00
    -1.031528657608125e+00     3.981087890367424e+00    -5.628679840203635e+00     1.228648086917695e+01    -2.026351069298394e+01     1.446522715144631e+01    -8.813188636453791e+00     7.076506621954125e+00    -2.072394705695320e+00
     4.634879870315812e-01    -1.831325599508534e+00     2.654698952465470e+00    -6.160821068499740e+00     1.060146007352575e+01    -8.029854846483220e+00     6.392766508432429e+00    -6.010612179322845e+00     2.063057315216251e+00
    -1.250000000000000e-01     5.000000000000000e-01    -7.343423985509674e-01     1.760073510670101e+00    -3.097858762962574e+00     2.422533358786554e+00    -2.260073510670101e+00     3.156875757337521e+00    -1.622207954610534e+00
                         0                         0                         0                         0                         0                         0                         0                         0     1.000000000000000e+00];

C{9} = [...
     1.000000000000000e+00                         0                         0                         0                         0                         0                         0                         0                         0                         0
    -1.605533287987095e+00     3.774211998692576e+00    -4.716740797303454e+00     4.577968182461101e+00    -2.713395656996982e+00     1.971397340113354e+00    -2.332591295939686e+00     1.532561987069195e+00    -5.977764550790732e-01     1.111111111111111e-01
     2.122176503411205e+00    -7.676792568735288e+00     1.340983093672841e+01    -1.429148858145430e+01     8.868338682929689e+00    -6.763498118355465e+00     8.129240531887660e+00    -5.411813499998133e+00     2.135767871840621e+00    -4.013833219967737e-01
    -2.326640972450079e+00     1.012137235455587e+01    -2.092139679608287e+01     2.545229269816214e+01    -1.697936551792863e+01     1.403096100240587e+01    -1.731834815429308e+01     1.179124223671129e+01    -4.748147968417145e+00     9.095042157476593e-01
     2.105423719790445e+00    -1.008486429311362e+01     2.291868204325486e+01    -3.060577203758421e+01     2.229852391359126e+01    -2.072434302889928e+01     2.667985546813741e+01    -1.884019117443428e+01     7.843511984147060e+00    -1.551093981633386e+00
    -1.551093981633386e+00     7.843511984147060e+00    -1.884019117443428e+01     2.667985546813741e+01    -2.072434302889928e+01     2.229852391359126e+01    -3.060577203758421e+01     2.291868204325486e+01    -1.008486429311362e+01     2.105423719790445e+00
     9.095042157476593e-01    -4.748147968417145e+00     1.179124223671129e+01    -1.731834815429308e+01     1.403096100240587e+01    -1.697936551792863e+01     2.545229269816214e+01    -2.092139679608287e+01     1.012137235455587e+01    -2.326640972450079e+00
    -4.013833219967737e-01     2.135767871840621e+00    -5.411813499998133e+00     8.129240531887660e+00    -6.763498118355465e+00     8.868338682929689e+00    -1.429148858145430e+01     1.340983093672841e+01    -7.676792568735288e+00     2.122176503411205e+00
     1.111111111111111e-01    -5.977764550790732e-01     1.532561987069195e+00    -2.332591295939686e+00     1.971397340113354e+00    -2.713395656996982e+00     4.577968182461101e+00    -4.716740797303454e+00     3.774211998692576e+00    -1.605533287987095e+00
                         0                         0                         0                         0                         0                         0                         0                         0                         0     1.000000000000000e+00];
               
C{10} = [...
     1.000000000000000e+00                         0                         0                         0                         0                         0                         0                         0                         0                         0                         0
    -1.544979959188385e+00     4.086345818906140e+00    -6.611846424776157e+00     9.427964041848320e+00    -1.103845256702311e+01     1.049819224368230e+01    -8.019905233312238e+00     4.803787616240077e+00    -2.148319131876897e+00     6.472135954999579e-01    -1.000000000000000e-01
     2.018847860326383e+00    -8.296233747028378e+00     1.817834028960906e+01    -2.825705194044166e+01     3.452194520484305e+01    -3.371039422183971e+01     2.623980586177789e+01    -1.594886753520937e+01     7.220681966037149e+00    -2.199291506782942e+00     3.433288798196412e-01
    -2.265301631738417e+00     1.129378002904672e+01    -2.907597803706018e+01     5.057949875185502e+01    -6.575031833364237e+01     6.683705072350899e+01    -5.356700636915316e+01     3.332055360198071e+01    -1.538443238845494e+01     4.769221601280011e+00    -7.570679476223936e-01
     2.193910620854298e+00    -1.214967871391707e+01     3.445234629517893e+01    -6.531655295207193e+01     9.117525525702634e+01    -9.770209684323544e+01     8.157504907941316e+01    -5.248149097398974e+01     2.495268757199056e+01    -7.946268368909144e+00     1.294458075279096e+00
    -1.828258850711915e+00     1.078393545250693e+01    -3.252048494833859e+01     6.554179286442887e+01    -9.723484192980722e+01     1.105157148238439e+02    -9.723484192980722e+01     6.554179286442887e+01    -3.252048494833859e+01     1.078393545250693e+01    -1.828258850711915e+00
     1.294458075279096e+00    -7.946268368909144e+00     2.495268757199056e+01    -5.248149097398974e+01     8.157504907941316e+01    -9.770209684323544e+01     9.117525525702634e+01    -6.531655295207193e+01     3.445234629517893e+01    -1.214967871391707e+01     2.193910620854298e+00
    -7.570679476223936e-01     4.769221601280011e+00    -1.538443238845494e+01     3.332055360198071e+01    -5.356700636915316e+01     6.683705072350899e+01    -6.575031833364237e+01     5.057949875185502e+01    -2.907597803706018e+01     1.129378002904672e+01    -2.265301631738417e+00
     3.433288798196412e-01    -2.199291506782942e+00     7.220681966037149e+00    -1.594886753520937e+01     2.623980586177789e+01    -3.371039422183971e+01     3.452194520484305e+01    -2.825705194044166e+01     1.817834028960906e+01    -8.296233747028378e+00     2.018847860326383e+00
    -1.000000000000000e-01     6.472135954999579e-01    -2.148319131876897e+00     4.803787616240077e+00    -8.019905233312238e+00     1.049819224368230e+01    -1.103845256702311e+01     9.427964041848320e+00    -6.611846424776157e+00     4.086345818906140e+00    -1.544979959188385e+00
                         0                         0                         0                         0                         0                         0                         0                         0                         0                         0     1.000000000000000e+00];
               
function testcode
% compute mixed derivatives of Z(s;E,F,G)
ru = randn(3,10); rv = randn(3,10);
E = dot(ru,ru); F = dot(ru,rv); G = dot(rv,rv);
s = -1;
[Zs,Zsd1,Zsd2,Zsd3,Zsd4,Zsd5,Zsd6,Zsd7,Zsd8,Zsd9,Zsd10] = epstein_zeta_d10(s,E,F,G);
Zd_new = {Zs,Zsd1,Zsd2,Zsd3,Zsd4,Zsd5,Zsd6,Zsd7,Zsd8,Zsd9,Zsd10};

% old ways of computing the first 4 orders of derivatives
[Zs,Zsd1,Zsd2,Zsd3,Zsd4,Zsd5,Zsd6,Zsd7] = epstein_zeta_d7(s,E,F,G);
Zd_old = {Zs,Zsd1,Zsd2,Zsd3,Zsd4,Zsd5,Zsd6,Zsd7};
% print err
for n = 1:8
    err = norm(Zd_new{n} - Zd_old{n})/norm(Zd_new{n});
    %disp(Zd_new{n})
    fprintf('relative error (n=%d): %.4g\n\n',n,err)
end