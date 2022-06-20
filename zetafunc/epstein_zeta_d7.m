function [Zs,Zsd1,Zsd2,Zsd3,Zsd4,Zsd5,Zsd6,Zsd7] = epstein_zeta_d7(s,E,F,G)
% Compute all the partial derivatives of the Epstein zeta function
% 
%   (d/dE)^(n-l)*(0.5*d/dF)^l*Z(s;E,F,G), l=0,...,n; n = 1,...,7
%
% by combining mixed derivatives
%
%   (a * d/dE + b * 0.5*d/dF)^l*Z(s;E,F,G),     l=0,...,n
%   (a * 0.5*d/dF + b * d/dG)^l*Z(s;E,F,G),     l=0,...,n
%
% for different values of a's and b's.

if nargin == 0, testcode; return; end % unit test

% compute mixed derivatives of Z(s;E,F,G)
N = 7;
ZdEF_mix = cell(N+1,1);
ZdFG_mix = cell(N+1,1);
[Zs,ZdEF_mix{1}] = epstein_zeta7(s,E,F,G,1,0,0);
for k = 1:N     % NOTE: only N+1 func calls
    ak = cos(k*pi/2/N); bk = sin(k*pi/2/N);
    [~,ZdEF_mix{k+1}] = epstein_zeta7(s,E,F,G,ak,bk/2,0);
    [~,ZdFG_mix{k+1}] = epstein_zeta7(s,E,F,G,0,ak/2,bk);
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
% The matrices C{n} = A{n}^-1 are precomputed symbolically and then
% converted to double. The code for precomputation is as follows.
% N = 7;
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

C{1} = eye(2);
C{2} = [...
   1.000000000000000                   0                   0
  -0.398736694441202   1.025716863272554  -0.626980168831352
                   0                   0   1.000000000000000];
C{3} = [...
   1.000000000000000                   0                   0                   0
  -0.852698671793288   1.232185281454791  -0.593389157216838   0.333333333333333
   0.333333333333333  -0.593389157216838   1.232185281454791  -0.852698671793288
                   0                   0                   0   1.000000000000000];
C{4} = [...
   1.000000000000000                   0                   0                   0                   0
  -0.838892351065567   1.665240867117520  -2.076521396572336   1.563662964936060  -0.313490084415676
   0.506668916411926  -1.419948675607013   3.541294073615151  -2.995972804956623   0.701291823869892
  -0.199368347220601   0.639524003844966  -2.076521396572336   2.589379828208613  -0.953014088260642
                   0                   0                   0                   0   1.000000000000000];
C{5} = [...
   1.000000000000000                   0                   0                   0                   0                   0
  -1.332192693694016   1.797583682973974  -1.156033494330103   0.921905948384995  -0.410286745309022   0.200000000000000
   1.145614477609766  -2.048955269424500   3.125338361580530  -2.702793191334143   1.319629874927662  -0.666096346847008
  -0.666096346847008   1.319629874927662  -2.702793191334143   3.125338361580530  -2.048955269424500   1.145614477609766
   0.200000000000000  -0.410286745309022   0.921905948384995  -1.156033494330103   1.797583682973974  -1.332192693694016
                   0                   0                   0                   0                   0   1.000000000000000];
C{6} = [...
   1.000000000000000                   0                   0                   0                   0                   0                   0
  -1.327516524364540   2.342345859533138  -1.384347597714891   1.110160578078347  -1.299903882909098   0.768254956987495  -0.208993389610451
   1.300507364567463  -3.357811703061212   3.260735384654701  -3.182886221671515   3.891145036394545  -2.377551088541974   0.665861227657991
  -0.908303771510135   2.674420455398924  -3.023354010675701   4.593273738177510  -6.201998062569174   4.089050141682962  -1.223088490504387
   0.447732720194580  -1.397290474022737   1.688519429242268  -3.182886221671515   5.463360991806979  -4.338072317580449   1.518635872030875
  -0.132912231480401   0.426349335896644  -0.531648925921603   1.110160578078347  -2.152602554702387   2.684251480623989  -1.403597682494590
                   0                   0                   0                   0                   0                   0   1.000000000000000];
C{7} = [...
   1.000000000000000                   0                   0                   0                   0                   0                   0                   0
  -1.317008497692849   2.885095622584175  -3.324733520149825   3.343345757329796  -2.666229271303400   1.601107277602765  -0.658504248846425   0.142857142857143
   1.404548994446753  -4.652479443148574   7.915678599239561  -8.876693479038767   7.484626452401385  -4.663216868454475   1.973496847608578  -0.439002832564283
  -1.217938435436870   4.789177676897046  -9.750709126798958  13.044137339108348 -12.147094243380247   8.097689058879354  -3.614323552018159   0.842729396668052
   0.842729396668052  -3.614323552018159   8.097689058879354 -12.147094243380247  13.044137339108348  -9.750709126798958   4.789177676897046  -1.217938435436870
  -0.439002832564283   1.973496847608578  -4.663216868454475   7.484626452401385  -8.876693479038767   7.915678599239561  -4.652479443148574   1.404548994446753
   0.142857142857143  -0.658504248846425   1.601107277602765  -2.666229271303400   3.343345757329796  -3.324733520149825   2.885095622584175  -1.317008497692849
                   0                   0                   0                   0                   0                   0                   0   1.000000000000000];

function testcode
% compute mixed derivatives of Z(s;E,F,G)
ru = randn(3,10); rv = randn(3,10);
E = dot(ru,ru); F = dot(ru,rv); G = dot(rv,rv);
s = -3;
[Zs,Zsd1,Zsd2,Zsd3,Zsd4,Zsd5,Zsd6,Zsd7] = epstein_zeta_d7(s,E,F,G);
Zd_new = {Zs,Zsd1,Zsd2,Zsd3,Zsd4,Zsd5,Zsd6,Zsd7};

% old ways of computing the first 4 orders of derivatives
[Zs,Zsd1,Zsd2,Zsd3,Zsd4] = epstein_zeta_d4_s2d1(s,E,F,G);
Zd_old = {Zs,Zsd1,Zsd2,Zsd3,Zsd4};
% print err
for n = 1:5
    err = norm(Zd_new{n} - Zd_old{n})/norm(Zd_new{n});
    %disp(Zd_new{n})
    fprintf('relative error (n=%d): %.4g\n\n',n,err)
end
keyboard