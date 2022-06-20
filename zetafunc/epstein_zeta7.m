function [S,Sd] = epstein_zeta7(s,E,F,G,L,M,N)
% Evaluate the Epstein zeta function Z(s) and its parametric derivatives.
% Z(s) is defined by
%       Z(s) = \sum_{(i,j)\neq(0,0)} (E*i^2+2*F*i*j+G*j^2)^(-s/2)
% for Re(s)>2, and analytically continued to the whole complex plane except
% a simple pole at s=2. See the manuscript [1], Appendix E, for more details.
%
% Output:
%   S      = Z(s)
%   Sd(k,:)= (L*{d/dE} + M*{d/dF} + N*{d/dG})^k Z(s),  k = 1,2,...,7
%
% [1] Wu, B., & Martinsson, P.G. (2020, arXiv:2007.02512). Corrected
%     Trapezoidal Rules for Boundary Integral Equations in Three
%     Dimensions.
%
% Bowei Wu, 2021/6 add 7th derivatives

if nargin == 0, test_epstein_zeta; return; end % no input for unit testing

assert(isequal(size(E),size(F)),'E and F must have the same size')
assert(isequal(size(E),size(G)),'E and G must have the same size')
if s==0 && nargout==1, S=-ones(size(E)); return; end % handle exception: Z(0)=-1.


% Preprocessing: make quadratic form has unit determinant
J = sqrt(E.*G-F.^2);
E = E./J; G = G./J; F = F./J;
if nargout > 1, L = L./J; M = M./J; N = N./J; end

% constants for computing derivatives
s1 = s/2; s2 = 1-s1;
if nargout > 1
    H   = (G.*L+E.*N-2*F.*M)/2;
    K   = L.*N-M.^2;
    Hk = computeHk(H,K);    % Hk(k,:) = (L*{d/dE} + M*{d/dF} + N*{d/dG})^k H
    DHk = bellpoly(Hk);
    [SQpSH1,SQmH1] = compute_SQ_p_sH(DHk,s1);
end

% Define the quadratic form
QA  = @(i,j) pi*(E*i^2+2*F*i*j+G*j^2);
if nargout > 1, QB = @(i,j) pi*(L*i^2+2*M*i*j+N*j^2); end

% Determine summation cutoffs based on decay of incomplete gamma func
lambda = min(min((E+G)/2 - sqrt((E-G).^2+4*F.^2)/2)); % min eigenvalue of Q
n = floor(sqrt(33/pi./lambda))+3;

% summation, exclude origin
if nargout > 1
    S = 0; Sd = 0; 
    for i = 0:n-1     % right-half ij-plane & j-axis
        for j = 1:ceil(sqrt(n^2-i^2))
            ii = [i,j]; jj = [j,-i]; % contribution from (i,j) & (j,-i)
            for k=1:2
                Q = QA(ii(k),jj(k)); R = QB(ii(k),jj(k));
                [G0,Gk] = computeGk(s1,Q);
                Qk = computeQk(SQmH1,Q,R);
                S = S + G0;
                Sd = Sd + computeDGk(bellpoly(Qk),Gk);
            end
        end
    end
else
    S = 0;
    for i = 0:n-1     % right-half ij-plane & j-axis
        for j = 1:ceil(sqrt(n^2-i^2))
            S = S + computeGk(s1,QA(i,j)) + computeGk(s1,QA(j,-i));
        end
    end
end

% Postprocessing: symmetry add-back & determinant scale-back
Cs1 = (pi./J).^s1 ./ gamma(s1); % scaling constant
S  = (2*S - 1 ./s1 - 1 ./s2).*Cs1;
% Recurrence for (D^k)Z(s)
if nargout > 1
    nck = my_binom10;
    Sd = 2*Sd.*Cs1 - SQpSH1.*S;
    for k = 2:7
        for i = 1:k-1
            %Sd(k,:) = Sd(k,:) - nchoosek(k,i)*SQpSH1(k-i,:).*Sd(i,:);
            Sd(k,:) = Sd(k,:) - nck(k+1,i+1)*SQpSH1(k-i,:).*Sd(i,:);
        end
    end
end
% handle exception: Z(0)=-1, dZ(0)=0, Z(2)=Inf, dZ(2)=Inf.
if s==0
    S(:)=-1;
    if nargout > 1, Sd(:)=0; end
end
if s==0
    S(:)=Inf;
    if nargout > 1, Sd(:)=Inf; end
end
end

function Hk = computeHk(H,K)
% compute (D^k)H, where D = L*{d/dE}+M*{d/dF}+N*{d/dG}
% Hk(k+1,:) = (D^k)H, k = 0,1,...,6

Hk = zeros(7,length(H));
Hk(1,:) = H;
Hk(2,:) = -2*H.^2 + K;
for k = 3:7
    Hk(k,:) =  -2*(k-1)*H.*Hk(k-1,:) - (k-1)*(k-2)*K.*Hk(k-2,:);
end
% Hk(3,:) =  -4*H.*Hk(2,:)- 2*K.*Hk(1,:);
% Hk(4,:) =  -6*H.*Hk(3,:)- 6*K.*Hk(2,:);
% Hk(5,:) =  -8*H.*Hk(4,:)-12*K.*Hk(3,:);
% Hk(6,:) = -10*H.*Hk(5,:)-20*K.*Hk(4,:);
% Hk(7,:) = -12*H.*Hk(6,:)-30*K.*Hk(5,:);

% Hp = cumprod(repmat(H,[7,1]),1); % H^p, p = 1,...,7
% Kp = cumprod(repmat(K,[3,1]),1); % H^p, p = 1,2,3
% Hk(3,:) =     8*Hp(3,:) -           6*H.*K;
% Hk(4,:) =   -48*Hp(4,:) +    48*Hp(2,:).*K -              6*Kp(2,:);
% Hk(5,:) =   384*Hp(5,:) -   480*Hp(3,:).*K +         120*H.*Kp(2,:);
% Hk(6,:) = -3840*Hp(6,:) +  5760*Hp(4,:).*K -  2160*Hp(2,:).*Kp(2,:) +     120*Kp(3,:);
% Hk(7,:) = 46080*Hp(7,:) - 80640*Hp(5,:).*K + 40320*Hp(3,:).*Kp(2,:) - 5040*H.*Kp(3,:);
end

function BX = bellpoly(X)
% Bell polynomials B_{n,m}(x_1,x_2,...,x_{n-m+1})
% input:    X(k,:) = x_k
% output:   BX{n}(m,:) = B_{n,m}(x_1,x_2,...), m=1,...,n
%
% Usage:
% 1. prepare for computing [(D+s*H)^k]1, where D=L*{d/dE}+M*{d/dF}+N*{d/dG}
%   and H = (G*L+E*N-2*F*M)/(2*sqrt(EG-F^2)).
%   Then input X = [H; (D^1)H; ...; (D^6)H], output BX are the components
%   associated with s^m in [(D+s*H)^k]1, m = 1,...,k
% See: compute_SQ_p_sH
%
% 2. prepare for computing (D^k)G0 given Gk, where G0 & Gk are incgamma
% values from computeGk.
% See: computeGk, computeDGk

X1k = repmat(X(1,:),[7,1]); X1k = cumprod(X1k,1); % x1^k, k = 1,...,7
X2k = repmat(X(2,:),[3,1]); X2k = cumprod(X2k,1); % x2^k, k = 1,...,3
BX = cell(7,1);
BX{1} = X(1,:);
BX{2} =[X(2,:); X1k(2,:)];
BX{3} =[X(3,:); 3*X(1,:).*X(2,:); X1k(3,:)];
BX{4} =[X(4,:); 4*X(1,:).*X(3,:)+3*X2k(2,:); 
         6*X1k(2,:).*X(2,:); X1k(4,:)];
BX{5} =[X(5,:); 5*X(1,:).*X(4,:)+10*X(2,:).*X(3,:);
         10*X1k(2,:).*X(3,:)+15*X(1,:).*X2k(2,:);
         10*X1k(3,:).*X(2,:); X1k(5,:)];
BX{6} =[X(6,:); 6*X(1,:).*X(5,:)+15*X(2,:).*X(4,:)+10*X(3,:).^2;
         15*X1k(2,:).*X(4,:)+60*X(1,:).*X(2,:).*X(3,:)+15*X2k(3,:);
         20*X1k(3,:).*X(3,:)+45*X1k(2,:).*X2k(2,:);
         15*X1k(4,:).*X(2,:); X1k(6,:)];
BX{7} =[X(7,:);
         7*X(1,:).*X(6,:)+21*X(2,:).*X(5,:)+35*X(3,:).*X(4,:);
         21*X1k(2,:).*X(5,:)+105*X(1,:).*X(2,:).*X(4,:)+...
         70*X(1,:).*X(3,:).^2+105*X(3,:).*X2k(2,:);
         35*X1k(3,:).*X(4,:)+210*X1k(2,:).*X(2,:).*X(3,:)+...
         105*X(1,:).*X2k(3,:);
         35*X1k(4,:).*X(3,:)+105*X1k(3,:).*X2k(2,:);
         21*X1k(5,:).*X(2,:); X1k(7,:)];
end

function [SQpSH1,SQmH1] = compute_SQ_p_sH(DHk,s)
% Compute [(D+s*H)^k]1 and [(D-H)^k]1 using precomputed vals from bellpoly.
% Here, D=L*{d/dE}+M*{d/dF}+N*{d/dG}, H=(G*L+E*N-2*F*M)/(2*sqrt(EG-F^2))
SQpSH1 = zeros(size(DHk{7})); SQmH1 = SQpSH1;
S = s.^(1:7);
M = (-1).^(1:7);
for k = 1:7
    SQpSH1(k,:) = S(1:k)*DHk{k};
    SQmH1(k,:) = M(1:k)*DHk{k};
end
end

function Qk = computeQk(SQmH1,Q,R)
% Compute (D^k)Q, k=1,2,...,7, using precomputed vals from compute_SQ_p_sH.
% Here, D=L*{d/dE}+M*{d/dF}+N*{d/dG}, Q=(E*i^2+2F*ij+G*j^2)/sqrt(EG-F^2)
% and R=(L*i^2+2M*ij+N*j^2)/sqrt(EG-F^2)
Qk = SQmH1.*Q;
Qk(1,:) = Qk(1,:) + R;
Qk(2:7,:) = Qk(2:7,:) + (2:7)'.*SQmH1(1:6,:).*R;
end

function [G0,Gk] = computeGk(s,Q)
% compute incomplete gamma func values
% G = incgamma(s,Q) + incgamma(1-s,Q)
% Gk(k,:) = incgamma(s+k,Q) + incgamma(1-s+k,Q), k=1,2,...,7
if nargout > 1
    [g1,g1p] = incgamma7(s,Q);
    [g2,g2p] = incgamma7(1-s,Q);
    Gk = g1p + g2p;
    G0 = g1 + g2;
else
    G0 = incgamma7(s,Q) + incgamma7(1-s,Q);
end
end

function DGk = computeDGk(DQk,Gk)
% Compute (D^k)G0 using precomputed vals from bellpoly.
% Here, D=L*{d/dE}+M*{d/dF}+N*{d/dG}, G0 & Gk are incgamma vals from
% computeGk.
DGk = zeros(size(DQk{7}));
M = (-1).^(1:7)'.*Gk;
for k = 1:7
    DGk(k,:) = sum(M(1:k,:).*DQk{k},1);
end
end

function nck = my_binom10
nck=[1     0     0     0     0     0     0     0     0     0     0
     1     1     0     0     0     0     0     0     0     0     0
     1     2     1     0     0     0     0     0     0     0     0
     1     3     3     1     0     0     0     0     0     0     0
     1     4     6     4     1     0     0     0     0     0     0
     1     5    10    10     5     1     0     0     0     0     0
     1     6    15    20    15     6     1     0     0     0     0
     1     7    21    35    35    21     7     1     0     0     0
     1     8    28    56    70    56    28     8     1     0     0
     1     9    36    84   126   126    84    36     9     1     0
     1    10    45   120   210   252   210   120    45    10     1];
end

function test_epstein_zeta
% unit testing

% gradient checking
ru = randn(3,1); rv = randn(3,1);
E = dot(ru,ru); F = dot(ru,rv); G = dot(rv,rv);
s = -1.5;
fprintf('s=%.1f, (E,F,G) = (%f,%f,%f), angle=pi*%f\n',s,E,F,G,acos(F/sqrt(E*G))/pi)

% 1. compute your func output
[S, Sd_EpG] = epstein_zeta7(s,E,F,G,1,0,1);
[~, Sd_EmG] = epstein_zeta7(s,E,F,G,1,0,-1);
[~, Sd_E] = epstein_zeta7(s,E,F,G,1,0,0);
[~, Sd_F] = epstein_zeta7(s,E,F,G,0,1,0);
[~, Sd_G] = epstein_zeta7(s,E,F,G,0,0,1);
[~, Sd_EpF] = epstein_zeta7(s,E,F,G,1,1,0);
[~, Sd_EmF] = epstein_zeta7(s,E,F,G,1,-1,0);
Sd2_EG = (Sd_EpG(2,:) - Sd_E(2,:) - Sd_G(2,:))/2;
Sd3_EEG = (Sd_EpG(3,:) - Sd_EmG(3,:) -2*Sd_G(3,:))/6;
Sd3_EFF = (Sd_EpF(3,:) + Sd_EmF(3,:) -2*Sd_E(3,:))/6;
Sd4_EEGG = (Sd_EpG(4,:)+Sd_EmG(4,:)-2*Sd_E(4,:)-2*Sd_G(4,:))/12;

% 2. finite diff approx
epsi = 1e-4;
% (d/dF) and (d^2/dF^2)
S_Fp = epstein_zeta7(s,E,F+epsi,G);
S_Fm = epstein_zeta7(s,E,F-epsi,G);
Sd1_F_approx = (S_Fp - S_Fm)/(2*epsi);
Sd2_FF_approx = (S_Fp - 2*S + S_Fm)/epsi^2;
% (d^2/dEdG), (d^3/dE^2dG) and (d^3/dEdF^2)
[S_EpGp,Sd_E_EpGp] = epstein_zeta7(s,E+epsi,F,G+epsi,1,0,0);
[S_EpGm,Sd_E_EpGm] = epstein_zeta7(s,E+epsi,F,G-epsi,1,0,0);
[S_EmGp,Sd_E_EmGp] = epstein_zeta7(s,E-epsi,F,G+epsi,1,0,0);
[S_EmGm,Sd_E_EmGm] = epstein_zeta7(s,E-epsi,F,G-epsi,1,0,0);
[~,Sd_E_Gp] = epstein_zeta7(s,E,F,G+epsi,1,0,0);
[~,Sd_E_Gm] = epstein_zeta7(s,E,F,G-epsi,1,0,0);
[~,Sd_G_Ep] = epstein_zeta7(s,E+epsi,F,G,0,0,1);
[~,Sd_G_Em] = epstein_zeta7(s,E-epsi,F,G,0,0,1);
[~,Sd_E_Fp] = epstein_zeta7(s,E,F+epsi,G,1,0,0);
[~,Sd_E_Fm] = epstein_zeta7(s,E,F-epsi,G,1,0,0);
Sd2_EG_approx = ((S_EpGp - S_EpGm)-(S_EmGp - S_EmGm))/(4*epsi^2);
Sd3_EEG_approx = ((Sd_E_EpGp(1,:) - Sd_E_EpGm(1,:))-(Sd_E_EmGp(1,:) - Sd_E_EmGm(1,:)))/(4*epsi^2);
Sd3_EEG_approx2 = (Sd_E_Gp(2,:) - Sd_E_Gm(2,:))/(2*epsi);
Sd3_EEG_approx3 = (Sd_G_Ep(1,:) - 2*Sd_G(1,:) + Sd_G_Em(1,:))/epsi^2;
Sd3_EFF_approx = (Sd_E_Fp(1,:) - 2*Sd_E(1,:) + Sd_E_Fm(1,:))/epsi^2;


% 3. verify func output against FDM approx
fprintf('\nverify your derivatives...\n')
fprintf('(d/dF)Z(s):\t\t[FDM approx, your output]=[%f, %f], the same?\n',Sd1_F_approx,Sd_F(1,:))
fprintf('(d^2/dF^2)Z(s):\t\t[FDM approx, your output]=[%f, %f], the same?\n',Sd2_FF_approx,Sd_F(2,:))
fprintf('(d^3/dE^2dG)Z(s):\t3 different FDM approxs: [%f, %f, %f], your output: %f\n',Sd3_EEG_approx,Sd3_EEG_approx2,Sd3_EEG_approx3,Sd3_EEG)

% 4. verify derivative equivalence: (d^2/dF^2)Z(s)/4 = (d^2/dEdG)Z(s) and (d^3/dE^2dG)Z(s) = (d^3/dEdF^2)Z(s)/4
% this test also implicitly verify correctness of the computed Z(s)
fprintf('\nderivatives equivalence check... (implicitly verify correctness of Z(s))\n')
fprintf('FDM approx:\t[(d^2/dF^2)Z(s)/4, (d^2/dEdG)Z(s)] = [%f, %f], the same?\n',Sd2_FF_approx/4,Sd2_EG_approx)
fprintf('your output:\t[(d^2/dF^2)Z(s)/4, (d^2/dEdG)Z(s)] = [%f, %f] the same? (should be the same as prev line, too)\n',Sd_F(2,:)/4,Sd2_EG)
fprintf('FDM approx:\t[(d^3/dE^2dG)Z(s), (d^3/dEdF^2)Z(s)/4] = [%f, %f] the same?\n',Sd3_EEG_approx,Sd3_EFF_approx/4)
fprintf('your output:\t[(d^3/dE^2dG)Z(s), (d^3/dEdF^2)Z(s)/4] = [%f, %f] the same? (should be the same as prev line, too)\n',Sd3_EEG,Sd3_EFF/4)
fprintf('your output:\t[(d^4/dE^2dG^2)Z(s), (d^4/dF^4)Z(s)/16] = [%f, %f] the same?\n',Sd4_EEGG,Sd_F(4,:)/16)
end
