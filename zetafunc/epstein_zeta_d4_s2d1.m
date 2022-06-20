function [Zs,Zsd1,Zsd2,Zsd3,Zsd4,Zs2,Zs2d1] = epstein_zeta_d4_s2d1(s,E,F,G)
% compute up to 4th derivative of Z(s;E,F,G) 
%     and up to 1st derivaitve of Z(s+2;E,F,G)
% Output:
%   Zs = Z(s)
%   Zsd1 = 1st derivs of Z(s)
%   Zsd2 = 2nd derivs of Z(s)
%   Zsd3 = 3rd derivs of Z(s)
%   Zsd4 = 4th derivs of Z(s)
%	Zs2 = Z(s+2)
%   Zs2d1 = 1st derivs of Z(s+2)

[Zs, d1E,   d2E,    d3E,    d4E, Zs2, d1Es2] = epstein_zeta_int(s,E,F,G,1,0,0);
[ ~, d1G,   d2G,    d3G,    d4G,   ~, d1Gs2] = epstein_zeta_int(s,E,F,G,0,0,1);
[ ~, d1F,   d2F,    d3F,    d4F,   ~, d1Fs2] = epstein_zeta_int(s,E,F,G,0,0.5,0);
[ ~,   ~, d2EpF,  d3EpF,  d4EpF] = epstein_zeta_int(s,E,F,G,1,.5,0);
[ ~,   ~,     ~,      ~, d4Ep2F] = epstein_zeta_int(s,E,F,G,1,1,0);
[ ~,   ~,     ~,  d3EmF,  d4EmF] = epstein_zeta_int(s,E,F,G,1,-.5,0);
[ ~,   ~, d2GpF,  d3GpF,  d4GpF] = epstein_zeta_int(s,E,F,G,0,.5,1);
[ ~,   ~,     ~,      ~, d4Gp2F] = epstein_zeta_int(s,E,F,G,0,1,1);
[ ~,   ~,     ~,  d3GmF,  d4GmF] = epstein_zeta_int(s,E,F,G,0,-.5,1);
Zsd1 = [d1E;
        d1F;
        d1G];
Zsd2 = [d2E;
       (d2EpF-d2E-d2F)/2;
        d2F; 
       (d2GpF-d2F-d2G)/2;
        d2G];
Zsd3 = [d3E;
       (d3EpF-d3EmF-2*d3F)/6;
       (d3EpF+d3EmF-2*d3E)/6;
        d3F;
       (d3GpF+d3GmF-2*d3G)/6;
       (d3GpF-d3GmF-2*d3F)/6;
        d3G];
Zsd4 = [d4E;
       (6*d4EpF-2*d4EmF-d4Ep2F-3*d4E+12*d4F)/24;
       (d4EpF+d4EmF-2*d4E-2*d4F)/12;
       (d4Ep2F-3*d4EpF-d4EmF+3*d4E-12*d4F)/24;
        d4F;
       (d4Gp2F-3*d4GpF-d4GmF+3*d4G-12*d4F)/24;
       (d4GpF+d4GmF-2*d4G-2*d4F)/12;
       (6*d4GpF-2*d4GmF-d4Gp2F-3*d4G+12*d4F)/24;
        d4G];
Zs2d1 = [d1Es2;
         d1Fs2;
         d1Gs2];
end