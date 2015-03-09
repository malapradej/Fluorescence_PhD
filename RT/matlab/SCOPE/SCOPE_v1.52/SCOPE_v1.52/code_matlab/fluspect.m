 function [leafopt] = fluspect(spectral,leafbio,optipar)
%
% function [leafopt] = fluspect(spectral,leafbio,optipar)
% calculates reflectance and transmittance spectra of a leaf using FLUSPECT, 
% plus four excitation-fluorescence matrices
%
% Authors: Wout Verhoef, Christiaan van der Tol (tol@itc.nl), Joris Timmermans, 
% Date: 2007
% Update from PROSPECT to FLUSPECT: January 2011 (CvdT)
%
%      Nov 2012 (CvdT) Output EF-matrices separately for PSI and PSII
%   31 Jan 2013 (WV)   Adapt to SCOPE v_1.40, using structures for I/O
%   30 May 2013 (WV)   Repair bug in s for non-conservative scattering
%   24 Nov 2013 (WV)   Simplified doubling routine
%   25 Nov 2013 (WV)   Restored piece of code that takes final refl and
%                      tran outputs as a basis for the doubling routine
%
% usage:
% [leafopt] = fluspect(spectral,leafbio,optipar)
% 
% inputs:
% Cab         = leafbio.Cab;
% Cw          = leafbio.Cw;
% Cdm         = leafbio.Cdm;
% Cs          = leafbio.Cs;
% N           = leafbio.N;
% fqe         = leafbio.fqe;
% 
% nr          = optipar.nr;
% Kdm         = optipar.Kdm;
% Kab         = optipar.Kan;
% Kw          = optipar.Kw;
% Ks          = optipar.Ks;
% phiI        = optipar.phiI;
% phiII       = optipar.phiII;
%
% outputs:
% refl          reflectance
% tran          transmittance
% Mb            backward scattering fluorescence matrix, I for PSI and II for PSII
% Mf            forward scattering fluorescence matrix,  I for PSI and II for PSII

%% parameters
% fixed parameters for the fluorescence module
ndub        = 15;           % number of doublings applied

% Fluspect parameters
Cab         = leafbio.Cab;
Cw          = leafbio.Cw;
Cdm         = leafbio.Cdm;
Cs          = leafbio.Cs;
N           = leafbio.N;
fqe         = leafbio.fqe;

nr          = optipar.nr;
Kdm         = optipar.Kdm;
Kab         = optipar.Kab;
Kw          = optipar.Kw;
Ks          = optipar.Ks;
phiI        = optipar.phiI;
phiII       = optipar.phiII;

%% PROSPECT calculations
Kall        = (Cab*Kab + Cdm*Kdm + Cw*Kw  + Cs*Ks)/N;

j           = find(Kall>0);               % Non-conservative scattering (normal case)
t1          = (1-Kall).*exp(-Kall);
t2          = Kall.^2.*expint(Kall);
tau         = ones(size(t1));
tau(j)      = t1(j)+t2(j);
kChlrel     = zeros(size(t1));
kChlrel(j)  = Cab*Kab(j)./(Kall(j)*N);

talf        = calctav(59,nr);
ralf        = 1-talf;
t12         = calctav(90,nr);
r12         = 1-t12;
t21         = t12./(nr.^2);
r21         = 1-t21;

% top layer
denom       = 1-r21.*r21.*tau.^2;
Ta          = talf.*tau.*t21./denom;
Ra          = ralf+r21.*tau.*Ta;

% deeper layers
t           = t12.*tau.*t21./denom;
r           = r12+r21.*tau.*t;

% Stokes equations to compute properties of next N-1 layers (N real)
% Normal case

D           = sqrt((1+r+t).*(1+r-t).*(1-r+t).*(1-r-t));
rq          = r.^2;
tq          = t.^2;
a           = (1+rq-tq+D)./(2*r);
b           = (1-rq+tq+D)./(2*t);

bNm1        = b.^(N-1);                  %
bN2         = bNm1.^2;
a2          = a.^2;
denom       = a2.*bN2-1;
Rsub        = a.*(bN2-1)./denom;
Tsub        = bNm1.*(a2-1)./denom;

s           = r./t;                             % Conservative scattering (CS)
s(j)        = 2*a(j)./(a(j).^2-1).*log(b(j));   % Normal case overwrites CS case 

%			Case of zero absorption
j           = find(r+t >= 1);

Tsub(j)     = t(j)./(t(j)+(1-t(j))*(N-1));
Rsub(j)	    = 1-Tsub(j);

% Reflectance and transmittance of the leaf: combine top layer with next N-1 layers
denom       = 1-Rsub.*r;
tran        = Ta.*Tsub./denom;
refl        = Ra+Ta.*Rsub.*t./denom;

leafopt.refl = refl;
leafopt.tran = tran;
leafopt.kChlrel = kChlrel;

t        =   tran;
r        =   refl;

I_rt     =   (r+t)<1;
D(I_rt)  =   sqrt((1 + r(I_rt) + t(I_rt)) .* ...
                  (1 + r(I_rt) - t(I_rt)) .* ...
                  (1 - r(I_rt) + t(I_rt)) .* ...
                  (1 - r(I_rt) - t(I_rt)));
a(I_rt)  =   (1 + r(I_rt).^2 - t(I_rt).^2 + D(I_rt)) ./ (2*r(I_rt));
b(I_rt)  =   (1 - r(I_rt).^2 + t(I_rt).^2 + D(I_rt)) ./ (2*t(I_rt));    
a(~I_rt) =   1;
b(~I_rt) =   1;

I_a      =   a>1;
s(I_a)   =   2.*a(I_a) ./ (a(I_a).^2 - 1) .* log(b(I_a));
s(~I_a)  =   r(~I_a) ./ t(~I_a);

k        =   (a-1) ./ (a+1) .* log(b);
kChl     =   kChlrel .* k;
k(j)        = 0;
kChl(j)     = 0;

if fqe > 0

    %% Fluorescence
    wle         = spectral.wlE';    % excitation wavelengths, transpose to column
    wlf         = spectral.wlF';    % fluorescence wavelengths, transpose to column
    wlp         = spectral.wlP;     % PROSPECT wavelengths, also a row vector

    minwle      = min(wle);
    maxwle      = max(wle);
    minwlf      = min(wlf);
    maxwlf      = max(wlf);

    % indices of wle and wlf within wlp

    Iwle        = find(wlp>=minwle & wlp<=maxwle);
    Iwlf        = find(wlp>=minwlf & wlp<=maxwlf);

    eps         = 2^(-ndub);

    % initialisations
    te          = 1-(k(Iwle)+s(Iwle)) * eps;   
    tf          = 1-(k(Iwlf)+s(Iwlf)) * eps;  
    re          = s(Iwle) * eps;
    rf          = s(Iwlf) * eps;

    sigmoid     = 1./(1+exp(-wlf/10)*exp(wle'/10));  % matrix computed as an outproduct
    [MfI,  MbI]  = deal(.5 * fqe(1) * ((.5*phiI( Iwlf))*eps) * kChl(Iwle)'.*sigmoid);
    [MfII, MbII] = deal(.5 * fqe(2) * ((.5*phiII(Iwlf))*eps) * kChl(Iwle)'.*sigmoid);

    Ih          = ones(1,length(te));     % row of ones
    Iv          = ones(length(tf),1);     % column of ones

    % Doubling routine
    
    for i = 1:ndub
        
        xe = te./(1-re.*re);  ten = te.*xe;  ren = re.*(1+ten);  
        xf = tf./(1-rf.*rf);  tfn = tf.*xf;  rfn = rf.*(1+tfn);
              
        A11  = xf*Ih + Iv*xe';           A12 = (xf*xe').*(rf*Ih + Iv*re');
        A21  = 1+(xf*xe').*(1+rf*re');   A22 = (xf.*rf)*Ih+Iv*(xe.*re)';
        
        MfnI   = MfI  .* A11 + MbI  .* A12;
        MbnI   = MbI  .* A21 + MfI  .* A22;
        MfnII  = MfII .* A11 + MbII .* A12;
        MbnII  = MbII .* A21 + MfII .* A22;
        
        te   = ten;  re  = ren;   tf   = tfn;   rf   = rfn;
        MfI  = MfnI; MbI = MbnI;  MfII = MfnII; MbII = MbnII;
        
    end

    leafopt.MbI  = MbI;
    leafopt.MbII = MbII;
    leafopt.MfI  = MfI;
    leafopt.MfII = MfII;
    
end    

return;

function tav = calctav(alfa,nr)

    rd          = pi/180;
    n2          = nr.^2;
    np          = n2+1;
    nm          = n2-1;
    a           = (nr+1).*(nr+1)/2;
    k           = -(n2-1).*(n2-1)/4;
    sa          = sin(alfa.*rd);

    b1          = (alfa~=90)*sqrt((sa.^2-np/2).*(sa.^2-np/2)+k);
    b2          = sa.^2-np/2;
    b           = b1-b2;
    b3          = b.^3;
    a3          = a.^3;
    ts          = (k.^2./(6*b3)+k./b-b/2)-(k.^2./(6*a3)+k./a-a/2);

    tp1         = -2*n2.*(b-a)./(np.^2);
    tp2         = -2*n2.*np.*log(b./a)./(nm.^2);
    tp3         = n2.*(1./b-1./a)/2;
    tp4         = 16*n2.^2.*(n2.^2+1).*log((2*np.*b-nm.^2)./(2*np.*a-nm.^2))./(np.^3.*nm.^2);
    tp5         = 16*n2.^3.*(1./(2*np.*b-nm.^2)-1./(2*np.*a-nm.^2))./(np.^3);
    tp          = tp1+tp2+tp3+tp4+tp5;
    tav         = (ts+tp)./(2*sa.^2);

return;