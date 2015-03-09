function [rad,profiles] = RTMf(spectral,rad,soil,leafopt,canopy,gap,angles,profiles,leafbio)  

% function 'RTMf' calculates the spectrum of fluorescent radiance in the
% observer's direction and also the TOC spectral hemispherical upward Fs flux
%
% Authors:  Wout Verhoef and Christiaan van der Tol (tol@itc.nl)
% Date:     12 Dec 2007
% Update:   26 Aug 2008 CvdT        small correction to matrices
%           07 Nov 2008 CvdT        changed layout
% Update:   19 Mar 2009 CvdT        major corrections: lines 95-96,
%                                   101-107, and 119-120.
% Update:    7 Apr 2009 WV & CvdT   major correction: lines 89-90, azimuth
%                                   dependence was not there in previous verions (implicit assumption of
%                                   azimuth(solar-viewing) = 0). This has been corrected
% Update:   May-June 2012 WV & CvdT Add calculation of hemispherical Fs
%                                   fluxes
% Update:   Jan-Feb 2013 WV         Inputs and outputs via structures for
%                                   SCOPE Version 1.40
%
% Table of contents of the function:
%   0       preparations
%       0.0     globals
%       0.1     initialisations
%       0.2     geometric quantities
%       0.3     solar irradiance factor and ext. in obs dir for all leaf angle/azumith classes
%   1.0     calculation of fluorescence flux in observation direction
%
% Usage: [rad] = RTMfH(spectral,rad,soil,leafopt,canopy,gap,angles,profiles)  
%
% The input and output are structures. These structures are further
% specified in a readme file.
%
% Input:
%   spectral    information about wavelengths and resolutions
%   rad         a large number of radiative fluxes: spectrally distributed 
%               and integrated, and canopy radiative transfer coefficients.     
%   soil        soil properties
%   leafopt     leaf optical properties
%   canopy      canopy properties (such as LAI and height)
%   gap         probabilities of direct light penetration and viewing
%   angles      viewing and observation angles
%   profiles      vertical profiles of fluxes 
%
% Output:
%   rad         a large number of radiative fluxes: spectrally distributed 
%               and integrated, and canopy radiative transfer coefficients.
%               Here, fluorescence fluxes are added
%% 0.0 globals

global constants

%% 0.1 initialisations
wlS          = spectral.wlS';       % SCOPE wavelengths, make column vectors
wlF          = spectral.wlF';       % Fluorescence wavelengths
wlE          = spectral.wlE';       % Excitation wavelengths
[dummy,iwlfi]    = intersect(wlS,wlE); %#ok<ASGLU>
[dummy,iwlfo]    = intersect(wlS,wlF); %#ok<ASGLU>
nwlfo        = length(iwlfo);
nwlfi        = length(iwlfi);
nl           = canopy.nlayers;
LAI          = canopy.LAI;
litab        = canopy.litab;
lazitab      = canopy.lazitab;
lidf         = canopy.lidf;
nli          = length(litab);
nlazi        = length(lazitab);
fqe1         = leafbio.fqe(1);
fqe2         = leafbio.fqe(2);

Ps           = gap.Ps;
Po           = gap.Po;
Pso          = gap.Pso;  

etah         = zeros(nl,1);
etau         = zeros(nli,nlazi,nl);

[Mb,Mf]             = deal(zeros(nwlfo,nwlfi));
[MpluEmin   ,...
    MpluEplu   ,...
    MminEmin   ,...
    MminEplu]       = deal(zeros(nwlfo,nl+1));
piLem_              = zeros(nwlfo,nl+1,2);
[piLs,piLd]         = deal(zeros(nl+1,1));
LoF_                = zeros(nwlfo,2);
Fhem_               = zeros(nwlfo,2);
Fiprofile           = zeros(nl+1,2);

%for speed-up the calculation only uses wavelength i and wavelength o part of the spectrum
Esunf_              = rad.Esun_(iwlfi);
Eminf_              = rad.Emin_(:,iwlfi)';          % transpose into [nwlfo,nl] matrix
Epluf_              = rad.Eplu_(:,iwlfi)';
iLAI                = LAI/nl;                       % LAI of a layer        [1]

%% 0.2 geometric quantities
rho                 = leafopt.refl(iwlfo);          % [nwlfo]     leaf/needle reflection
tau                 = leafopt.tran(iwlfo);          % [nwlfo]     leaf/needle transmission
MbI                 = leafopt.MbI;
MbII                = leafopt.MbII;
MfI                 = leafopt.MfI;
MfII                = leafopt.MfII;

rs                  = soil.refl;                    % [nwlfo]     soil reflectance

deg2rad             = constants.deg2rad;
tto                 = angles.tto;
tts                 = angles.tts;
psi                 = angles.psi;

cos_tto             = cos(tto*deg2rad);             % cos observation zenith angle
sin_tto             = sin(tto*deg2rad);             % sin observation zenith angle

cos_tts             = cos(tts*deg2rad);             % cos solar angle
sin_tts             = sin(tts*deg2rad);             % sin solar angle

cos_ttli            = cos(litab*deg2rad);           % cos leaf inclinaation angles
sin_ttli            = sin(litab*deg2rad);           % sin leaf inclinaation angles
cos_phils           = cos(lazitab*deg2rad);         % cos leaf azimuth angles rel. to sun azi
cos_philo           = cos((lazitab-psi)*deg2rad);   % cos leaf azimuth angles rel. to viewing azi

%% 0.3 geometric factors for all leaf angle/azumith classes
cds                 = cos_ttli*cos_tts*ones(1,36) + sin_ttli*sin_tts*cos_phils;  % [nli,nlazi]
cdo                 = cos_ttli*cos_tto*ones(1,36) + sin_ttli*sin_tto*cos_philo;  % [nli,nlazi]
fs                  = cds/cos_tts;                                               % [nli,nlazi]
absfs               = abs(fs);                                                   % [nli,nlazi]
fo                  = cdo/cos_tto;                                               % [nli,nlazi]
absfo               = abs(fo);                                                   % [nli,nlazi]
fsfo                = fs.*fo;                                                    % [nli,nlazi]
absfsfo             = abs(fsfo);                                                 % [nli,nlazi]
foctl               = fo.*(cos_ttli*ones(1,36));                                 % [nli,nlazi]
fsctl               = fs.*(cos_ttli*ones(1,36));                                 % [nli,nlazi]
ctl2                = cos_ttli.^2*ones(1,36);                                    % [nli,nlazi]

%% 1.0 calculation of fluorescence flux in observation direction

% fluorescence efficiencies from ebal, after default fqe has been applied

etahi = profiles.etah;
etaui = profiles.etau;

% fluorescence matrices and efficiencies for PSI and PSII
[Emin_,Eplu_] = deal(zeros(61,211));
for PS = 2:-1:1
    switch PS
        case 1, Mb = MbI;  Mf = MfI;    etah(:) = 1;          etau(:) = 1;
        case 2, Mb = MbII; Mf = MfII;   etah(:) = etahi(:); etau(:) = etaui(:);
    end

    Mplu = 0.5*(Mb+Mf);    % [nwlfo,nwlfi]  
    Mmin = 0.5*(Mb-Mf);    % [nwlfo,nwlfi]

    % in-products: we convert incoming radiation to a fluorescence spectrum using the matrices.
    % resolution assumed is 1 nm
      
    for j = 1:nl+1 % nr of layers
        MpluEmin(:,j)   = Mplu * Eminf_(:,j);          % [nwlfo,nl+1] Emin_f = PAR
        MpluEplu(:,j)   = Mplu * Epluf_(:,j);          % [nwlfo,nl+1]
        MminEmin(:,j)   = Mmin * Eminf_(:,j);          % [nwlfo,nl+1]
        MminEplu(:,j)   = Mmin * Epluf_(:,j);          % [nwlfo,nl+1]
    end
    MpluEsun            = Mplu * Esunf_;               % integration by inproduct
    MminEsun            = Mmin * Esunf_;               % integration by inproduct

    xdd2 = mean(ctl2'*lidf);

    % we calculate the spectrum for all individual leaves, sunlit and
    % shaded
    
    [Fmin,Fplu] = deal(zeros(nwlfo,nl+1));
    [G1,G2]     = deal(zeros(nl+1,1));
       
    for i = 1 : nwlfo
        
        wEsuni          = absfsfo * MpluEsun(i) + fsfo  * MminEsun(i);         % [nli,nlazi]
        sfEsuni         = absfs   * MpluEsun(i) - fsctl * MminEsun(i);         % [nli,nlazi]
        sbEsuni         = absfs   * MpluEsun(i) + fsctl * MminEsun(i);         % [nli,nlazi]
        
        for j = 1:nl
            
            sigfEmini   = MpluEmin(i,j)   - ctl2 * MminEmin(i,j);              % [nli,nlazi]
            sigfEplui   = MpluEplu(i,j+1) - ctl2 * MminEplu(i,j+1);            % [nli,nlazi]
            sigbEmini   = MpluEmin(i,j)   + ctl2 * MminEmin(i,j);              % [nli,nlazi]
            sigbEplui   = MpluEplu(i,j+1) + ctl2 * MminEplu(i,j+1);            % [nli,nlazi]
            vbEmini     = absfo * MpluEmin(i,j)   + foctl * MminEmin(i,j);     % [nli,nlazi]
            vfEplui     = absfo * MpluEplu(i,j+1) - foctl * MminEplu(i,j+1);   % [nli,nlazi]
            
            piLs(j)     = mean((etau(:,:,j) .*(wEsuni + vbEmini + vfEplui))'*lidf);
            piLd(j)     = etah(j) * mean((vbEmini + vfEplui)'*lidf);
            
            Fsmin       = mean((etau(:,:,j).*(sfEsuni + sigfEmini + sigbEplui))'*lidf);
            Fsplu       = mean((etau(:,:,j).*(sbEsuni + sigbEmini + sigfEplui))'*lidf);
            Fdmin       = etah(j) * mean((sigfEmini + sigbEplui)'*lidf);
            Fdplu       = etah(j) * mean((sigbEmini + sigfEplui)'*lidf);
            Fmin(i,j+1) = Ps(j+1) * Fsmin + (1-Ps(j+1)) * Fdmin;
            Fplu(i,j)   = Ps(j)   * Fsplu + (1-Ps(j))   * Fdplu;
            
            Emin_       = Emin_+Fmin';
            Eplu_       = Eplu_+Fplu';
        end

        %   Fmin(i,1)       = 0;   % No flurescence from the sky and from the soil
        %   Fplu(i,nl+1)    = 0;
  
        piLtoti         = iLAI*sum(Pso.*piLs + (Po-Pso).*piLd);  % Total canopy pi * F in obs dir for this PS
        LoF_(i,PS)      = piLtoti/pi;
        
        piLem           = Ps.*piLs + (1-Ps).*piLd; % CvT added 18 Nov 2013
        piLem_(i,:,PS)  = piLem(:);       

        t2   = xdd2 * (rho(i)-tau(i))/2;
        att  = 1-(rho(i)+tau(i))/2+t2;
        sig  = (rho(i)+tau(i))/2+t2;
        m    = sqrt(att^2-sig^2);
        rinf = (att - m)/sig;
        fac  = 1 - m * iLAI;
        facs = (rs(i)-rinf)/(1-rs(i)*rinf);
             
        G1(1) = 2; Gnew = 0;                 % (to make sure we will enter the loop)
        while abs(Gnew-G1(1)) > 1e-3
            G1(1) = Gnew;
            for j=2:nl+1
                G1(j) = fac * G1(j-1) + (Fmin(i,j)+rinf*Fplu(i,j-1)) * iLAI;
            end
            G2(nl+1) = G1(nl+1) * facs;
            for j=nl:-1:1
                G2(j) = fac * G2(j+1) + (rinf*Fmin(i,j+1)+Fplu(i,j)) * iLAI;
            end
            Gnew = -rinf * G2(1);
        end
        
        Fhem_(i,PS) = (rinf*G1(1)+G2(1))/(1-rinf^2);
    end
    for ilayer = 1:nl+1
        Fiprofile(ilayer,PS) = 0.001 * Sint(Fplu(:,ilayer),spectral.wlF);
    end
end
rad.LoF_  = LoF_(:,1)  + LoF_(:,2);
rad.LoF1_ = LoF_(:,1);
rad.LoF2_ = LoF_(:,2);
rad.Fhem_ = Fhem_(:,1) + Fhem_(:,2);

profiles.fluorescence   = Fiprofile(:,1) + Fiprofile(:,2);
profiles.fluorescenceEm = zeros(nl+1,1);
for i = 1:nl+1
    profiles.fluorescenceEm(i) = 0.001 * Sint(sum(piLem_(:,i,:),3),spectral.wlF)';
end

rad.Eoutf = 0.001 * Sint(sum(Fhem_,2),spectral.wlF); 
rad.Eminf_ = Emin_;
rad.Epluf_ = Eplu_;
%ipl

%Esun   = pi.*(rad.T4+rad.LoF_);
%Esky   = pi./(1-t3.*rad.rdd).*(rad.T5+rad.T4.*rad.t3.*rad.rsd+rad.Fhem_.*rad.t3);

%rad.L_BOA  =(rad.rso.*Esun+rad.rdo.*Esky)./pi; 
%rad.L_adj = ((rad.T5.*rad.rdd+rad.T4.*rad.rsd)+rad.Fhem_).*rad.t7./(1-rad.rdd.*rad.t3);
%rad.L_TOA = rad.Lp0 + L_BOA.*rad.t6+ L_adj.*rad.t7;


