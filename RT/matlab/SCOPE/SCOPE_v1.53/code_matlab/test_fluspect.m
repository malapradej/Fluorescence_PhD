path_input      = '../../data/input/';     
%% characteristic spectra for chemicals (input for FLUSPECT)
opticoef    = load([path_input,'fluspect_parameters/','optipar_fluspect.txt']);  % file with leaf spectral parameters

% Optical coefficient data used by fluspect
optipar.nr    = opticoef(:,2);
optipar.Kdm   = opticoef(:,3);
optipar.Kab   = opticoef(:,4);
optipar.Kw    = opticoef(:,5);
optipar.Ks    = opticoef(:,6);
optipar.phiI  = opticoef(:,7);
optipar.phiII = opticoef(:,8);

%% spectrals settings
[spectral] = define_bands;

%% leaf chemistry
leafbio.Cab         = 70;
leafbio.Cw          = .001;
leafbio.Cdm         = .1;
leafbio.Cs          = .0;
leafbio.N           = 1.4;
leafbio.fqe         = [.02/5 .02];



[leafopt] = fluspect(spectral,leafbio,optipar);       
