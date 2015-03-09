function [soil,leafbio,canopy,meteo,angles,xyt] = select_input(V,vi,canopy,options,xyt,soil)

soil.spectrum= V(15).Val(vi(15));
soil.rss = V(16).Val(vi(16));
soil.rs_thermal = V(17).Val(vi(17));
soil.cs = V(18).Val(vi(18));
soil.rhos = V(19).Val(vi(19));
soil.CSSOIL  = V(42).Val(vi(42));
soil.lambdas = V(20).Val(vi(20));
soil.rbs  = V(43).Val(vi(43));
soil.SMC  = V(53).Val(vi(53));

leafbio.Cab = V(1).Val(vi(1));
leafbio.Cdm  = V(2).Val(vi(2));
leafbio.Cw = V(3).Val(vi(3));
leafbio.Cs   = V(4).Val(vi(4));
leafbio.N = V(5).Val(vi(5));
leafbio.Vcmo  = V(8).Val(vi(8));
leafbio.m  = V(9).Val(vi(9));
leafbio.Type = V(10).Val(vi(10));
leafbio.Tparam  = V(13).Val(:); % this is correct (: instead of 13)
leafbio.fqe(2)   = V(14).Val(vi(14));
leafbio.Rdparam  = V(12).Val(vi(12));

leafbio.rho_thermal   = V(6).Val(vi(6));
leafbio.tau_thermal  = V(7).Val(vi(7));

leafbio.Tyear = V(54).Val(vi(54));
leafbio.beta   = V(55).Val(vi(55));
leafbio.kNPQs   = V(56).Val(vi(56));
leafbio.qLs  = V(57).Val(vi(57));
leafbio.stressfactor  = V(58).Val(vi(58));


canopy.LAI  = V(21).Val(vi(21));
canopy.hc  = V(22).Val(vi(22));
canopy.LIDFa = V(25).Val(vi(25));
canopy.LIDFb  = V(26).Val(vi(26)); % this is correct (25 instead of 26)
canopy.leafwidth  = V(27).Val(vi(27));
canopy.rb   = V(37).Val(vi(37));
canopy.Cd  = V(38).Val(vi(38));
canopy.CR = V(39).Val(vi(39));
canopy.CD1  = V(40).Val(vi(40));
canopy.Psicor  = V(41).Val(vi(41));
canopy.rwc  = V(44).Val(vi(44));
canopy.kV = V(11).Val(vi(11));

meteo.zo  = V(23).Val(vi(23));
meteo.d   = V(24).Val(vi(24));
meteo.z  = V(28).Val(vi(28));
meteo.Rin   = V(29).Val(vi(29));
meteo.Ta = V(30).Val(vi(30));
meteo.Rli  = V(31).Val(vi(31));
meteo.p  = V(32).Val(vi(32));
meteo.ea  = V(33).Val(vi(33));
meteo.u   = V(34).Val(vi(34));
meteo.Ca = V(35).Val(vi(35));
meteo.Oa  = V(36).Val(vi(36));

xyt.startDOY = V(45).Val(vi(45));
xyt.endDOY = V(46).Val(vi(46));
xyt.LAT = V(47).Val(vi(47));
xyt.LON = V(48).Val(vi(48));
xyt.timezn = V(49).Val(vi(49));

angles.tts = V(50).Val(vi(50));
angles.tto = V(51).Val(vi(51));
angles.psi = V(52).Val(vi(52));

%% derived input
if options.soil_heat_method ==1
    soil.GAM =  Soil_Inertia1(soil.SMC);
else
    soil.GAM  = Soil_Inertia0(soil.cs,soil.rhos,soil.lambdas);
end
if options.calc_rss_rbs
    [soil.rss,soil.rbs] = calc_rssrbs(soil.SMC,canopy.LAI,soil.rbs);
end

if leafbio.Type,
    leafbio.Type = 'C4';
else
    leafbio.Type = 'C3';
end
canopy.hot  = canopy.leafwidth/canopy.hc;
if options.calc_zo
    [canopy.zo,canopy.d ]  = zo_and_d(soil,canopy);
end
leafbio.fqe(1)   = leafbio.fqe(2)/5;

