clear all; close all; clc;

% -- CONSTELLATION LETTERS --
%   G	GPS       GPS
%   R	GLONASS   GLO
%   E	Galileo   GAL

addpath('Observations/');
addpath('Misc/');

observationFile         = 'RINEX/BCLN00ESP_R_20182870000_01D_30S_MO.rnx';
navigationFileGalileo   = 'RINEX/BCLN00ESP_R_20182870000_01D_EN.rnx';
navigationFileGPS       = 'RINEX/BCLN00ESP_R_20182870000_01D_GN.rnx';

constellations   =  ["GPS", "GLO", "GAL"];
Nepoch           =  10;       

[ephEN, ionoEN]  =  getNavRINEX(navigationFileGalileo);
[ephGN, ionoGN]  =  getNavRINEX(navigationFileGPS);

fid              = fopen(observationFile);
[Nobs, Obs_types, year, Rin_vers]  =   anheader(fid);

for epoch = 1:1:Nepoch
    [GPS.pr{epoch}, GPS.TOW{epoch}, GPS.sats{epoch}] = ...
        getPR_epoch0(fid, year, Obs_types, Nobs, Rin_vers, constellations(1));
    
%     [Obs, svn] = getObsRINEX(fid, length(GPS.sats{epoch}), Nobs, Rin_vers, constellations(1), colc);
end