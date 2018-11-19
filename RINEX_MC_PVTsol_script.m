% Analyze RINEX files in order to compute the PVT solution
%
clearvars;
close all;
clc;
fclose('all');
%
addpath 'Corrections';
addpath 'Corrections/Control_segment';
addpath 'Corrections/Prop_Effects';
addpath 'Ephemeris';
addpath 'Misc';
addpath 'Observations';

%-  Setting Parameters
%--     Satellite constellations to be used for the PVT computation
const       =   ["GPS", "GAL"];
%--     Enable/disable corrections
enab_corr   =   true;
%--     Weight given to the PVT solution of every constellation
weight      =   [0.5, 0.5];
%--     Number of constellations used
Nconst      =   length(const);
%--     Navigation RINEX files
NavFile     =   ["RINEX/BCLN00ESP_R_20182870000_01D_GN.rnx", "RINEX/BCLN00ESP_R_20182870000_01D_EN.rnx"];
%--     Observation RINEX file
ObsFile     =   "RINEX/BCLN00ESP_R_20182870000_01D_30S_MO.rnx";
%--     Number of epochs to be analyzed (max. 2880)
Nepoch      =   1;                
%--     Number of unknowns of the PVT solution
Nsol        =   4;                  
%--     Number of iterations used to obtain the PVT solution
Nit         =   8;                   
%--     Reference position (check RINEX file or website of the station)
PVTr        =   [4788065.1430, 167551.1700, 4196354.9920];   %FIXME: add reference time
%--     Preliminary guess for PVT solution 
PVT0        =   [0 0 0 0];
%--     Speed of light (for error calculations)
c           =   299792458;       %   Speed of light (m/s)
%--     Number of satellites for every epoch
Nsat        =   zeros(Nepoch, Nconst);
%--     Time corrections mean for every epoch
Tcorr       =   zeros(Nepoch, Nconst);
%--     Propagation corrections mean for every epoch
Pcorr       =   zeros(Nepoch, Nconst);


%-  Initialization Parameters
eph     =   cell(1, Nconst);
iono    =   cell(1, Nconst);
for con = 1:Nconst
    %--     Get ephemerides and iono information from navigation message
    [eph{con}, iono{con}] =  getNavRINEX(NavFile(con));
end

%
PVT         =   nan(Nepoch,Nsol,Nconst);        %   PVT solution for every constellation
PVTf        =   zeros(Nepoch, Nsol);            %   Final PVT solution
GDOP        =   zeros(Nepoch,Nconst);           %   Gdop
PDOP        =   zeros(Nepoch,Nconst);           %   Pdop
TOW         =   nan(Nepoch,1);             %   Time Of the Week (TOW)
% G           =   cell(1, Nepoch);              %   Array of geometry matrixes as cells
pos_llh     =   nan(Nepoch, 3, Nconst);         %   Position in Latitude, Longitude and Height

%- Obtain the PVT for every constellation
for con = 1:Nconst
        %--     Open Observation RINEX file and take the header information
    %   Nobs:       # of observables (integer), check RINEX version
    %   Obs_types:  List of observation types (string), check RINEX version
    %   Rin_vers:   RINEX version (integer: 2 or 3)
    fid         =   fopen(ObsFile);
    [Nobs, Obs_types, year, Rin_vers]  =   anheader(fid);

    %-  Sequentially read the Observation file and compute the PVT solution
    for  epoch = 1:Nepoch
        %--     Get the observed pseudoranges for all satellites in view of 
        %       a given constellation at the next epoch
        %
        %   pr:     Pseudoranges at given TOW for satellites in sats
        %           (Nsatx1)
        %   TOW:    Time Of the Week (TOW) of the next epoch    
        %   sats:   Satellites in view  
        %
        %--  
        [pr, TOW(epoch), sats] = getPR_epoch0(fid, year, Obs_types, Nobs, Rin_vers, const(con));
        Nsat(epoch, con) = length(sats);
    
        %--     Compute the PVT solution at the next epoch
        [PVT(epoch, :, con), A, Tcorr(epoch, con), Pcorr(epoch, con)]  = ...  
            PVT_recLS(pr, sats, TOW(epoch), eph{con}, iono{con}, Nit, PVT0, enab_corr);

        G           = inv(A'*A);      % Geometry matrix computation

        G_diag= diag(G);
        GDOP(epoch) = sqrt(sum(G_diag));
        PDOP(epoch) = sqrt(G_diag(1) + G_diag(2) + G_diag(3));
        
        PVT0 = PVT(epoch, :, con);
    end
    PVTf(:, :)   = PVTf(:, :) + weight(con) .* PVT(:, :, con);
    fclose(fid);
end

%
%
%-  Show results

fprintf(' ==== RESULTS ==== \n')
Nmov            =   5;

pos_mean        =   nanmean(PVTf(:,1:3),1);
posllh_mean     =   rad2deg(xyz2llh(pos_mean));
t_err           =   PVTf(:,4)/c;
t_err_mean      =   mean(t_err);
mu_mov          =   movmean(PVTf(:,1:3),[Nmov-1 0],1);
spread          =   nanstd(PVTf(:,1:3), 0, 1);
p_err           =   sqrt((PVTr(1:3) - PVTf(:, 1:3)).^2);
p_err_mov       =   PVTr(1:3) - mu_mov;
rms             =   ((PVTr(1) - PVTf(:,1)).^2 + (PVTr(2) - PVTf(:,2)).^2 + (PVTr(3) - PVTf(:,3)).^2).^0.5;

% -------------------------------------------------------------------------
fprintf('\nMean position as computed from %2.0f epochs:', Nepoch);
fprintf('\n\nX: %12.3f  Y: %12.3f  Z: %12.3f\n', pos_mean(1), pos_mean(2), pos_mean(3));
fprintf('\nLat.: %12.3fº   Long.: %12.3fº   Height: %12.3f m\n', posllh_mean(1), posllh_mean(2), posllh_mean(3));
fprintf('\nMean time error as computed from %2.0f epochs: %G seconds\n', Nepoch, t_err_mean);
fprintf('\nstd (m) of each position coordinate:');
fprintf('\nX: %2.2f Y: %2.2f Z:%2.2f\n', spread(1), spread(2), spread(3));
if Nepoch==1    % Print error position when there's only one epoch
    fprintf('\nError in position as computed from 1 epoch:');
    fprintf('\nX: %2.2f m Y: %2.2f m Z:%2.2f m \n', p_err(1), p_err(2), p_err(3));
end
% -------------------------------------------------------------------------

if Nepoch > 1   % Plot evolutions for all epochs
    % -- Errors (RMS) in X-Y-Z
    fig = figure('DefaultAxesFontSize', 12); plot(TOW, p_err);
    legend('X error', 'Y error', 'Z error');
    xlabel('Time of the Week (s)');
    ylabel('Error for each coordinate (m)');
    title(sprintf('Evolution of the errors in the different axes (GPS weight = %1.2f, GAL weight = %1.2f)\n', weight(1), weight(2)));
    filename = sprintf('Capt/mult_1/err_XYZ_%u_%u_%u.jpg', 100*weight(1), 100*weight(2), enab_corr);
    saveas(fig, filename);
    %
    % -- RMS evolution 
    fig = figure('DefaultAxesFontSize', 12); plot(TOW, rms);
    xlabel('Time of the Week (s)');
    ylabel('Root Mean Square error (m)');
    title(sprintf('RMS evolution (GPS weight = %1.2f, GAL weight = %1.2f)\n', weight(1), weight(2)));
    filename = sprintf('Capt/mult_1/rms_%u_%u_%u.jpg', 100*weight(1), 100*weight(2), enab_corr);
    saveas(fig, filename);
    %
    % -- Error in time
    figure('DefaultAxesFontSize', 12); plot(TOW, t_err);
    xlabel('Time of the Week (s)');
    ylabel('Error in time (s)');
    title(sprintf('Evolution of the bias in time (GPS weight = %1.2f, GAL weight = %1.2f)\n', weight(1), weight(2)));
    filename = sprintf('Capt/mult_1/tbias_%u_%u_%u.jpg', 100*weight(1), 100*weight(2), enab_corr);
    saveas(fig, filename);
    %
    % -- # satellites in view
    figure('DefaultAxesFontSize', 12); plot(TOW, Nsat);
    legend(const);
    xlabel('Time of the Week (s)');
    ylabel('Number of satellites');
    title(sprintf('Evolution of number of satellites in view (GPS weight = %1.2f, GAL weight = %1.2f)\n', weight(1), weight(2)));
    filename = sprintf('Capt/mult_1/nsat_%u_%u_%u.jpg', 100*weight(1), 100*weight(2), enab_corr);
    saveas(fig, filename);
    %
    % -- RMS evolution & # sal¡tellites
    figure('DefaultAxesFontSize', 12);
    yyaxis left;
    plot(TOW, rms);
    ylabel('Root Mean Square error (m)');
    yyaxis right;
    plot(TOW, Nsat);
    legend(const);
    ylabel('Number of satellites');
    xlabel('Time of the Week (s)');
    title(sprintf('RMS evolution (GPS weight = %1.2f, GAL weight = %1.2f)\n', weight(1), weight(2)));
    filename = sprintf('Capt/mult_1/rms-nsat_%u_%u_%u.jpg', 100*weight(1), 100*weight(2), enab_corr);
    saveas(fig, filename);
    %
    % -- GDOP evolution & # sal¡tellites
    figure('DefaultAxesFontSize', 12); 
    yyaxis left;
    plot(TOW, GDOP);
    legend(const);
    ylabel('GDOP');
    yyaxis right;
    plot(TOW, Nsat);
    legend(const);
    ylabel('Number of satellites');
    xlabel('Time of the Week (s)');
    title(sprintf('GDOP evolution (GPS weight = %1.2f, GAL weight = %1.2f)\n', weight(1), weight(2)));
    filename = sprintf('Capt/mult_1/gdop-nsat_%u_%u_%u.jpg', 100*weight(1), 100*weight(2), enab_corr);
    saveas(fig, filename);
    %
    % -- PDOP evolution & # sal¡tellites
    figure('DefaultAxesFontSize', 12); 
    yyaxis left;
    plot(TOW, PDOP);
    legend(const);
    ylabel('PDOP');
    yyaxis right;
    plot(TOW, Nsat);
    legend(const);
    ylabel('Number of satellites');
    xlabel('Time of the Week (s)');
    title(sprintf('PDOP evolution (GPS weight = %1.2f, GAL weight = %1.2f)\n', weight(1), weight(2)));
    filename = sprintf('Capt/mult_1/pdop-nsat_%u_%u._%u.jpg', 100*weight(1), 100*weight(2), enab_corr);
    saveas(fig, filename);
    %
    % -- Time & Propagation corrections evolution
    figure('DefaultAxesFontSize', 12);
    yyaxis left;
    plot(TOW, Tcorr);
    legend(const);
    ylabel('Time corrections (s)');
    yyaxis right;
    plot(TOW, Pcorr);
    legend(const);
    ylabel('Propagation corrections (m)');
    xlabel('Time of the Week (s)');
    title(sprintf('Time & Propagation corrections evolution (GPS weight = %1.2f, GAL weight = %1.2f)\n', weight(1), weight(2)));
    filename = sprintf('Capt/mult_1/tcorr-pcorr_%u_%u_%u.jpg', 100*weight(1), 100*weight(2), enab_corr);
    saveas(fig, filename);
end



