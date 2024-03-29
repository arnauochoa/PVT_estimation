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
%--     Number of constellations used
Nconst      =   length(const);
%--     Navigation RINEX files
NavFile     =   ["RINEX/BCLN00ESP_R_20182870000_01D_GN.rnx", "RINEX/BCLN00ESP_R_20182870000_01D_EN.rnx"];
%--     Observation RINEX file
ObsFile     =   "RINEX/BCLN00ESP_R_20182870000_01D_30S_MO.rnx";
%--     Number of epochs to be analyzed (max. 2880)
Nepoch      =   1000;                
%--     Number of unknowns of the PVT solution
Nsol        =   5;                  
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
Tcorr       =   zeros(Nepoch, 1);
%--     Propagation corrections mean for every epoch
Pcorr       =   zeros(Nepoch, 1);


%-  Initialization Parameters
eph     =   cell(1, Nconst);
iono    =   cell(1, Nconst);
for con = 1:Nconst
    %--     Get ephemerides and iono information from navigation message
    [eph{con}, iono{con}] =  getNavRINEX(NavFile(con));
end

%
PVT         =   nan(Nepoch,Nsol);        %   PVT solution for every constellation
GDOP        =   zeros(Nepoch, 1);           %   Gdop
PDOP        =   zeros(Nepoch, 1);           %   Pdop
TOW         =   nan(Nepoch,1);             %   Time Of the Week (TOW)
% G           =   cell(1, Nepoch);              %   Array of geometry matrixes as cells
pos_llh     =   nan(Nepoch, 3, Nconst);         %   Position in Latitude, Longitude and Height

%- Obtain pr and sats vectors for each constellation and epoch
pr_arr      = cell(Nepoch, 1);
sats_arr    = cell(Nepoch, 1);
const_arr   = cell(Nepoch, 1);
for con = 1:Nconst
    fid     =   fopen(ObsFile);
    [Nobs, Obs_types, year, Rin_vers]  =   anheader(fid);
    for  epoch = 1:Nepoch
        [pr, TOW(epoch), sats]  = getPR_epoch0(fid, year, Obs_types, Nobs, Rin_vers, const(con));
        Nsat(epoch, con)        = length(sats);
        
        pr_arr{epoch}           = [pr_arr{epoch}, pr];
        sats_arr{epoch}         = [sats_arr{epoch}, sats];
        const_arr{epoch}        = [const_arr{epoch}, con * ones(1, Nsat(epoch, con))];
    end
    fclose(fid);
end

for  epoch = 1:Nepoch
    %--     Compute the PVT solution at the next epoch
        [PVT(epoch, :), A, Tcorr(epoch), Pcorr(epoch)]  = ...
            MC_PVT_recLS(pr_arr{epoch}, sats_arr{epoch}, TOW(epoch), eph, ...
            iono, Nit, PVT0, const_arr{epoch}, const, enab_corr);
        
        G           = inv(A'*A);      % Geometry matrix computation

        G_diag= diag(G);
        GDOP(epoch) = sqrt(sum(G_diag));
        PDOP(epoch) = sqrt(G_diag(1) + G_diag(2) + G_diag(3));
        
        PVT0 = PVT(epoch, :);
end

%
%
%-  Show results

fprintf(' ==== RESULTS ==== \n')
Nmov            =   5;

pos_mean        =   nanmean(PVT(:,1:3), 1);
posllh_mean     =   rad2deg(xyz2llh(pos_mean));
t_err           =   PVT(:,4)/c;
t_err_mean      =   mean(t_err);
t_dif_err       =   PVT(:,5)/c;
t_dif_err_mean  =   mean(t_dif_err);
mu_mov          =   movmean(PVT(:,1:3),[Nmov-1 0], 1);
spread          =   nanstd(PVT(:,1:3), 0, 1);
p_err           =   sqrt((PVTr(1:3) - PVT(:, 1:3)).^2);
p_err_mov       =   PVTr(1:3) - mu_mov;
rms             =   ((PVTr(1) - PVT(:,1)).^2 + (PVTr(2) - PVT(:,2)).^2 + (PVTr(3) - PVT(:,3)).^2).^0.5;

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
    title(sprintf('Evolution of the errors in the different axes (Multiconst)\n'));
    filename = sprintf('Capt/mult_2/err_XYZ_%u.jpg', Nepoch);
    saveas(fig, filename);
    %
    % -- RMS evolution 
    fig = figure('DefaultAxesFontSize', 12); plot(TOW, rms);
    xlabel('Time of the Week (s)');
    ylabel('Root Mean Square error (m)');
    title(sprintf('RMS evolution (Multiconst)\n'));
    filename = sprintf('Capt/mult_2/rms_%u_%u.jpg', Nepoch, enab_corr);
    saveas(fig, filename);
    %
    % -- Error in time
    fig = figure('DefaultAxesFontSize', 12); plot(TOW, t_err);
    xlabel('Time of the Week (s)');
    ylabel('Error in time (s)');
    title(sprintf('Evolution of the bias in time (Multiconst)\n'));
    filename = sprintf('Capt/mult_2/tbias_%u_%u.jpg', Nepoch, enab_corr);
    saveas(fig, filename);
    %
    % -- Time difference between constellations
    fig = figure('DefaultAxesFontSize', 12); plot(TOW, t_dif_err);
    xlabel('Time of the Week (s)');
    ylabel('Error in time (s)');
    title(sprintf('Evolution of the time difference between %s and %s (Multiconst)\n', const(1), const(2)));
    filename = sprintf('Capt/mult_2/t_dif_bias_%u_%u.jpg', Nepoch, enab_corr);
    saveas(fig, filename);
    %
    % -- # satellites in view
    fig = figure('DefaultAxesFontSize', 12); plot(TOW, Nsat);
    legend(const);
    xlabel('Time of the Week (s)');
    ylabel('Number of satellites');
    title(sprintf('Evolution of number of satellites in view (Multiconst)\n'));
    filename = sprintf('Capt/mult_2/nsat_%u_%u.jpg', Nepoch, enab_corr);
    saveas(fig, filename);
    %
    % -- RMS evolution & # sal¡tellites
    fig = figure('DefaultAxesFontSize', 12);
    yyaxis left;
    plot(TOW, rms);
    ylabel('Root Mean Square error (m)');
    yyaxis right;
    plot(TOW, Nsat);
    legend(const);
    ylabel('Number of satellites');
    xlabel('Time of the Week (s)');
    title(sprintf('RMS evolution (Multiconst)\n'));
    filename = sprintf('Capt/mult_2/rms-nsat_%u._%u.jpg', Nepoch, enab_corr);
    saveas(fig, filename);
    %
    % -- GDOP evolution & # sal¡tellites
    fig = figure('DefaultAxesFontSize', 12); 
    yyaxis left;
    plot(TOW, GDOP);
    ylabel('GDOP');
    yyaxis right;
    plot(TOW, Nsat);
    ylabel('Number of satellites');
    xlabel('Time of the Week (s)');
    title(sprintf('GDOP evolution (Multiconst)\n'));
    filename = sprintf('Capt/mult_2/gdop-nsat_%u_%u.jpg', Nepoch, enab_corr);
    saveas(fig, filename);
    %
    % -- PDOP evolution & # sal¡tellites
    fig = figure('DefaultAxesFontSize', 12); 
    yyaxis left;
    plot(TOW, PDOP);
    ylabel('PDOP');
    yyaxis right;
    plot(TOW, Nsat);
    ylabel('Number of satellites');
    xlabel('Time of the Week (s)');
    title(sprintf('PDOP evolution (Multiconst)\n'));
    filename = sprintf('Capt/mult_2/pdop-nsat_%u_%u.jpg', Nepoch, enab_corr);
    saveas(fig, filename);
    %
    % -- Time & Propagation corrections evolution
    fig = figure('DefaultAxesFontSize', 12);
    yyaxis left;
    plot(TOW, Tcorr);
    ylabel('Time corrections (s)');
    yyaxis right;
    plot(TOW, Pcorr);
    ylabel('Propagation corrections (m)');
    xlabel('Time of the Week (s)');
    title(sprintf('Time & Propagation corrections evolution (Multiconst)\n'));
    filename = sprintf('Capt/mult_2/tcorr-pcorr_%u_%u.jpg', Nepoch, enab_corr);
    saveas(fig, filename);
end
    
