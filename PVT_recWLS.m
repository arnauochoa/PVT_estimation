function    [PVT, A, tcorr, Pcorr, X, mask_sats, sats_el]  =   PVT_recWLS(pr, svn, TOW, eph, iono, Nit, PVT0, enab_corr, threshold)
% PVT_recLS:    Computation of the receiver position at time TOW from  
%               pseudoranges (pr) and ephemerides information (eph). 
%               Implementation using the iterative Least-Squares principle 
%               after linearization of the navigation equations. We use the
%               svn Id (svn) to get the satellite coordinates corresponding 
%               to the pseudoranges given by the input (pr).
%               
% Input:            
%           pr:         1xNsat vector with the observed pseudoranges at time 
%                       TOW for the Nsat satellites with ID given by svn
%           svn:        1xNsat vector containing the ID of the satellites from
%                       which the pseudoranges in pr has been measured.
%           TOW:        TOW of the epoch that the pseudoranges in pr has been 
%                       measured
%           eph:        Matrix with the ephemerides information needed to
%                       obtain pseudorange corrections and the coordinates of
%                       the satellites
%           iono:       Ionosphere parameters to get ionosphere corrections
%           Nit:        # of iterations applied in the iterative LS algorithm
%                       used to obtain the PVT solution
%           PVT0:       Initial guess of the user position used for the
%                       linearization process of the navigaiton equations
%           enab_corr:  Flag for enable/disable corrections
%           threshold:  Threshold angle for rejecting satellites
%
% Output:   PVT:        Nsolx1 vector with the estimated PVT solution 
%           A:          NsatxNsol matrix with the geometry information. This is
%                       the geometry matrix studied in lectures: "PVT = A*pr"
%           tcorr:      Satellite clock correction in [sec.]
%           Pcorr:      Corrections in [meters] for the propagation effects. 
%                       Mainly ionosphere and troposphere corrections
%                   
%
    %-  Initialize parameters
    c       =   299792458;       %   Speed of light (m/s)
    Nsat    =   length(pr);      %   Number of satellites
    tcorr   =   zeros(Nsat, 1);  %   Satellite clock corrections    
    Pcorr   =   zeros(Nsat, 1);  %   Propagation effects corrections  
    X       =   zeros(3, Nsat);  %   Satellite coordinates
    sats_el =   zeros(Nsat, 1);  %   Satellites elevations (rad)
    
    %-  Reject satellites under the threshold elevation angle
    mask    =   [];
    for sat = 1:Nsat
        [X(:, sat), tcorr(sat)]  =   getCtrl_corr(eph, svn(sat), TOW, pr(sat));
        alpha           =   elevation_angle(PVT0(1:3)', X(:, sat));
        sats_el(sat)    = alpha;
        if alpha < threshold, mask = [mask, sat]; end
    end
    mask_sats   =   length(mask);
    pr(mask)    =   [];             % Deletes elemets of pr given by the indexes in mask
    svn(mask)   =   [];             % Deletes elemets of svn given by the indexes in mask
    Nsat        =   length(pr);     %   Number of satellites update
    
    %
    %-  Iterative LS to compute the PVT solution
    for iter = 1:Nit
        %-- Get initial position at epoch given by TOW (only 1st iteration)
        if (iter == 1)
            A       =   zeros(Nsat, 4);
            p       =   zeros(Nsat, 1);
            PVT     =   PVT0;
        end
        %   Fill the measurement vector and geometry matrix with the
        %   information of each satellite
        for sat = 1:Nsat
            %--     Get satellite coordinates and corrections for control 
            %       segment effects (e.g. satellite coordinates and satellite 
            %       clock offset)       
            %   X:     3x1 matrix with the cartesian coordinates of 
            %          the satellite with Id given by svn(sat) at time TOW
            %   tcorr: Clock correction (sec.) for the satellite with Id 
            %          given by svn(sat) at time TOW.      
            %
            if (iter == 1)  % Only 1st iteration
                [X(:, sat), tcorr(sat)]  =   getCtrl_corr(eph, svn(sat), TOW, pr(sat));
                % pr(sat) = pr(sat) + tcorr * c;  % Clock correction
            end
            %
            %--     Get corrections for propagation efects (e.g. iono, tropo,
            %       earth rotation, relativistic effects). Unit meters
            Pcorr     =   getProp_corr(X(:, sat), PVT0, iono, TOW);  % Iono + Tropo correction
            % TODO: falta añadir efecto relativista

            %--     Get corrected pseudorange (rho_c = rho - cor)
            if enab_corr         % If corrections are enabled
                corr          =   Pcorr + c * tcorr(sat);            %   Total correction factor in meters (TBD)
            else
                corr = 0;
            end
            pr_c          =   pr(sat) + corr;   %   Corrected pseudorange
            % pr_if = iono_corr(pr_c1, pr_c2);

            if (~isnan(pr_c))    % Fill as long as there is C1 measurement,
                                 % otherwise discard the measurement.
                %--     Fill the measurement vector (rhoc_i - d0_i)
                d0  = sqrt((X(1,sat) - PVT(1))^2 + (X(2,sat) - PVT(2))^2 + (X(3,sat) - PVT(3))^2); %TODO: shorten
                p(sat)       =   pr_c - d0;
                %
                %--     Fill the geometry matrix (TBD)
                ax          =   -(X(1, sat) - PVT(1)) / d0;   % dx/d0
                ay          =   -(X(2, sat) - PVT(2)) / d0;   % dy/d0
                az          =   -(X(3, sat) - PVT(3)) / d0;   % dz/d0
                A(sat, :)   =   [ax ay az 1];
                % 
            end
        end
        %--     Get the LS estimate of PVT at iteration iter
        d               =   pinv(A) * p;         % PVT update ("correction") !!!
        PVT(1:3)        =   PVT(1:3) + d(1:3)';  % Update the PVT coords.
        PVT(4)          =   d(4);        % Update receiver clock offset
        %
    end
    tcorr = mean(tcorr);
    Pcorr = mean(Pcorr);
end