function    [Obs,svn]     =   getObsRINEX(fid, Nsvn, Nobs, vers, const, colc)
% getObs:   Get the observations of the selected epoch in a RINEX file.
%           Reads observations of Nsvn satellites. The selected epoch is
%           given by the position of the fid. It controls the RINEX version
%           given by vers = {2,3} for RINEX 2.XX or 3.XX, resp. For version
%           3, the function also outputs the SVN (string) in the given
%           epoch.
%
% Modified from Kai Borre 09-13-96

    global lin

    
    % Get number ob observation lines for each epoch/svn
    if( vers == 2)  % RINEX 2.XX (only 5 types of observation per line)
        Obs     =   nan(Nsvn,Nobs);
        Nlin    =   ceil(Nobs/5);
        svn     =   [];
    else            % RINEX 3.XX (All observation in the same line/svn)
        Nlin    =   1;
        switch const
            case 'GPS'
                strc    =   'G';            % Constellation letter for GPS
        end
    end
    %
    %-  Read line for each satellite
    for ii = 1:Nsvn
        if( vers == 2 )
            for jj = 1:Nlin
                lin     =   fgetl(fid);     % Get line jj for SVN ii
                init    =   (jj-1)*5;
                if( jj == Nlin )
                    Ntmp    =   Nobs - init;
                else
                    Ntmp    =   5;
                end
                %-  Store each observable (e.g. Phase, Code, SNR, ...) in a diferent
                %   column (kk+init) of SVN ii
                for kk = 1:Ntmp
                    Obs(ii,kk+init)     =   str2double(lin(1+16*(kk-1):16*kk-2));
                end
            end        
        else
            lin             =   fgetl(fid);     % Get line for SVN ii
            %
            if( strcmp(lin(1),strc) )
                svn(ii)     =   str2num(lin(2:3));       % Check format RINEX 3.XX
                %-  Store each observable (e.g. Phase, Code, SNR, ...) in a diferent
                %   column (kk+init) of SVN ii
                init        =   3+1+16*(colc-1);
                fin         =   16*colc-2+3;
                pr          =   str2num(lin(init:fin));
                if isempty(pr)
                    Obs(ii)     =   NaN;
                else
                    Obs(ii)     =   pr;
                end
            end
        end
    end

end