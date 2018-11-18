function    [pr,TOW,sats]  =   getPR_epoch0(fido,year,Obs_types,Nobs,vers,const)
% getPR_epoch0:     Get the observed pseudoranges of the given epoch in the
%                   RINEX file (fido) taken at year (given by the header).
%                   The number of the observations in the RINEX file must
%                   be provided (Nobs) as well as the RINEX version. The
%                   function outputs the pseudorange (pr) at the given TOW
%                   of the epoch (TOW) for all the satellites in view
%                   (sats) of constellation (const).
%
% Output:   pr:     1xNsat vector with the pseudorange observable (C1) at
%                   the given epoch for the Nsat satellites in view of
%                   constellation (const).
%
%           TOW:    Time Of the Week (TOW). Time GPS in sec.
%           sats:   1xNsat vector with the integer identifiers for the Nsat
%                   satellites in view of constellation (const) at time TOW


    %--     Find the next epoch in the observation file and get the
    %       corresponding Time Of the Week (TOW), the # of total satellites
    %       (for all constellations)in the epoch (Nsat), and the list of 
    %       satellites from all constellations available (sats)
    %--
    if( vers == 2 )
        [TOW,Nsat,sats,~,~]     =   fepoch_0(fido,year);
    else
        [TOW,Nsat,~,~,~]           =   fepoch_0(fido,year);
    end

    %
    %--     Get the observed C1 pseudorange at the next epoch for all the
    %       satellites of constellation (const). 
    %--
    if( vers == 2 )
        [pr,~]                 =   getObs(fido,Nsat,Nobs,vers,const,Obs_types);
        %---    Find the corresponding satellites to the given
        %       constellation in (const)
        switch const
            case 'GPS'
                let     =   'G';
            case 'GLO'
                let     =   'R';
        end
        idx             =   strfind(sats,let);
        sats            =   str2double(sats((idx+1):(idx+2)));
        idx             =   (idx+2)/3;
        pr              =   pr(idx);
        %
    else
        %--     Get the observed C1 pseudorange (some column of Obs) 
        %---    Get number of observations corresponding to the desired 
        %       constelation 'const'
        Nobsc           =   Nobs.(const);   
%         [tmp,~]         =   obs_type_find(Obs_types,{const},'L1');
%         colc            =   tmp.(const).C1;
        [pr,sats]       =   getObs(fido,Nsat,Nobsc,vers,const,Obs_types);
    end
        
    %
    
end