function    [eph,iono]  =   getNavRINEX(NavRINEXFile)
% get_eph:  Reads a RINEX Navigation file and reformats the data into a
%           matrix with 22 rows with the ephemerides information and a 
%           column for each satellite/epoch. The code is done for RINEX
%           2.XX, check properties of RINEX 3.XX.
%
% Input:            
%           NavRINEXFile:   File name of the Rinex Navigation file
%
% Output:   eph:            22xNeph matrix with the ephemerides information
%                           for the satellites contained in the RINEX file. 
%                           Units are either seconds, meters, or radians  
%           iono:           8x1 vector containing ionosphere parameters
%
% Modified from Kai Borre 04-18-96
% Copyright (c) by Kai Borre
% $Revision: 1.0 $  $Date: 1997/09/24  $

    fide        =   fopen(NavRINEXFile);
    head_lines  =   0;
    iono        =   zeros(8,1);
    %
    %-  Read the header. Skip all the header (look for 'END OF HEADER' 
    %   string) except the ionosphere labels, which are read
    answer      =   [];
    while(isempty(answer))  
       head_lines   =   head_lines+1;   %   Count the number of header lines
       line         =   fgetl(fide);    %   Read next line
       %
       %    Version flag
       vers_found   =   ~isempty(strfind(line,'RINEX VERSION / TYPE'));
       if(vers_found)
           version  =   str2double(line(1:9));   %   Get RINEX version
       end
       %    Iono info flag
       iono_found   =   (~isempty(strfind(line,'ION ALPHA')) ||...
                        ~isempty(strfind(line,'IONOSPHERIC CORR')));
       if(iono_found)   %   If the ionosphere parameters label was found
           %--  Save the 8 ionosphere parameters
           %
           %---     Get first 4 parameters
           if( version >= 3 )
               data     =   textscan(line(5:end),'%f%f%f%f%*[^\n]');
               id       =   1;
           else
               data     =   textscan(line(1:end),'%f%f%f%f%*[^\n]');
               id       =   0;
           end
           const        =   line(1:3);  % Get constellation
           %
            if(strcmp(const,'GPS'))     % Get iono param for GPS
                iono(1) = data{1};
                iono(2) = data{2};
                iono(3) = data{3};
                iono(4) = data{4};
                line    = [];
                while isempty(line)
                    line = fgetl(fide);
                    head_lines  =   head_lines + 1;
                end
                data = textscan(line(5:end),'%f%f%f%f%*[^\n]');
                if ~isempty(data(4))
                    iono(5) =   data{1};
                    iono(6) =   data{2};
                    iono(7) =   data{3};
                    iono(8) =   data{4};
                else
                    iono    =   zeros(8,1);
                end
            else                        % We do not use iono param gor GAL
                iono        =   [];
            end
       end
       answer       =   strfind(line,'END OF HEADER');
    end
    %
    %-  Get number of ephemerides/satellites (Read the whole file and count the
    %   lines, the move the file pointer to the begining of the file)
    noeph   =   0;
    line    =   fgetl(fide);
    while(line ~= -1)
       noeph = noeph+1;
       line = fgetl(fide);
    %    if line == -1, break;  end
    end;
    noeph   =   noeph/8;    % Each ephemerides contains 8 lines
    frewind(fide);          % Come back to the begining of the file
    % Skip the header again
    for ii = 1:head_lines 
        line = fgetl(fide); 
    end
    %
    %-  Set aside memory for the input
    svprn       =   zeros(1,noeph);
%     weekno      =   zeros(1,noeph);
%     t0c         =   zeros(1,noeph);
    tgd         =   zeros(1,noeph);
%     aodc        =   zeros(1,noeph);
%     toe         =   zeros(1,noeph);
    af2         =   zeros(1,noeph);
    af1         =   zeros(1,noeph);
    af0         =   zeros(1,noeph);
%     aode        =   zeros(1,noeph);
    deltan      =   zeros(1,noeph);
    M0          =   zeros(1,noeph);
    ecc         =   zeros(1,noeph);
    roota       =   zeros(1,noeph);
    toe         =   zeros(1,noeph);
    cic         =   zeros(1,noeph);
    crc         =   zeros(1,noeph);
    cis         =   zeros(1,noeph);
    crs         =   zeros(1,noeph);
    cuc         =   zeros(1,noeph);
    cus         =   zeros(1,noeph);
    Omega0      =   zeros(1,noeph);
    omega       =   zeros(1,noeph);
    i0          =   zeros(1,noeph);
    Omegadot    =   zeros(1,noeph);
    idot        =   zeros(1,noeph);
%     accuracy    =   zeros(1,noeph);
%     health      =   zeros(1,noeph);
%     fit         =   zeros(1,noeph);
    %
    %- Get the data taking into account the RINEX format (from Kai Borre)
    for ii = 1:noeph
        line         =   fgetl(fide);               %   1st line
        svprn(ii)    =   str2num(line(id+(1:2)));
%        year         =   line(3:6);
%        month        =   line(7:9);
%        day          =   line(10:12);
%        hour         =   line(13:15);
%        minute       =   line(16:18);
%        second       =   line(19:22);
        af0(ii)    = str2num(line(id+(23:41)));
        af1(ii)    = str2num(line(id+(42:60)));
        af2(ii)    = str2num(line(id+(61:79))); 
       %
        line         =   fgetl(fide);               %   2nd line  
%        IODE         =   line(4:22);
        crs(ii)    = str2num(line(id+(23:41)));
        deltan(ii) = str2num(line(id+(42:60)));
        M0(ii)     = str2num(line(id+(61:79)));
       %
        line         =   fgetl(fide);               %   3rd line
        cuc(ii)    = str2num(line(id+(4:22)));
        ecc(ii)    = str2num(line(id+(23:41)));
        cus(ii)    = str2num(line(id+(42:60)));
        roota(ii)  = str2num(line(id+(61:79)));
       %
       line         =   fgetl(fide);                %   4th line
       toe(ii)      =   str2num(line(id+(4:22)));
        if (strcmp(const,'GAL') && toe(ii) > 2500) 
            toe(ii) = toe(ii) - 1024; 
        end
       cic(ii)      =   str2num(line(id+(23:41)));
       Omega0(ii)   =   str2num(line(id+(42:60)));
       cis(ii)      =   str2num(line(id+(61:79)));
       %
       line         =   fgetl(fide);                %   5th line   
       i0(ii)       =   str2num(line(id+(4:22)));
       crc(ii)      =   str2num(line(id+(23:41)));
       omega(ii)    =   str2num(line(id+(42:60)));
       Omegadot(ii) =   str2num(line(id+(61:79)));
       %
       line         =   fgetl(fide);                %   6th line 
       idot(ii)     =   str2num(line(id+(4:22)));
%        codes        =   str2num(line(23:41));
%        weekno       =   str2num(line(id+42:60));
%        L2flag       =   str2num(line(61:79));
       %
       line         =   fgetl(fide);                %   7th line 
%        svaccur      =   str2num(line(4:22));
%        svhealth     =   str2num(line(23:41));
       tgd(ii)      =   str2num(line(id+(42:60)));
%        iodc         =   line(61:79);
       %
       line         =   fgetl(fide);                %   8th line 
%        tom(ii)      =   str2num(line(id+4:22));
    end
    status = fclose(fide);
    %
    %-  Description of variable eph.
    eph(1,:)    =   svprn;
    eph(2,:)    =   af2;
    eph(3,:)    =   M0;
    eph(4,:)    =   roota;
    eph(5,:)    =   deltan;
    eph(6,:)    =   ecc;
    eph(7,:)    =   omega;
    eph(8,:)    =   cuc;
    eph(9,:)    =   cus;
    eph(10,:)   =   crc;
    eph(11,:)   =   crs;
    eph(12,:)   =   i0;
    eph(13,:)   =   idot;
    eph(14,:)   =   cic;
    eph(15,:)   =   cis;
    eph(16,:)   =   Omega0;
    eph(17,:)   =   Omegadot;
    eph(18,:)   =   toe;
    eph(19,:)   =   af0;
    eph(20,:)   =   af1;
    eph(21,:)   =   toe;
    eph(22,:)   =   tgd;

end
%
%%%%%%%%% end get_eph.m %%%%%%%%%


















