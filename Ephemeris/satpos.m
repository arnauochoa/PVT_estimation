function [satp,satv]    =   satpos(t,eph)
%SATPOS   Calculation of X,Y,Z coordinates and velocity at time t
%         for given ephemeris eph

% Based on Kai Borre 04-09-96
% Copyright (c) by Kai Borre
% $Revision: 1.0 $  $Date: 1997/09/26  $

GM = 3.986008e14;             % earth's universal gravitational
% parameter m^3/s^2
Omegae_dot = 7.2921151467e-5; % earth rotation rate, rad/s

%  Units are either seconds, meters, or radians
%  Assigning the local variables to eph
M0      =   eph(3);
roota   =   eph(4);
deltan  =   eph(5);
ecc     =   eph(6);
omega   =   eph(7);
cuc     =   eph(8);
cus     =   eph(9);
crc     =  eph(10);
crs     =  eph(11);
i0      =  eph(12);
idot    =  eph(13);
cic     =  eph(14);
cis     =  eph(15);
Omega0  =  eph(16);
Omegadot=  eph(17);
toe     =  eph(18);

% Procedure for coordinate calculation
A = roota*roota;
tk = check_t(t-toe);
n0 = sqrt(GM/A^3);
n = n0+deltan;
M = M0+n*tk;
M = rem(M+2*pi,2*pi);
E = M;
for i = 1:10
   E_old = E;
   E = M+ecc*sin(E);
   dE = rem(E-E_old,2*pi);
   if abs(dE) < 1.e-12
      break;
   end
end
E = rem(E+2*pi,2*pi);
v = atan2(sqrt(1-ecc^2)*sin(E), cos(E)-ecc);
phi = v+omega;
phi = rem(phi,2*pi);
u = phi              + cuc*cos(2*phi)+cus*sin(2*phi);
r = A*(1-ecc*cos(E)) + crc*cos(2*phi)+crs*sin(2*phi);
i = i0+idot*tk       + cic*cos(2*phi)+cis*sin(2*phi);
Omega = Omega0+(Omegadot-Omegae_dot)*tk-Omegae_dot*toe;
Omega = rem(Omega+2*pi,2*pi);
x1 = cos(u)*r;
y1 = sin(u)*r;
satp(1,1) = x1*cos(Omega)-y1*cos(i)*sin(Omega);
satp(2,1) = x1*sin(Omega)+y1*cos(i)*cos(Omega);
satp(3,1) = y1*sin(i);

%- Compute satellite velocity
Ek          =   E;
fk          =   v;
phik        =   phi;
uk          =   u;
x1k         =   x1;
y1k         =   y1;
ik          =   i;
xk          =   satp(1);
yk          =   satp(2);
%
Omegak      =   Omega;
Ek_dot      =   n/(1-ecc*cos(Ek));
fk_dot      =   sin(Ek)*Ek_dot*(1+ecc*cos(fk)) / ((1-cos(Ek)*ecc)*sin(fk));
phik_dot    =   fk_dot;
uk_dot      =   phik_dot + 2*(cus*cos(2*phik)-cuc*sin(2*phik))*phik_dot;
rk_dot      =   A*ecc*sin(Ek)*Ek_dot + 2*(crs*cos(2*phik)-crc*sin(2*phik))*phik_dot;
ik_dot      =   idot + 2*(cis*cos(2*phik)-cic*sin(2*phik))*phik_dot;
Omegak_dot  =   Omegadot - Omegae_dot;
x1k_dot     =   rk_dot*cos(uk) - y1k*uk_dot;
y1k_dot     =   rk_dot*sin(uk) + x1k*uk_dot;
xk_dot      =   x1k_dot*cos(Omegak) - y1k_dot*cos(ik)*sin(Omegak) + y1k*sin(ik)*sin(Omegak)*ik_dot - yk*Omegak_dot;
yk_dot      =   x1k_dot*sin(Omegak) + y1k_dot*cos(ik)*cos(Omegak) - y1k*sin(ik)*ik_dot*cos(Omegak) + xk*Omegak_dot;
zk_dot      =   y1k_dot*sin(ik) + y1k*cos(ik)*ik_dot;

satv        =   [xk_dot; yk_dot; zk_dot];



%%%%%%%%% end satpos.m %%%%%%%%%
