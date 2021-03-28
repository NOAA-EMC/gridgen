function [LAT,LON]=cart2ll_polarstereo(x,y,earth_radius,eccentricity,north)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Ali Abdolali EMC/NCEP/NOAA ali.abdolali@noaa.gov 26, March 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%This script transforms map coordinates to lat/lon data for a polar stereographic system
%Equations are taken form: Map Projections - A Working manual - by J.P. Snyder. 1987 

%earth_radius=6378137.0; %radius of ellipsoid, WGS84
%eccentricity=0.08181919; %eccentricity, WGS84
%north=1;%if south, north=0

if north==1
LAT_C=70; %latitude of true scale in the north hemisphere
else
LAT_C=-70; %latitude of true scale in the south hemisphere
end

LON_C=0; %meridian along positive Y axis

%convert to radians
LAT_C=deg2rad(LAT_C);
LON_C=deg2rad(LON_C);

%plus in the north or minus in the south

    pm=sign(LAT_C); 
    LAT_C=pm*LAT_C;
    LON_C=pm*LON_C;
    x=pm*x;
    y=pm*y;

THETA=(sqrt(x.^2+y.^2))*(tan(pi/4-LAT_C/2)/((1-eccentricity*sin(LAT_C))...
    /(1+eccentricity*sin(LAT_C)))^(eccentricity/2))/(earth_radius*...
    (cos(LAT_C)/sqrt(1-eccentricity^2*(sin(LAT_C))^2)));

% find LAT with a series intiated from, a first guess for phi =pi/2 - 2 * 
%atan(t) with a threshold of pi*1e-8

LAT=(pi/2 - 2 * atan(THETA))+(eccentricity^2/2 + 5*eccentricity^4/24 + ...
    eccentricity^6/12 + 13*eccentricity^8/360)*sin(2*(pi/2 - 2 * atan(THETA)))...
    + (7*eccentricity^4/48 + 29*eccentricity^6/240 + ...
    811*eccentricity^8/11520)*sin(4*(pi/2 - 2 * atan(THETA)))...
    + (7*eccentricity^6/120+81*eccentricity^8/1120)*sin(6*(pi/2 - 2 * atan(THETA)))...
    + (4279*eccentricity^8/161280)*sin(8*(pi/2 - 2 * atan(THETA)));

LON=LON_C + atan2(x,-y);
% correct the signs and phase
LAT=pm*LAT;
% LAT_alt=pm*LAT_alt;
LON=pm*LON;
LON=mod(LON+pi,2*pi)-pi; %want longitude in the range -pi to pi
%convert back to degrees
LAT=rad2deg(LAT);
LON=rad2deg(LON);

LON(LON<0)=LON(LON<0)+360; %make sure longitude is from 0-360 deg

