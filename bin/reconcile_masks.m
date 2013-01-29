function [m1_out,m2_out] = reconcile_masks(m1,lon1,lat1,m2,lon2,lat2 )

% -------------------------------------------------------------------------
%|                                                                        |
%|                    +----------------------------+                      |
%|                    | GRIDGEN          NOAA/NCEP |                      |
%|                    |                            |                      |
%|                    | Last Update :  23-Oct-2012 |                      |
%|                    +----------------------------+                      |
%|                     Distributed with WAVEWATCH III                     |
%|                                                                        |
%|                 Copyright 2009 National Weather Service (NWS),         |
%|  National Oceanic and Atmospheric Administration.  All rights reserved.|
%|                                                                        |
%| DESCRIPTION                                                            |
%|  The aim of this function is to reconcile the masks of two grids with  |
%|  similar resolution in the areas of overlap. This is only needed for   |
%|  grids of similar resolution with overlapping domains ( where the      |
%|  the multi-grid version of WAVEWATCH tries to reconcile the solutions) |
%|                                                                        |
%|  [m1_out,m2_out] = reconcile_masks(m1,lon1,lat1,m2,lon2,lat2)          |
%|                                                                        |
%| INPUT                                                                  |
%|  m1           : Mask for grid1                                         |
%|  lon1         : Longitude (x) array for grid1                          |
%|  lat1         : Latittude (y) array for grid1                          |
%|  m2           : Mask for grid2                                         |
%|  lon2         : Longitude (x) array for grid2                          |
%|  lat2         : Latittude (y) array for grid2                          |
%|                                                                        |
%| OUTPUT                                                                 |
%|  m1_out       : Re-conciled mask for grid1                             |
%|  m2_out       : Re-conciled mask for grid2                             |
%|                                                                        |
% -------------------------------------------------------------------------

%@@@ Determine overlap points

[Nx1,Ny1] = size(m1);
[Nx2,Ny2] = size(m2);

if (Nx1*Ny1 <= Nx2*Ny2)
    px = [lon2(1) lon2(end) lon2(end) lon2(1) lon2(1)];
    py = [lat2(1) lat2(1) lat2(end) lat2(end) lat2(1)];
    [x,y] = meshgrid(lon1,lat1);
    lonb = lon2;
    latb = lat2;
    dx = lon1(2) - lon1(1);
    dy = lat1(2) - lat1(1);
    Nt = Nx1*Ny1;
    mt = m1;
    mb = m2;
else
    px = [lon1(1) lon1(end) lon1(end) lon1(1) lon1(1)];
    py = [lat1(1) lat1(1) lat1(end) lat1(end) lat1(1)];
    [x,y] = meshgrid(lon2,lat2);
    lonb = lon1;
    latb = lat1;
    dx = lon2(2) - lon2(1);
    dy = lat2(2) - lat2(1);
    Nt = Nx2*Ny2;
    mt = m2;
    mb = m1;
end

%@@@ Determine region of overlap

inout = inpolygon(x,y,px,py);
loc = find(inout > 0);
N = length(loc);

fprintf(1,' Found %d per cent of grid overlap points\n',round(N*100./Nt));

%@@@ Check masks for overlap points

for i = 1:N
    lon = x(loc(i));
    lat = y(loc(i));
    [xmin,jpos] = min(abs(lon-lonb));
    [ymin,ipos] = min(abs(lat-latb));
    if (max(abs(xmin),abs(ymin)) > min(abs(dx),abs(dy)))
        fprintf(1,' Error ! Closest points are too far ! \n');
        return;
    end;
    
    if (mt(loc(i)) ~= mb(ipos,jpos))
        mt(loc(i)) = 0;
        mb(ipos,jpos) = 0;
    end;
end;

if (Nx1*Ny1 <= Nx2*Ny2)
    m1_out = mt;
    m2_out = mb;
else
    m1_out = mb;
    m2_out = mt;
end;

return;

