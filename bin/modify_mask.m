function m_new = modify_mask(m,lon,lat,mb,lonb,latb,igl,px,py)

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
%| This routine was developed for WAVEWATCH III v3.10 or higher           | 
%| (multi-grid version) where the traditional mask with values of 0 and 1 | 
%| for land and water are modified to give values of 0,1,2 and 3 for cells| 
%| that are on land, water, boundary or to be ignored, respectively. These|
%| masks are needed for grids that are nested with a larger grid from     |
%| which they get their boundary information                              |
%|                                                                        |
%| m_new = modify_mask(m,lon,lat,mb,lonb,latb,igl,[px,py])                |
%|                                                                        |
%| INPUT                                                                  |
%|   m     : 2D land sea mask for grid                                    |
%|   lon   : longitude (x) coordinates of grid                            |
%|   lat   : lattitude (y) coordinates of grid                            |
%|   mb    : 2D mask for base grid (grid with which boundary data is      |
%|           exchanged)                                                   |
%|   lonb  : longitude (x) coordinates of base grid                       |
%|   latb  : lattitude (y) coordinates of base grid                       |
%|   igl   : flag indicating if base grid is global (1) or not (0)        |
%|   px,py : Optional x,y coordinates of polygon defining the region of   |
%|           active computation in the grid. If ommitted the entire grid  |
%|           is used.                                                     |
%|                                                                        |
%| OUTPUT                                                                 |
%|   m_new : New 2D mask file with values ranging from 0-3.               |
%|              0 -> Active (within the computational region) dry cells   |
%|              1 -> Active wet cells                                     |
%|              2 -> Boundary cells (for data exchange)                   |
%|              3 -> Inactive cells                                       |
%|                                                                        |
% ------------------------------------------------------------------------

%@@@ Determine the number of inputs

narg = nargin;
edge_boundary = 1;
m_new = m;

%if (narg == 9)   
 %   px = varargin(1);
 %   py = varargin(2);
    edge_boundary = 0;
    in_points = inpolygon(lon,lat,px,py);
    m_new(in_points == 0) = 3;
%end;

fprintf(1,'Finished initial allocation ....\n');  

[Ny,Nx] = size(lon);

% moving along the rows

for j = 1:Ny
    if ( m_new(j,1) == 1)
        m_new(j,1) = 2;
    end;
    if ( m_new(j,Nx) == 1)
        m_new(j,Nx) = 2;
    end;
end;

% moving along the columns

for k = 1:Nx
    if ( m_new(1,k) == 1)
        m_new(1,k) = 2;
    end;    
    if ( m_new(Ny,k) == 1)
        m_new(Ny,k) = 2;
    end;    
end;
    
% Identifying boundaries inside the domain

[rows,cols] = find(m_new ==3);
Nr = length(rows);
for i = 1:Nr
    ny_down = max(1,rows(i)-1);
    ny_up = min(Ny,rows(i)+1);
    nx_left = max(1,cols(i)-1);
    nx_right = min(Nx,cols(i)+1);
    found_wet = 0;
    j=rows(i);
    k=cols(i);
    if (m_new(j,nx_left) == 1 || m_new(j,nx_right) == 1 || m_new(ny_down,k) == 1 || m_new(ny_up,k) == 1)
        found_wet = 1;
    end;
   % for j = ny_down:ny_up
   %     for k = nx_left:nx_right
   %         if (m_new(j,k) == 1)
   %             found_wet = 1;
   %         end;
   %     end;
   %end;
    if (found_wet == 1)
        m_new(rows(i),cols(i)) = 2;
    end;
end;

clear rows cols;

% Locating the correct values in the base    

[rows,cols] = find(m_new == 2);

Nr = length(rows);
[Nyb,Nxb] = size(lonb);

dxb = lonb(1,2)-lonb(1,1);
dyb = latb(2,1)-latb(1,1);

for i = 1:Nr
    
    y = lat(rows(i),cols(i));
    x = lon(rows(i),cols(i));
    
    ry = (y-latb(1,1))/dyb;
    jy = 1 + floor(ry);
    ry = ry - (jy-1);
    if (ry < 0)
        ry = 1 - ry;
        jy = jy-1;
    end;
    if (jy == 0 && abs(ry-1) < 0.05) 
        jy = 1;
        ry = 0;
    end;
    if (jy == Nyb && abs(ry) < 0.05)
        jy = jy-1;
        ry = 1;
    end;
    if (jy < 1 || jy >= Nyb)
        m_new(rows(i),cols(i)) = 3;
        continue;
    end;

    rx = (x-lonb(1,1))/dxb;
    jx = 1 + floor(rx);
    rx = rx - (jx-1);
 
    if (igl ~= 1) 
        if (jx == 0 && abs(rx-1) < 0.05)
            jx = 1;
            rx = 0.;
        end;
        if (jx == Nxb && abs(rx) < 0.05)
            jx = jx-1;
            rx = 1;
        end;
        if (jx < 1 || jx >= Nxb)
            m_new(rows(i),cols(i)) = 3;
            continue;
        end;
    else
        jx = 1 + mod(jx-1,Nxb);
    end;

    jx1 = jx;
    if (igl == 1) 
        jx2 = 1 + mod(jx,Nxb);
    else
        jx2 = jx + 1;
    end;

    jy1 = jy;
    jy2 = jy+1;

    flagok = (abs(mb(jy1,jx1)) == 1 || abs(mb(jy1,jx1)) == 2 || ...
        (1-rx)*(1-ry) < 0.05 ) & (abs(mb(jy1,jx2)) == 1 || ...
        abs(mb(jy1,jx2)) == 2 || rx*(1-ry) < 0.05 ) & ...
        (abs(mb(jy2,jx1)) == 1 || abs(mb(jy2,jx1)) == 2 || ...
        (1-rx)*ry < 0.05 ) & (abs(mb(jy2,jx2)) == 1 || ...
        abs(mb(jy2,jx2)) == 2 || rx*ry < 0.05 ); 
               
    if (flagok == 0)
        m_new(rows(i),cols(i)) = 3;
        continue;
    end;
end;


return;

