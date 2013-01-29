function mask = clean_mask(lon,lat,mask,bound_ingrid,lim)

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
%| This function checks all the wet cells in a 2D mask array and          |  
%| determines if they lie outside the boundary polygons or not            |
%|                                                                        |
%| mask = clean_mask(lon,lat,mask,bound_ingrid,lim)                       |
%|                                                                        |
%| INPUT                                                                  |
%|  lon          : An array of length Nx consisting of the longitude (x)  |
%|                 coordinates of the grid                                |
%|  lat          : An array of length Ny consisting of the lattitude (y)  |
%|                 coordinates of the grid                                |
%|  mask         : A 2D mask array of size (Nx,Ny) of initial mask values |
%|                 with wet cells having a flag of 1 and dry cells a flag |
%|                 value of 0. The initial mask can be all set to 1 or for|
%|		           better efficiency be based on the bathymetry data      |
%|  bound_ingrid : Data structure array of boundary polygons that lie     |
%|                 inside the grid                                        |
%|  lim          : Fraction value between 0 and 1 that define the cut-off |
%|                 limit for the proportion of a cell that has to lie     | 
%|                 inside boundary polygons for the cell to be marked dry |
%|                                                                        |
%|  OUTPUT                                                                |
%|    mask       : New land/sea mask array that is generated after        |
%|                 checking with the boundary polygons                    |
%|                                                                        |
% -------------------------------------------------------------------------

%@@@ Determine array limits

N1 = length(bound_ingrid);
Ny = length(lat);
Nx = length(lon);

dx = lon(2)-lon(1);
dy = lat(2) - lat(1);

%@@@ Initialize 2D array specifying proportion of cell inside boundary(ies)

mask_obstr = zeros(Ny,Nx);
mask_status = zeros(Ny,Nx);
mask_points = cell(Ny,Nx);
mask_x = cell(Ny,Nx);
mask_y = cell(Ny,Nx);

coord_start = 0;
coord_end = 360;

itmp = 0;

%@@@ Loop through all the boundaries

for i = 1:N1

    %@@@ Determine limits of boundary

    west = bound_ingrid(i).west-dx;
    south = bound_ingrid(i).south-dy;
    east = bound_ingrid(i).east+dx;
    north = bound_ingrid(i).north+dy;
    
    if (west < coord_start)
        west = west+360;
    end;
    if (east > coord_end)
        east = east - 360;
    end;

    %@@@ Determine the longitude and lattitude cells that lie within the 
    %@@@ boundary range
    
    lat_loc = find(lat >= south & lat <= north);
    
    if (west < east)
        lon_loc = find(lon >= west & lon <= east);
    else        
        loc1 = find(lon >= west & lon <= coord_end);
        loc2 = find(lon >= coord_start & lon <= east);
        lon_loc = [loc1 loc2];
    end;
    
    Nx = length(lon_loc);
    Ny = length(lat_loc);

    %@@@ Loop through all the cells that lie within the boundary range

    cell_count=0;
   
    for j = 1:Ny 
        for k = 1:Nx
            xp = lon_loc(k);
            yp = lat_loc(j);
            cell_count = cell_count+1;
              
            %@@@ Check if cell is within a boundary only if it is a wet cell
            
            if (mask(yp,xp) == 1)
                
                %@@@ Determine the edges of the cell accounting for cell 
                %@@@ wrapping around
                if (mask_status(yp,xp) == 0)
                    
                    lon_start = lon(xp)-dx/2.0;
                    lat_start = lat(yp)-dy/2.0;
                    lon_end = lon(xp)+dx/2.0;
                    lat_end = lat(yp)+dy/2.0;
                    if (lon_start < coord_start)
                        lon_start = lon_start+360;
                    end;
                    if (lon_end > coord_end)
                        lon_end = lon_end-360;
                    end;
                
                    %@@@ Check to see what part of cell is inside domain to
                    %@@@ proceed               
                    
                    ytt = lat_start:dy/5:lat_end; 
                    if (lon_start < lon_end)
                        xtt = lon_start:dx/5:lon_end;
                    else
                        xtta = lon_start:dx/5:coord_end;
                        xttb = (coord_start+dx/5):dx/5:lon_end;
                        xtt = [xtta xttb];
                    end;
 
                    [xtt2,ytt2] = meshgrid(xtt,ytt);
                    xtt2 = reshape(xtt2,numel(xtt2),1);
                    ytt2 = reshape(ytt2,numel(ytt2),1);
                    mask_x(yp,xp) = {xtt2};
                    mask_y(yp,xp) = {ytt2};
                    mask_points(yp,xp) = {zeros(size(xtt2))};
                    mask_status(yp,xp) = 1;
                end;
                   
                clear xt yt status loc inout;
                    
                xt = cell2mat(mask_x(yp,xp));
                yt = cell2mat(mask_y(yp,xp));
                status = cell2mat(mask_points(yp,xp));
                Na = numel(xt);
                    
                loc = find(status == 0);
                if (~isempty(loc))
                    inout = inpolygon(xt(loc),yt(loc),bound_ingrid(i).x,...
                        bound_ingrid(i).y);
                    status(loc) = inout;
                    clear loc;
                end;
                    
                %@@@ Update mask_obstr. This variable mantains  
                %@@@ an estimate for proportion of cell that is 
                %@@@ covered by boundaries

                loc = find(status > 0);
                mask_obstr(yp,xp) = round(length(loc)/Na*10)/10;

                %@@@ If proportion exceeds user specified limit 
                %@@@ then mark the cell dry

                if (mask_obstr(yp,xp) >= lim) 
                    mask(yp,xp) = 0;
                else
                    mask_points(yp,xp) = {status};
                end;                              
                
            end; %@@@ Corresponds to wet cell check

        end;  %@@@ Corresponds to for loop of longitudes within the boundary

    end;     %@@@ Corresponds to for loop of lattitudes within the boundary

    %@@@ Counter to update proportion of land sea mask clean up

    itmp_prev = itmp;
    itmp = floor(i/N1*100);
    if (mod(itmp,5) == 0 && itmp_prev ~= itmp)
        fprintf(1,'Completed %d per cent of land sea mask clean up\n',itmp);
    end;

end; %@@@ Corresponds to for loop of all the boundaries

return;
