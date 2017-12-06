function mask = clean_mask(x,y,mask,bound_ingrid,lim,offset)

% -------------------------------------------------------------------------
%|                                                                        |
%|                    +----------------------------+                      |
%|                    | GRIDGEN          NOAA/NCEP |                      |
%|                    |                            |                      |
%|                    | Last Update :  29-Mar-2013 |                      |
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
%| mask = clean_mask(x,y,mask,bound_ingrid,lim,offset)                    |
%|                                                                        |
%| INPUT                                                                  |
%|  x            : A 2D array specifying the longitudes of each cell      |
%|  y            : A 2D array specifying the longitudes of each cell      |
%|  mask         : A 2D mask array of size (Nx,Ny) of initial mask values |
%|                 with wet cells having a flag of 1 and dry cells a flag |
%|                 value of 0. The initial mask can be all set to 1 or for|
%|		   better efficiency be based on the bathymetry data      |
%|  bound_ingrid : Data structure array of boundary polygons that lie     |
%|                 inside the grid                                        |
%|  lim          : Fraction value between 0 and 1 that define the cut-off |
%|                 limit for the proportion of a cell that has to lie     | 
%|                 inside boundary polygons for the cell to be marked dry |
%|  offset       : An additional width around the boundaries to check if  |
%|                 wet cells are being crossed by the boundary. Typically |
%|                 set it to largest grid resolution                      |
%|                                                                        |
%|  OUTPUT                                                                |
%|    mask       : New land/sea mask array that is generated after        |
%|                 checking with the boundary polygons                    |
%|                                                                        |
% -------------------------------------------------------------------------

%@@@ Determine array limits

N1 = length(bound_ingrid);
[Ny,Nx] = size(x);

%@@@ Initialize 2D array specifying proportion of cell inside boundary(ies)

mask_obstr = zeros(Ny,Nx);
mask_status = zeros(Ny,Nx);
mask_points = cell(Ny,Nx);
mask_x = cell(Ny,Nx);
mask_y = cell(Ny,Nx);

itmp = 0;

%@@@ Loop through all the boundaries

for i = 1:N1

    %@@@ Determine limits of boundary

    west = bound_ingrid(i).west-offset;
    south = bound_ingrid(i).south-offset;
    east = bound_ingrid(i).east+offset;
    north = bound_ingrid(i).north+offset;

    %@@@ Determine the longitude and lattitude cells that lie within the 
    %@@@ boundary range
    
    py = [south south north north south]; 
    px = [west east east west west];
    
    clear in_bnd cell_loc row_pos column_pos;
    in_bnd = inpolygon(x,y,px,py);
    cell_loc = find(in_bnd > 0);
    [row_pos,column_pos] = ind2sub(size(in_bnd),cell_loc);
    
    %@@@ Loop through all the cells that lie within the boundary range
    
    N_cell = length(cell_loc);
    
    for cell_indx = 1:N_cell
        
        k = row_pos(cell_indx);
        j = column_pos(cell_indx);
        
        %@@@ Check if cell is within a boundary only if it is a wet cell
            
        if (mask(k,j) == 1)
            
            %@@@ distribute points inside a cell (only have to do this once)
            
            if (mask_status(k,j) == 0)
                
                [c1,c2,c3,c4,~,~] = compute_cellcorner(x,y,j,k,Nx,Ny);
                xmin = min([c1(1) c2(1) c3(1) c4(1)]);
                xmax = max([c1(1) c2(1) c3(1) c4(1)]);
                ymin = min([c1(2) c2(2) c3(2) c4(2)]);
                ymax = max([c1(2) c2(2) c3(2) c4(2)]);
                
                xtt = linspace(xmin,xmax,6);
                ytt = linspace(ymin,ymax,6);
                [xtt2,ytt2] = meshgrid(xtt,ytt);
                px1 = [c4(1) c1(1) c2(1) c3(1) c4(1)]';
                py1 = [c4(2) c1(2) c2(2) c3(2) c4(2)]';
                clear loc in_cell;
                in_cell = inpolygon(xtt2,ytt2,px1,py1);
                loc = find(in_cell > 0);
                xval = xtt2(loc);
                yval = ytt2(loc);
                
                mask_x(k,j) = {xval};
                mask_y(k,j) = {yval};
                mask_points(k,j) = {zeros(size(xval))};
                mask_status(k,j) = 1;
                
            end;
            
            %@@@ Check if cell points within boundary
            
            clear xt yt status loc inout;
                    
            xt = cell2mat(mask_x(k,j));
            yt = cell2mat(mask_y(k,j));
            status = cell2mat(mask_points(k,j));
            Na = numel(xt);
            loc = find(status == 0);
            if (~isempty(loc))
                inout = inpolygon(xt(loc),yt(loc),bound_ingrid(i).x,...
                    bound_ingrid(i).y);
                status(loc) = inout;
                clear loc;
            end;
            
            %@@@ Update mask_obstr. This variable mantains an estimate for 
            %@@@ proportion of cell that is covered by boundaries

            loc = find(status > 0);
            mask_obstr(k,j) = round(length(loc)/Na*10)/10;

            %@@@ If proportion exceeds user specified limit then mark the cell dry

            if (mask_obstr(k,j) >= lim) 
                mask(k,j) = 0;
            else
                mask_points(k,j) = {status};
            end;
            
        end; %@@@ Corresponds to wet cell check

    end;  %@@@ Corresponds to loop of cells within the boundary

    %@@@ Counter to update proportion of land sea mask clean up

    itmp_prev = itmp;
    itmp = floor(i/N1*100);
    if (mod(itmp,5) == 0 && itmp_prev ~= itmp)
        fprintf(1,'Completed %d per cent of land sea mask clean up\n',itmp);
    end;

end; %@@@ Corresponds to for loop of all the boundaries

return;
