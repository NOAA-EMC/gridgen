function [sx,sy] = create_obstr(x,y,bound,mask,offset_left,offset_right)

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
%| This routine generates the 2D obstruction grid in x and y given a 2D   |
%| mask and set of boundary polygons. Obstructions are only generated for |
%| wet cells and obstructions for cells on either side of a dry cell are  |
%| also set to 0 (to prevent spurious suppression of swell near the coast)|  
%| The routine allows for the possibility of curvilinear coordinates and  |
%| locally rotates the grid to align the coordinates in lat/lon space with|
%| local p/q space which is determined by the 2D x and y matrices         |
%|                                                                        |
%| [sx,sy] = create_obstr_curv(x,y,bound,mask,coords,offset_left,...      |
%|                                                   offset_right)        |
%|                                                                        |
%| INPUT                                                                  |
%|  x            : A 2D array specifying the longitudes of each cell      | 
%|  y            : A 2D array specifying the lattitudes of each cell      | 
%|  bound        : Data structure array of boundary polygons              |
%|  mask         : 2D array of size (Nx,Ny) that determines land/sea mask |
%|  offset_left  : Flag to determine if neighbor to the left/down in x/y  |
%|                 should be considered. (0/1 = no/yes)                   |
%|  offset_right : Similar for neighbor to the right/up in x/y            |
%|                                                                        |
%| OUTPUT                                                                 |
%|  sx,sy        : 2D obstruction grids of size (Nx,Ny) for obstructions  |
%|                 in x and y. Values range from 0 for no obstruction to 1| 
%|                 for full obstruction                                   | 
%|                                                                        |
% -------------------------------------------------------------------------

%@@@ Initialize variables

[Ny,Nx] = size(x);
sx = zeros(Ny,Nx);
sy = sx;
loc = find(mask == 0);
sx(loc) = 0;
sy(loc) = 0;
clear loc;

cell = repmat(struct('px',[],'py',[],'angle',[],'width',[],'height',[],...
    'nx',[],'ny',[],'south_lim',[],'north_lim',[],'east_lim',[],...
    'west_lim',[],'bndx',[],'bndy',[]),Ny,Nx);

cell_bnd = mask;

loc_wet = find(mask ~= 0);
N_wet = length(loc_wet);
Nb = Nx*Ny;

fprintf(1,' Total Number of cells = %d\n',Nb);
fprintf(1,'   Number of wet cells = %d\n',N_wet);

%@@@ Set up the cells

for j = 1:Nx
    for k = 1:Ny
            
        [c1,c2,c3,c4,wdth,hgt] = compute_cellcorner(x,y,j,k,Nx,Ny);
        
        cell(k,j).px = [c4(1) c1(1) c2(1) c3(1) c4(1)]';
        cell(k,j).py = [c4(2) c1(2) c2(2) c3(2) c4(2)]';
        cell(k,j).angle = atan2(c1(2)-c4(2),c1(1)-c4(1));
       
        %cell(k,j).width = sqrt((c1(1)-c4(1))^2+(c1(2)-c4(2))^2);
        %cell(k,j).height = sqrt((c2(1)-c1(1))^2+(c2(2)-c1(2))^2);
        cell(k,j).width = wdth;
        cell(k,j).height = hgt;
        cell(k,j).nx = 0;
        cell(k,j).ny = 0;
        cell(k,j).south_lim = [];
        cell(k,j).north_lim = [];
        cell(k,j).east_lim = [];
        cell(k,j).west_lim = [];
        cell(k,j).bndx = [];
        cell(k,j).bndy = [];
        
    end;
end;

N = length(bound);

%@@@ Preparing the boundaries

fprintf(1,'Preparing the boundaries \n');
itmp = 0;

bnd_x = [];
bnd_y = [];
bnd_indx = [];

for i = 1:N
   bnd_x = [bnd_x;bound(i).x];
   bnd_y = [bnd_y;bound(i).y];
   bnd_indx = [bnd_indx;i*ones(bound(i).n,1)];
   itmp_prev = itmp;
    itmp = floor(i/N*100);
    if (mod(itmp,10) == 0 && itmp_prev ~= itmp)
        fprintf(1,' Completed %d per cent of boundaries \n',itmp);
    end;
end;

%@@@ Loop through the wet cells and determine the boundaries that are
%@@@ within

fprintf(1,'Loop through the wet cells to identify boundaries \n');
itmp = 0;

for indx_wet = 1:N_wet
    
    cell_pos = loc_wet(indx_wet);
    
    %@@@ Rotation Matrix and origin (for each cell the first corner is taken as
    %@@@ the origin about which the rotation takes place)
    
    RM = [cos(cell(cell_pos).angle) -sin(cell(cell_pos).angle) ; ...
        sin(cell(cell_pos).angle) cos(cell(cell_pos).angle)];
    x0 = cell(cell_pos).px(1);
    y0 = cell(cell_pos).py(1);
    px = cell(cell_pos).px;
    py = cell(cell_pos).py;  
    
    %@@@ Find the boundaries inside the cell
    
    clear Nbnds bnds in_box;
    in_box = inpolygon(bnd_x,bnd_y,px,py);
    bnds = unique(bnd_indx(in_box > 0));
    Nbnds = length(bnds);
    cell_bnd(cell_pos) = Nbnds;
    
    for i = 1:Nbnds
        
        indx_bnd = bnds(i);
        clear in_box2 in_box_coords;
        in_box2 = inpolygon(bound(indx_bnd).x,bound(indx_bnd).y,px,py);
        in_box_coords = find(in_box2 > 0);
        
        if (~isempty(in_box_coords))
            xt = bound(indx_bnd).x(in_box_coords);
            yt = bound(indx_bnd).y(in_box_coords);
                   
            % rotate points within the cell (including the corners)
            
            clear tmp
            tmp = [(xt-x0) (yt-y0)]*RM;
            xt = tmp(:,1);
            yt = tmp(:,2);
            clear tmp;
            tmp = [(px-x0) (py-y0)]*RM;
            pxt = tmp(:,1);
            pyt = tmp(:,2);
            
            % Boundary segments
            
            south_limit = min(yt)/cell(cell_pos).height;
            north_limit = max(yt)/cell(cell_pos).height;
            west_limit = min(xt)/cell(cell_pos).width;
            east_limit = max(xt)/cell(cell_pos).width;
            
            south_limit = max(0.0,south_limit);
            south_limit = min(1.0,south_limit);
            north_limit = max(0.0,north_limit);
            north_limit = min(1.0,north_limit);
            west_limit = max(0.0,west_limit);
            west_limit = min(1.0,west_limit);
            east_limit = max(0.0,east_limit);
            east_limit = min(1.0,east_limit);
            
            cell(cell_pos).nx = cell(cell_pos).nx+1;
            idx = cell(cell_pos).nx;
            
            cell(cell_pos).south_lim(idx) = south_limit;
            cell(cell_pos).north_lim(idx) = north_limit;
            cell(cell_pos).bndx(idx) = indx_bnd;
            cell(cell_pos).ny = cell(cell_pos).ny+1;
            idx = cell(cell_pos).ny;
            cell(cell_pos).east_lim(idx) = east_limit;
            cell(cell_pos).west_lim(idx) = west_limit;
            cell(cell_pos).bndy(idx) = indx_bnd;
            
        end;
    end;
    
    itmp_prev = itmp;
    itmp = floor(indx_wet/N_wet*100);
    if (mod(itmp,5) == 0 && itmp_prev ~= itmp)
        fprintf(1,' Completed %d per cent\n',itmp);
    end;
   
end;

%@@@ Loop through all the wet cells with boundaries and move boundary segments 
%@@@ that are part of the same boundary and cross neighboring cells. This is 
%@@@ done to prevent double counting in building the obstruction grids. Only the
%@@@  immediate neighbors are checked when moving the segments

loc_bnd = find(cell_bnd ~= 0);
N_bnd = length(loc_bnd);
[row_bnd,column_bnd] = ind2sub(size(cell_bnd),loc_bnd);
fprintf(1,'Number of wet cells enclosing boundaries = %d\n',N_bnd);

fprintf(1,'Computing subgrid obstruction masks for these cells \n');

for indx_bnd = 1:N_bnd
    
    j = column_bnd(indx_bnd);
    k = row_bnd(indx_bnd);
    
    %@@@ First check the neighbors in x direction
            
    if (j < Nx)
        jj = j+1;

        %@@@ Check to see that boundary segments are present in both cells

        if (cell(k,j).nx~=0 && cell(k,jj).nx~=0)

            %@@@ Save information to temporary variables

            set1 = cell(k,j);
            set2 = cell(k,jj);
            found_common = 0;

            %@@@ Loop through boundary segments and move segments  
            %@@@ of common boundaries to the cell with the larger 
            %@@@ segment

            for l = 1:set1.nx
                for m = 1:set2.nx
                    if (set1.bndx(l) == set2.bndx(m))
                        if ((set1.north_lim(l)-set1.south_lim(l)) >= ...
                                (set2.north_lim(m)-set2.south_lim(m)))
                            set1.north_lim(l) = max([set1.north_lim(l) ...
                                set2.north_lim(m)]);
                            set1.south_lim(l) = min([set1.south_lim(l) ...
                                set2.south_lim(m)]);
                            for n = (m+1):set2.nx
                                set2.bndx(n-1) = set2.bndx(n);
                                set2.north_lim(n-1) = set2.north_lim(n);
                                set2.south_lim(n-1) = set2.south_lim(n);
                            end;
                            set2.nx = set2.nx-1;
                        else
                            set2.north_lim(m) = max([set1.north_lim(l) ...
                                set2.north_lim(m)]);
                            set2.south_lim(m) = min([set1.south_lim(l) ...
                                set2.south_lim(m)]);
                            for n = (l+1):set1.nx
                                set1.bndx(n-1) = set1.bndx(n);
                                set1.north_lim(n-1) = set1.north_lim(n);
                                set1.south_lim(n-1) = set1.south_lim(n);
                            end;
                            set1.nx = set1.nx-1;
                        end;
                        found_common = 1;
                        break;
                    end;
                end;
            end;

            %@@@ Write cell information back from temporary variables 
            %@@@ if common boundaries were found

            if (found_common == 1)
                cell(k,j).bndx = [];
                cell(k,j).north_lim = [];
                cell(k,j).south_lim = [];
                cell(k,jj).bndx = [];
                cell(k,jj).north_lim = [];
                cell(k,jj).south_lim = [];
                        
                cell(k,j).nx = set1.nx;
                cell(k,j).bndx = set1.bndx(1:set1.nx);
                cell(k,j).north_lim = set1.north_lim(1:set1.nx);
                cell(k,j).south_lim = set1.south_lim(1:set1.nx);
                        
                cell(k,jj).nx = set2.nx;
                cell(k,jj).bndx = set2.bndx(1:set2.nx);
                cell(k,jj).north_lim = set2.north_lim(1:set2.nx);
                cell(k,jj).south_lim = set2.south_lim(1:set2.nx);
                
            end;

        end;

    end;
            
    %@@@ Now check the neighbors in y direction (comments are omitted
    %@@@ because the operations are same as for x)
            
    if (k < Ny)
        kk = k+1;

        if (cell(k,j).ny~=0 && cell(kk,j).ny~=0)
            set1 = cell(k,j);
            set2 = cell(kk,j);
            found_common = 0;

            for l = 1:set1.ny
                for m = 1:set2.ny

                    if (set1.bndy(l) == set2.bndy(m))
                        if ((set1.east_lim(l)-set1.west_lim(l)) >= ...
                                (set2.east_lim(m)-set2.west_lim(m)))
                            set1.east_lim(l) = max([set1.east_lim(l) ...
                                set2.east_lim(m)]);
                            set1.west_lim(l) = min([set1.west_lim(l) ...
                                set2.west_lim(m)]);
                            for n = (m+1):set2.ny
                                set2.bndy(n-1) = set2.bndy(n);
                                set2.east_lim(n-1) = set2.east_lim(n);
                                set2.west_lim(n-1) = set2.west_lim(n);
                            end;
                            set2.ny = set2.ny-1;
                        else
                            set2.east_lim(m) = max([set1.east_lim(l) ...
                                set2.east_lim(m)]);
                            set2.west_lim(m) = min([set1.west_lim(l) ...
                                set2.west_lim(m)]);
                            for n = (l+1):set1.ny
                                set1.bndy(n-1) = set1.bndy(n);
                                set1.east_lim(n-1) = set1.east_lim(n);
                                set1.west_lim(n-1) = set1.west_lim(n);
                            end;
                            set1.ny = set1.ny-1;
                        end;
                        found_common = 1;
                        break;
                    end;
                end;
            end;

            if (found_common == 1)
                cell(k,j).bndy = [];
                cell(k,j).east_lim = [];
                cell(k,j).west_lim = [];
                cell(kk,j).bndy = [];
                cell(kk,j).east_lim = [];
                cell(kk,j).west_lim = [];
            
                cell(k,j).ny = set1.ny;
                cell(k,j).bndy = set1.bndy(1:set1.ny);
                cell(k,j).east_lim = set1.east_lim(1:set1.ny);
                cell(k,j).west_lim = set1.west_lim(1:set1.ny);
                        
                cell(kk,j).ny = set2.ny;
                cell(kk,j).bndy = set2.bndy(1:set2.ny);
                cell(kk,j).east_lim = set2.east_lim(1:set2.ny);
                cell(kk,j).west_lim = set2.west_lim(1:set2.ny);                       
            end;

        end;

    end; %@@@ End of if condition to check neighbor in y direction        
    
end;  %@@@ For loop along all the wet cells with boundaries


%@@@ Loop through the cells a second time to reduce overlapping segments 
%@@@ into individual non-overlapping segments

for indx_bnd = 1:N_bnd
    
    j = column_bnd(indx_bnd);
    k = row_bnd(indx_bnd);

    %@@@ First for x obstruction

    if (cell(k,j).nx ~= 0)

    %@@@ Reduction of overlapping segments only done if there are more than
    %@@@ 1 segment

        if (cell(k,j).nx > 1)
            n_segs = cell(k,j).nx;
            baseseg_n = cell(k,j).north_lim;
            baseseg_s = cell(k,j).south_lim;
            cell(k,j).north_lim = [];
            cell(k,j).south_lim = [];
            ind_segs=0;
            indseg_n = [];
            indseg_s = [];

            %@@@ Loop till all the segments in the cell done

            while (n_segs > 0)
                overlap_found = 0;                       
                if (n_segs > 1)
                    for l = 2:n_segs
                        if (baseseg_n(1) >= baseseg_s(l) && ...
                                baseseg_s(1) <= baseseg_n(l))
                            baseseg_n(1) = max([baseseg_n(1) baseseg_n(l)]);
                            baseseg_s(1) = min([baseseg_s(1) baseseg_s(l)]);
                            overlap_found = 1;
                            if (l == n_segs)
                                n_segs = n_segs-1;
                            else
                                for m = l+1:n_segs
                                    baseseg_n(m-1) = baseseg_n(m);
                                    baseseg_s(m-1) = baseseg_s(m);
                                end;
                                n_segs = n_segs-1;
                            end;
                            break;
                        end;
                    end;
                end;
                        
                if (n_segs == 1)
                    ind_segs = ind_segs+1;
                    indseg_n(ind_segs) = baseseg_n(1);
                    indseg_s(ind_segs) = baseseg_s(1);
                    n_segs = n_segs-1;
                else
                    if (overlap_found == 0)
                        ind_segs = ind_segs+1;
                        indseg_n(ind_segs) = baseseg_n(1);
                        indseg_s(ind_segs) = baseseg_s(1);
                        for l = 2:n_segs                                    
                            baseseg_n(l-1) = baseseg_n(l);
                            baseseg_s(l-1) = baseseg_s(l);
                        end;
                        n_segs = n_segs-1;
                    end;
                end;
            end;

            %@@@ Store the unique segments back in the cell variables

            cell(k,j).nx = ind_segs;
            cell(k,j).north_lim = indseg_n;
            cell(k,j).south_lim = indseg_s;

        end; %@@@ Corresponds to check for multiple segments in x

   end; %@@@ Corresponds to check for obstruction in x

   %@@@ Now check for overlapping segments in y 

   if (cell(k,j).ny ~= 0)

        if (cell(k,j).ny > 1)

            n_segs = cell(k,j).ny;
            baseseg_n = cell(k,j).east_lim;
            baseseg_s = cell(k,j).west_lim;
            cell(k,j).east_lim = [];
            cell(k,j).west_lim = [];
            ind_segs=0;
            indseg_n = [];
            indseg_s = [];

            while (n_segs > 0)
                overlap_found = 0;
                if (n_segs > 1)
                    for l = 2:n_segs
                        if (baseseg_n(1) >= baseseg_s(l) && ... 
                                baseseg_s(1) <= baseseg_n(l))
                            baseseg_n(1) = max([baseseg_n(1) baseseg_n(l)]);
                            baseseg_s(1) = min([baseseg_s(1) baseseg_s(l)]);
                            overlap_found = 1;
                            if (l == n_segs)
                                n_segs = n_segs-1;
                            else
                                for m = l+1:n_segs
                                    baseseg_n(m-1) = baseseg_n(m);
                                    baseseg_s(m-1) = baseseg_s(m);
                                end;
                                n_segs = n_segs-1;
                            end;
                            break;
                        end;
                    end;
                end;
                        
                if (n_segs == 1)
                    ind_segs = ind_segs+1;
                    indseg_n(ind_segs) = baseseg_n(1);
                    indseg_s(ind_segs) = baseseg_s(1);
                    n_segs = n_segs-1;
                else
                    if (overlap_found == 0)
                        ind_segs = ind_segs+1;
                        indseg_n(ind_segs) = baseseg_n(1);
                        indseg_s(ind_segs) = baseseg_s(1);
                        for l = 2:n_segs                                    
                            baseseg_n(l-1) = baseseg_n(l);
                            baseseg_s(l-1) = baseseg_s(l);
                        end;
                        n_segs = n_segs-1;
                    end;
                end;
            end;

            cell(k,j).ny = ind_segs;
            cell(k,j).east_lim = indseg_n;
            cell(k,j).west_lim = indseg_s;

        end;

   end; %@@@ End of overlapping segment computation in y        
    
end; %@@@ End of for loop along wet cells with boundaries


%@@@ Final loop through the cell rows and columns to construct the 
%@@@ obstruction grid accounting for neighboring cell information


for indx_bnd = 1:N_bnd
    
    j = column_bnd(indx_bnd);
    k = row_bnd(indx_bnd);

    %@@@ Computing x obstruction

    if (cell(k,j).nx ~=0)

        n_segs = cell(k,j).nx;
        baseseg_n = cell(k,j).north_lim;
        baseseg_s = cell(k,j).south_lim;
    
        no_boundary = 0;    

        %@@@ Compare boundary segments in cell with that of cell 
        %@@@ to left (if applicable)

        for off = 1:offset_left

            jj = j-off;
            if (jj >= 1)
                if (cell(k,jj).nx ~= 0)
                    set1 = cell(k,jj); 

                    %@@@ remove segments in shadow of previous cell  
                            
                    shadow_flags = zeros(size(baseseg_n));

                    for m = 1:n_segs
                        for l = 1:set1.nx
                            if (set1.north_lim(l) >= baseseg_n(m) ...
                                    && set1.south_lim(l) <= baseseg_s(m))
                                shadow_flags(m) = 1;
                                break;
                            end;
                        end;                                
                    end;

                    loc = find(shadow_flags == 0);
                    if (isempty(loc))
                        no_boundary = 1;
                        n_segs = 0;
                        baseseg_n = [];
                        baseseg_s = [];
                    elseif (length(loc) < n_segs)
                        tmp_n = baseseg_n;
                        tmp_s = baseseg_s;
                        baseseg_n = [];
                        baseseg_s = [];
                        baseseg_n = tmp_n(loc);
                        baseseg_s = tmp_s(loc);
                        n_segs = length(baseseg_n);
                        clear tmp_n tmp_s;
                    end;
                    clear loc shadow_flags;
                            
                    %@@@ Remove segments in previous cell that are
                    %@@@ shadows of present cell, 

                    if (~no_boundary)

                        shadow_flags = zeros(size(set1.north_lim));

                        for m = 1:set1.nx
                            for l = 1:n_segs
                                if (set1.north_lim(m) <= baseseg_n(l) ...
                                        && set1.south_lim(m) >= baseseg_s(l))
                                    shadow_flags(m) = 1;
                                    break;
                                end;
                            end;                                
                        end;

                        loc = find(shadow_flags == 0);
                        if (~isempty(loc))
                            if (length(loc) < set1.nx)
                                tmp_n = set1.north_lim;
                                tmp_s = set1.south_lim;
                                set1.north_lim = [];
                                set1.south_lim = [];
                                set1.north_lim = tmp_n(loc);
                                set1.south_lim = tmp_s(loc);
                                set1.nx = length(set1.north_lim);
                                clear tmp_n tmp_s;
                            end;
                            n_segs = n_segs+set1.nx;
                            baseseg_n = [baseseg_n set1.north_lim];
                            baseseg_s = [baseseg_s set1.south_lim];
                        end;
                        clear loc shadow_flags;
                    end;                           
                end;
            end;
        end; %@@@ End of search of segments in the left cell

        %@@@ Repeat for cell to the right (if applicable)

        if (~no_boundary)

            for off = 1:offset_right
                jj = j+off;
                if (jj <= Nx)
                    if (cell(k,jj).nx ~= 0)
                        set1 = cell(k,jj); 

                        %@@@ See if segments lie in shadow zone to 
                        %@@@ segments in right cell

                        shadow_flags = zeros(size(baseseg_n));
                        for m = 1:n_segs
                            for l = 1:set1.nx
                                if (set1.north_lim(l) >= baseseg_n(m) ...
                                        && set1.south_lim(l) <= baseseg_s(m))
                                    shadow_flags(m) = 1;
                                    break;
                                end;
                            end;                                
                        end;

                        %@@@ Set no_boundary flag if all segments removed

                        loc = find(shadow_flags == 0);
                        if (isempty(loc))
                            no_boundary = 1;
                            n_segs = 0;
                            baseseg_n = [];
                            baseseg_s = [];
                        elseif (length(loc) < n_segs)
                            tmp_n = baseseg_n;
                            tmp_s = baseseg_s;
                            baseseg_n = [];
                            baseseg_s = [];
                            baseseg_n = tmp_n(loc);
                            baseseg_s = tmp_s(loc);
                            n_segs = length(baseseg_n);
                            clear tmp_n tmp_s;
                        end;

                        clear loc shadow_flags;

                        %@@@ Remove segments from right cell that 
                        %@@@ are in the shadow zone 
                            
                        if (~no_boundary)
                            shadow_flags = zeros(size(set1.north_lim));
                            for m = 1:set1.nx
                                for l = 1:n_segs
                                    if (set1.north_lim(m) <= baseseg_n(l) ...
                                            && set1.south_lim(m) >= baseseg_s(l))
                                        shadow_flags(m) = 1;
                                        break;
                                    end;
                                end;                                
                            end;

                            loc = find(shadow_flags == 0);
                            if (~isempty(loc))
                                if (length(loc) < set1.nx)
                                    tmp_n = set1.north_lim;
                                    tmp_s = set1.south_lim;
                                    set1.north_lim = [];
                                    set1.south_lim = [];
                                    set1.north_lim = tmp_n(loc);
                                    set1.south_lim = tmp_s(loc);
                                    set1.nx = length(set1.north_lim);
                                    clear tmp_n;
                                    clear tmp_s;
                                end;
                                n_segs = n_segs+set1.nx;
                                baseseg_n = [baseseg_n set1.north_lim];
                                baseseg_s = [baseseg_s set1.south_lim];
                            end;
                            clear loc;
                            clear shadow_flags;
                        end;                           
                    end;
                end;
            end;
        end; %@@@ End of check of cells to the right
                
        %@@@ Now ready to build obstruction grid from the total 
        %@@@ set of segments

        if (~no_boundary)

            %@@@ Obstruction grid straightforward if n_segs == 1

            if (n_segs==1)
                
                sx(k,j) = baseseg_n(1)-baseseg_s(1);

            else    

                %@@@ Remove overlapping segments

                ind_segs=0;
                indseg_n = [];
                indseg_s = [];
                
                while (n_segs > 0)
                    overlap_found = 0;
                    if (n_segs > 1)
                        for l = 2:n_segs
                            if (baseseg_n(1) >= baseseg_s(l) && ...
                                    baseseg_s(1) <= baseseg_n(l))
                                baseseg_n(1) = max([baseseg_n(1) baseseg_n(l)]);
                                baseseg_s(1) = min([baseseg_s(1) baseseg_s(l)]);
                                overlap_found = 1;
                                if (l == n_segs)
                                    n_segs = n_segs-1;
                                else
                                    for m = l+1:n_segs
                                        baseseg_n(m-1) = baseseg_n(m);
                                        baseseg_s(m-1) = baseseg_s(m);
                                    end;
                                    n_segs = n_segs-1;
                                end;
                                break;
                            end;
                        end;
                    end;
                        
                    if (n_segs == 1)
                        ind_segs = ind_segs+1;
                        indseg_n(ind_segs) = baseseg_n(1);
                        indseg_s(ind_segs) = baseseg_s(1);
                        n_segs = n_segs-1;
                    else
                        if (overlap_found == 0)
                            ind_segs = ind_segs+1;
                            indseg_n(ind_segs) = baseseg_n(1);
                            indseg_s(ind_segs) = baseseg_s(1);
                            for l = 2:n_segs                                    
                                baseseg_n(l-1) = baseseg_n(l);
                                baseseg_s(l-1) = baseseg_s(l);
                            end;
                            n_segs = n_segs-1;
                        end;
                    end;
                end;
      
                %@@@ Compute the obstruction values from the 
                %@@@ independant segments

                for l = 1:ind_segs
                    sx(k,j) = sx(k,j) + (indseg_n(l)-indseg_s(l));
                end;

                cell(k,j).indsegs_north = indseg_n;
                cell(k,j).indsegs_south = indseg_s;
                clear baseseg_n baseseg_s baseseg_bound;

            end; %@@@ End of if statement to check number of segments

        end; %@@@ End of if statement to check boundaries in cell
 
    end; %@@@ End of if statement to check obstructions along x
            
    %@@@ Now computing obstruction in y in a similar fashion to x

    if (cell(k,j).ny~=0)

        n_segs = cell(k,j).ny;
        baseseg_n = cell(k,j).east_lim;
        baseseg_s = cell(k,j).west_lim;
                
        no_boundary = 0;    

        %@@@ First check with boundaries in cell below (if applicable)

        for off = 1:offset_left

            kk = k-off;

            if (kk >= 1)
                if (cell(kk,j).ny ~= 0)
                    set1 = cell(kk,j); 
                            
                    %@@@ Like in x obstruction, remove boundaries 
                    %@@@ in shadow zone and combine segments
                           
                    shadow_flags = zeros(size(baseseg_n));
                    for m = 1:n_segs
                        for l = 1:set1.ny
                            if (set1.east_lim(l) >= baseseg_n(m) && ...
                                    set1.west_lim(l) <= baseseg_s(m))
                                shadow_flags(m) = 1;
                                break;
                            end;
                        end;                                
                    end;

                    loc = find(shadow_flags == 0);

                    if (isempty(loc))
                        no_boundary = 1;
                        n_segs = 0;
                        baseseg_n = [];
                        baseseg_s = [];
                    elseif (length(loc) < n_segs)
                        tmp_n = baseseg_n;
                        tmp_s = baseseg_s;
                        baseseg_n = [];
                        baseseg_s = [];
                        baseseg_n = tmp_n(loc);
                        baseseg_s = tmp_s(loc);
                        n_segs = length(baseseg_n);
                        clear tmp_n tmp_s;
                    end;

                    clear loc;
                    clear shadow_flags
                            
                    if (~no_boundary)
                        shadow_flags = zeros(size(set1.east_lim));
                        for m = 1:set1.ny
                            for l = 1:n_segs
                                if (set1.east_lim(m) <= baseseg_n(l) && ...
                                        set1.west_lim(m) >= baseseg_s(l))
                                    shadow_flags(m) = 1;
                                    break;
                                end;
                            end;                                
                        end;
                        loc = find(shadow_flags == 0);
                        if (~isempty(loc))
                            if (length(loc) < set1.ny)
                                tmp_n = set1.east_lim;
                                tmp_s = set1.west_lim;
                                set1.east_lim = [];
                                set1.west_lim = [];
                                set1.east_lim = tmp_n(loc);
                                set1.west_lim = tmp_s(loc);
                                set1.ny = length(set1.east_lim);
                                clear tmp_n;
                                clear tmp_s;
                            end;
                            n_segs = n_segs+set1.ny;
                            baseseg_n = [baseseg_n set1.east_lim];
                            baseseg_s = [baseseg_s set1.west_lim];
                        end;
                        clear loc;
                        clear shadow_flags;
                    end;                          
                end;
            end;
        end;  %@@@ End of accounting for boundaries in cell below

        %@@@ Now moving to boundaries in the cell above (if applicable)
                
        if (~no_boundary)
            
            for off = 1:offset_right
                kk = k+off;
                if (kk <= Ny)
                    if (cell(kk,j).ny ~= 0)
                        set1 = cell(kk,j); 

                        %@@@ Again follow the same procedure to merge boundaries 
                                
                        shadow_flags = zeros(size(baseseg_n));
                        for m = 1:n_segs
                            for l = 1:set1.ny
                                if (set1.east_lim(l) >= baseseg_n(m) && ...
                                        set1.west_lim(l) <= baseseg_s(m))
                                    shadow_flags(m) = 1;
                                    break;
                                end;
                            end;                                
                        end;

                        loc = find(shadow_flags == 0);
                        if (isempty(loc))
                            no_boundary = 1;
                            n_segs = 0;
                            baseseg_n = [];
                            baseseg_s = [];
                        elseif (length(loc) < n_segs)
                            tmp_n = baseseg_n;
                            tmp_s = baseseg_s;
                            baseseg_n = [];
                            baseseg_s = [];
                            baseseg_n = tmp_n(loc);
                            baseseg_s = tmp_s(loc);
                            n_segs = length(baseseg_n);
                            clear tmp_n tmp_s;
                        end;
                        clear loc shadow_flags;
                                
                        if (~no_boundary)
                            shadow_flags = zeros(size(set1.east_lim));
                            for m = 1:set1.ny
                                for l = 1:n_segs
                                    if (set1.east_lim(m) <= baseseg_n(l) && ...
                                            set1.west_lim(m) >= baseseg_s(l))
                                        shadow_flags(m) = 1;
                                        break;
                                    end;
                                end;                                
                            end;

                            loc = find(shadow_flags == 0);
                            if (~isempty(loc))
                                if (length(loc) < set1.ny)
                                    tmp_n = set1.east_lim;
                                    tmp_s = set1.west_lim;
                                    set1.east_lim = [];
                                    set1.west_lim = [];
                                    set1.east_lim = tmp_n(loc);
                                    set1.west_lim = tmp_s(loc);
                                    set1.ny = length(set1.east_lim);
                                    clear tmp_n tmp_s;
                                end;
                                n_segs = n_segs+set1.ny;
                                baseseg_n = [baseseg_n set1.east_lim];
                                baseseg_s = [baseseg_s set1.west_lim];
                            end;
                            clear loc shadow_flags;

                        end;                           
                    end;
                end;
            end;
        end;  %@@@ End of accounting for cell above
                
        %@@@ Building obstruction grid in y  
                
        if (~no_boundary)

            if (n_segs==1)
            
                sy(k,j) = baseseg_n(1)-baseseg_s(1);
            
            else    

                ind_segs=0;
                indseg_n = [];
                indseg_s = [];

                while (n_segs > 0)
                    overlap_found = 0;
                    if (n_segs > 1)
                        for l = 2:n_segs
                            if (baseseg_n(1) >= baseseg_s(l) && baseseg_s(1) ...
                                    <= baseseg_n(l))
                                baseseg_n(1) = max([baseseg_n(1) baseseg_n(l)]);
                                baseseg_s(1) = min([baseseg_s(1) baseseg_s(l)]);
                                overlap_found = 1;
                                if (l == n_segs)
                                    n_segs = n_segs-1;
                                else
                                    for m = l+1:n_segs
                                        baseseg_n(m-1) = baseseg_n(m);
                                        baseseg_s(m-1) = baseseg_s(m);
                                    end;
                                    n_segs = n_segs-1;
                                end;
                                break;
                            end;
                        end;
                    end;
                        
                    if (n_segs == 1)
                        ind_segs = ind_segs+1;
                        indseg_n(ind_segs) = baseseg_n(1);
                        indseg_s(ind_segs) = baseseg_s(1);
                        n_segs = n_segs-1;
                    else
                        if (overlap_found == 0)
                            ind_segs = ind_segs+1;
                            indseg_n(ind_segs) = baseseg_n(1);
                            indseg_s(ind_segs) = baseseg_s(1);
                            for l = 2:n_segs                                    
                                baseseg_n(l-1) = baseseg_n(l);
                                baseseg_s(l-1) = baseseg_s(l);
                            end;
                            n_segs = n_segs-1;
                        end;
                    end;
                end;
      
                for l = 1:ind_segs
                    sy(k,j) = sy(k,j) + indseg_n(l)-indseg_s(l);
                end;
                cell(k,j).indsegs_east = indseg_n;
                cell(k,j).indsegs_west = indseg_s;
                clear baseseg_n baseseg_s baseseg_bound ;
            end;
        end;
        
    end; %@@@ End of obstruction grid in y                                 

    %@@@ Setting the obstruction grid to zero if neighboring cells are
    %@@@ dry cells to prevent spurious swell attenuation near the coast 

     if (j < Nx && mask(k,j+1) == 0)
        sx(k,j) = 0;
    end;
    if (j > 1 && mask(k,j-1) == 0)
        sx(k,j) = 0;
    end;
    if (k < Ny && mask(k+1,j) == 0)
        sy(k,j) = 0;
    end;
    if (k > 1 && mask(k-1,j) == 0)
        sy(k,j) = 0;
    end;

end; %@@@ End of wet cell (with enclosed boundaries) loop

return;
