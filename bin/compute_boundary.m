function [bound_ingrid,Nb] = compute_boundary(coord,bound,varargin)

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
%| Computes the shoreline polygons from the GSHHS database that lie within| 
%| the grid domain, properly accounting for polygons that cross the domain| 
%| The routine has been designed to work with coastal polygons by default.|
%| That can be changed by using a different boundary flag. See GSHHS      |
%| documentation (or below) for the meaning of the different flags        |
%|                                                                        |
%| [bound_ingrid,Nb] = compute_boundary(coord,bound,[bflg])               |
%|                                                                        |
%| INPUT                                                                  |
%|   coord : An array defining the corner points of the grid              |
%|           coord(1) = Lattitude (y) of lower left hand corner           |
%|           coord(2) = Longitude (x) of lower left hand corner           |
%|           coord(3) = Lattitude (y) of upper right hand corner          |
%|           coord(4) = Longitude (x) of upper right hand corner          |
%|   bound : A data structure array of the basic polygons (The GSHHS      |
%|           polygons are stored as mat files with several different      | 
%|           resolutions and the user should ensure that the files have   |
%|           been loaded before using this routine).                      |
%|                                                                        |
%|           The different available files are --                         |
%|             coastal_bound_ful.mat    -- Full resolution                |
%|                                         (188606 polygons)              |
%|		       coastal_bound_high.mat   -- High resolution        |
%|                                         (0.2 km; 153539 polygons)      |
%|		       coastal_bound_inter.mat  -- Intermediate resolution|
%|                                         (1 km; 41523 polygons)         |
%|		       coastal_bound_low.mat    -- Low resolution         |
%|                                         (5 km; 10769 polygons)         |
%|		       coastal_bound_coarse.mat -- Coarse resolution      |
%|                                         (25 km; 1866 polygons)         |
%|                                                                        |
%|    bflg : Optional definition of flag type from the gshhs boundary     |
%|           database (1 = land; 2 = lake margin; 3 = in-lake island).    |
%|           If left blank, defaults to land (1).                         |
%|                                                                        |
%|	    Alternatively, a separate list of user defined polygons can   |
%|          be generated having the same fields as bound. One such list is|
%|	    "optional_coastal_polygons.mat" which is also distributed with|
%|          reference data. This is an ever growing list of water bodies  |
%|          that see very little wave action and for most practical       | 
%|          purposes can be masked out as land.                           |
%|                                                                        | 
%|          See also optional_bound.m which shows how the optional coastal|
%|          polygons are used                                             |
%|                                                                        | 
%| OUTPUT                                                                 |
%|  bound_ingrid : Subset data structure array of polygons that lie inside|
%|                 the grid                                               |
%|  Nb           : Total number of polygons found that lie inside the grid|
%|                                                                        |
% -------------------------------------------------------------------------
narg=nargin;

% Determine if third input variable present (requesting inland features)
if narg == 3
  bflg = varargin{1}; % Highest level admitted to provide shoreline data
elseif narg == 2
  bflg = 1; % Normal application only continental bounds
else
  disp('Too many input arguments, exiting')
  return;
end  

lat_start = coord(1);
lon_start = coord(2);
lat_end = coord(3);
lon_end = coord(4);

%@@@ Definitions

%@@@ Minimum distance between points (to avoid round off errors from points
%@@@ too close to each other)

eps = 1e-5;

%@@@ Maximum distance between points. This is being defined so that we do
%@@@ not have large gaps between subsequent points of the final boundary
%@@@ polygon

MAX_SEG_LENGTH = 0.25;

%@@@ Polygon defining the bounding grid. Bounding grid is defined in the 
%@@@ counter clockwise direction

px = [lon_start lon_end lon_end lon_start lon_start];
py = [lat_start lat_start lat_end lat_end lat_start];

%@@@ Slope and intercepts for each of the 4 lines of the bounding box
 
for i = 1:4
    if (px(i+1)==px(i))
        m_grid(i)=inf;
        c_grid(i)=0;
    else
        p = polyfit(px(i:i+1),py(i:i+1),1);
        m_grid(i) = p(1);
        c_grid(i) = p(2);
    end;
    box_length(i) = sqrt((px(i+1)-px(i))^2+(py(i+1)-py(i))^2);
    norm(i,1) = (py(i+1)-py(i))./box_length(i);
    norm(i,2) = -(px(i+1)-px(i))./box_length(i);
end;
norm(end+1,1) = norm(1,1);
norm(end+1,2) = norm(1,2);

%@@@ Initializing variables

N = length(bound);
in_coord = 1;
itmp = 0;

%@@@ Loop through all the boundaries in the database

for i = 1:N
    
    %@@@ Limit boundaries to coastal type only. This flag needs to be 
    %@@@ changed if interested in other boundaries. See GSHHS documentation 
    %@@@ for boundary type flags
   
    if (bound(i).level == bflg )
                                                     
        %@@@ Determine if boundary lies completely outside the domain
        
        if (bound(i).west > lon_end || bound(i).east < lon_start ...     
              || bound(i).south > lat_end || bound(i).north < lat_start) 
            in_grid = 0;                                                 
        else
            in_grid = 1;
        end;

        lev1 = bound(i).level;

        %@@@ Determine if boundary lies completely inside the domain

        if (bound(i).west >= lon_start && bound(i).east <= lon_end && ...
               bound(i).south >= lat_start && bound(i).north <= lat_end)
            inside_grid = 1;
        else
            inside_grid = 0;
        end;

        %@@@ Ignore boundaries outside the domain

        if (in_grid)

            %@@@ Modify boundaries that are not completely inside domain
 
            if (~inside_grid)

                %@@@ Determine the points of the boundary that are 
                %@@@ inside/on/outside the bounding box
 
                [in_points,on_points] = inpolygon(bound(i).x,bound(i).y,px,py);
                
                loc1 = find(in_points == 1);
                loc2 = find(on_points == 1);

                %@@@ Ignore points that lie on the domain but neighboring 
                %@@@ points do not
                
                for j = 1:length(loc2)
                    if (loc2(j) == 1)
                        p1 = bound(i).n;
                        p2 = loc2(j)+1;
                    elseif (loc2(j) == bound(i).n)
                        p1 = loc2(j)-1;
                        p2 = 1;
                    else
                        p1 = loc2(j)-1;
                        p2 = loc2(j)+1;
                    end;
                    if (in_points(p1) == 0 && in_points(p2) == 0)
                        in_points(loc2(j)) = 0;
                    end;
                end;                

                %@@@ Points of domain in the boundary
                
                domain_inb = inpolygon(px,py,bound(i).x,bound(i).y);
                loc_t = find(domain_inb == 1);
                domain_inb_lth = length(loc_t);

                if (isempty(loc1))
                    if (domain_inb_lth == length(px))
                        bound_ingrid(in_coord).x = px;
                        bound_ingrid(in_coord).y = py; 
                        bound_ingrid(in_coord).n = length(px);
                        bound_ingrid(in_coord).west = lon_start;
                        bound_ingrid(in_coord).east = lon_end;
                        bound_ingrid(in_coord).north = lat_end;
                        bound_ingrid(in_coord).south = lat_start;
                        bound_ingrid(in_coord).height = lon_end - lon_start;
                        bound_ingrid(in_coord).width = lat_end - lat_start;
                        bound_ingrid(in_coord).level = lev1;
                        in_coord = in_coord + 1;
                    end;
                    clear loc_t domain_inb;
                end;
                
                %@@@ Loop through only if there are points inside the domain
                
                if (~isempty(loc1)) 
                
                    n = bound(i).n;

                    %@@@ Flag the points where the boundary moves from in 
                    %@@@ to out of the domain as well as out to in

                    in2out_count = 1;
                    out2in_count = 1;
                
                    out2in = [];
                    in2out = [];

                    for j = 1:bound(i).n-1
                        if (in_points(j) > 0 && in_points(j+1) == 0)
                            in2out(in2out_count) = j;
                            in2out_count = in2out_count+1;
                        end;
                        if (in_points(j) == 0 && in_points(j+1) > 0)
                            out2in(out2in_count) = j;
                            out2in_count = out2in_count+1;
                        end;
                    end;

                    in2out_count = in2out_count-1;
                    out2in_count = out2in_count-1;
                    if (in2out_count ~= out2in_count)
                        fprintf(1,'Error: mismatch in grid crossings, check boundary %d !! \n',i);
                        return;
                    end;                
               
                    %@@@ Crossing points are oriented to make sure we start 
                    %@@@ from out to in

                    if (in_points(1) > 0)
                        in2out_tmp = in2out;
                        for j = 1:in2out_count-1
                            in2out(j) = in2out_tmp(j+1);
                        end;
                        in2out(in2out_count) = in2out_tmp(1);
                    end;
                    
                    clear in2out_tmp;
                                  
                    %@@@ For each in2out and out2in find a grid intersecting point
                                    
                    clear in2out_gridbox in2out_gridboxdist in2out_xcross in2out_ycross;
                    clear out2in_gridbox out2in_gridboxdist out2in_xcross out2in_ycross;
                    
                    in2out_gridbox = zeros(in2out_count,1);
                    out2in_gridbox = zeros(out2in_count,1);
                    in2out_gridboxdist = zeros(in2out_count,1);
                    out2in_gridboxdist = zeros(out2in_count,1);
                
                    for j = 1:in2out_count

                        if(on_points(in2out(j)) == 1)

                            x1 = bound(i).x(in2out(j));
                            y1 = bound(i).y(in2out(j));

                            for k = 1:4

                                if (isinf(m_grid(k)))
                                    g = abs(x1-px(k));
                                else
                                    g = abs(m_grid(k)*x1+c_grid(k)-y1);
                                end;

                                if (g <= eps)
                                   in2out_gridbox(j) = k;
                                   in2out_gridboxdist(j) = k-1 + ...
                                       sqrt((px(k)-x1)^2+(py(k)-y1)^2)/box_length(k);
                                   break;
                                end;

                            end;

                            in2out_xcross(j) = NaN;
                            in2out_ycross(j) = NaN;

                        else

                            x1 = bound(i).x(in2out(j));
                            x2 = bound(i).x(in2out(j)+1);
                            y1 = bound(i).y(in2out(j));
                            y2 = bound(i).y(in2out(j)+1);
                        
                            if (x2==x1)
                                m = inf;
                                c = 0;
                            else
                                if (abs(x2-x1) > 90)
                                    x2 = x2+360;
                                    if (abs(x2-x1) > 90)
                                        x2 = x2-720;
                                    end;
                                end;
                                p = polyfit([x1 x2],[y1 y2],1);
                                m = p(1);
                                c = p(2);
                            end;
                            d = sqrt((x1-x2)^2+(y1-y2)^2);
           
                            for k = 1:4
                                if (m ~= m_grid(k));
                                    if (~isinf(m) && ~isinf(m_grid(k)))
                                        x = (c_grid(k)-c)/(m-m_grid(k));
                                        y = m*x+c;
                                    elseif (isinf(m))
                                        x = x1;
                                        y = m_grid(k)*x+c_grid(k);
                                    else
                                        x = px(k);
                                        y = m*x+c;
                                    end;
                                
                                    d1 = sqrt((x1-x)^2+(y1-y)^2);
                                    d2 = sqrt((x-x2)^2+(y-y2)^2);
                                    if (abs(1-(d1+d2)/d) < 0.001)
                                        in2out_gridbox(j) = k;
                                        break;
                                    end;
                                end;
                            end;
                        
                            in2out_xcross(j) = x;
                            in2out_ycross(j) = y;
                            in2out_gridboxdist(j) = k-1 + ...
                                sqrt((px(k)-x)^2+(py(k)-y)^2)/box_length(k);

                        end; %@@@ corresponds to if(on_points(in2out(j)) == 1)
                        
                        if (on_points(out2in(j)+1) == 1)

                            x1 = bound(i).x(out2in(j)+1);
                            y1 = bound(i).y(out2in(j)+1);

                            for k = 1:4
                                if (isinf(m_grid(k)))
                                    g = abs(x1-px(k));
                                else
                                    g = abs(m_grid(k)*x1+c_grid(k)-y1);
                                end;
                                if (g <= eps)
                                   out2in_gridbox(j) = k;
                                   out2in_gridboxdist(j) = k-1 + ...
                                       sqrt((px(k)-x1)^2+(py(k)-y1)^2)/box_length(k);
                                   break;
                                end;
                            end;

                            out2in_xcross(j) = NaN;
                            out2in_ycross(j) = NaN;

                        else

                            x1 = bound(i).x(out2in(j));
                            x2 = bound(i).x(out2in(j)+1);
                            y1 = bound(i).y(out2in(j));
                            y2 = bound(i).y(out2in(j)+1);
                        
                            if (x2==x1)
                                m = inf;
                                c = 0;
                            else
                                if (abs(x2-x1) > 90)
                                    x1 = x1+360;
                                    if (abs(x2-x1) > 90)
                                        x1 = x1-720;
                                    end;
                                end;
                                p = polyfit([x1 x2],[y1 y2],1);
                                m = p(1);
                                c = p(2);
                            end;
                            d = sqrt((x1-x2)^2+(y1-y2)^2);
           
                            for k = 1:4
                                if (m ~= m_grid(k))
                                    if (~isinf(m) && ~isinf(m_grid(k)))
                                        x = (c_grid(k)-c)/(m-m_grid(k));
                                        y = m*x+c;
                                    elseif (isinf(m))
                                        x = x1;
                                        y = m_grid(k)*x+c_grid(k);
                                    else
                                        x = px(k);
                                        y = m*x+c;
                                    end;
                                    d1 = sqrt((x1-x)^2+(y1-y)^2);
                                    d2 = sqrt((x-x2)^2+(y-y2)^2);
                                    if (abs(1-(d1+d2)/d) < 0.001)
                                        out2in_gridbox(j) = k;
                                        break;
                                    end;
                                end;
                            end;
                            out2in_xcross(j) = x;
                            out2in_ycross(j) = y;
                            out2in_gridboxdist(j) = k-1 + ...
                                sqrt((px(k)-x)^2+(py(k)-y)^2)/box_length(k);

                        end; %@@@ corresponds to if(on_points(out2in(j)) == 1)

                    end;     %@@@ end of j loop for all the intersection points
                
                    %@@@ Loop through the intersection points 

                    if (in2out_count > 0)

                        subseg_acc = zeros(in2out_count,1);
                        crnr_acc = 1 - domain_inb;
     
                        while(~isempty(find(subseg_acc == 0,1)))
                            
                            %@@@ Starting from the closest unaccounted
                            %@@@ segment
                            
                            min_pos = 0;
                            min_val = 4;
                            for j = 1:in2out_count
                                if (subseg_acc(j) == 0)
                                    if (out2in_gridboxdist(j) <= min_val)
                                        min_val = out2in_gridboxdist(j);
                                        min_pos = j;
                                    end;
                                end;
                            end;
                            
                            j = min_pos;   
                            
                            bound_ingrid(in_coord).x = [];
                            bound_ingrid(in_coord).y = [];
                            bound_ingrid(in_coord).n = 0;
                            bound_ingrid(in_coord).east = 0;
                            bound_ingrid(in_coord).west = 0;
                            bound_ingrid(in_coord).north = 0;
                            bound_ingrid(in_coord).south = 0;
                            bound_ingrid(in_coord).height = 0;
                            bound_ingrid(in_coord).width = 0;
                            
                            bound_x = [];
                            bound_y = [];

                            if (~isnan(out2in_xcross(j)))
                                bound_x = out2in_xcross(j);
                                bound_y = out2in_ycross(j);
                            end;

                            if ((out2in(j)+1) <= in2out(j))
                                bound_x = [bound_x;bound(i).x((out2in(j)+1):in2out(j))];
                                bound_y = [bound_y;bound(i).y((out2in(j)+1):in2out(j))];                       
                            else                               
                                bound_x = [bound_x;bound(i).x((out2in(j)+1):n);...
                                                  bound(i).x(2:in2out(j))];
                                bound_y = [bound_y;bound(i).y((out2in(j)+1):n);...
                                                  bound(i).y(2:in2out(j))];                              
                            end;

                            if (~isnan(in2out_xcross(j)))
                                bound_x = [bound_x;in2out_xcross(j)];
                                bound_y = [bound_y;in2out_ycross(j)];
                            end;

                            close_bound=0; %@@@ Flag initializing close boundary
                                
                            starting_edge = out2in_gridbox(j);
                            ending_edge = in2out_gridbox(j);
                                
                            subseg_acc(j) = 1;
                            
                            %@@@ Find the next closest segment going
                            %@@@ anti-clockwise
                            
                            seg_index = j;                           
                            
                            
                            while (close_bound == 0)
                                
                                %@@@ Check if last segment and see if can
                                %@@@ proceed counter clockwise 
                            
                                if (isempty(find(subseg_acc == 0,1)))
                                    for k = in2out_gridbox(seg_index):4
                                        if (domain_inb(k+1) == 1 && crnr_acc(k+1) == 0)
                                            bound_x = [bound_x;px(k+1)];
                                            bound_y = [bound_y;py(k+1)];
                                            crnr_acc(k+1) = 1;
                                            ending_edge = k;
                                        else
                                            close_bound = 1;
                                            break;
                                        end;
                                    end;
                                
                                    if (close_bound == 0)
                                        for k = 1:(in2out_gridbox(seg_index)-1)
                                            if (domain_inb(k+1) == 1 && crnr_acc(k+1) == 0)
                                                bound_x = [bound_x;px(k+1)];
                                                bound_y = [bound_y;py(k+1)];
                                                crnr_acc(k+1) = 1;
                                                ending_edge = k;
                                            else
                                                close_bound = 1;
                                                break;
                                            end;
                                        end;
                                        close_bound = 1;
                                    end;
                                
                                else
                                                                                                    
                                    curr_seg = seg_index;                           
                                    kstart = in2out_gridbox(curr_seg);
                                    start_dist = in2out_gridboxdist(curr_seg);
                                    min_pos = 0;
                                    min_val = 4.0;
                                
                                    %@@@ Check all segments

                                    for k1 = 1:in2out_count                                        
                                        if ((out2in_gridboxdist(k1)-start_dist) > eps ... 
                                            && (out2in_gridboxdist(k1)-start_dist) < min_val)
                                            min_pos = k1;
                                            min_val = out2in_gridboxdist(k1)-start_dist;
                                        end;                                        
                                    end;
                                
                                    if (min_pos == 0)              
                                        for k1 = 1:out2in_count
                                            if (out2in_gridboxdist(k1) < min_val)
                                                min_pos = k1;
                                                min_val = out2in_gridboxdist(k1);
                                            end;
                                        end;                                                                                    
                                    end;
                            
                                    if (subseg_acc(min_pos) == 1)                                                              
                                        close_bound = 1;          
                                        ending_edge = in2out_gridbox(curr_seg);
                                    else 
                                        kend = out2in_gridbox(min_pos);
                                        x_mid = [];
                                        y_mid = [];

                                        %@@@ If the boundary polygon crosses 
                                        %@@@ the grid domain along different
                                        %@@@ domain edges then include the 
                                        %@@@ common grid domain corner points

                                        
                                        if (kstart ~= kend)
                                            if (kstart < kend)
                                                for k1 = kstart:(kend-1)
                                                    if (domain_inb(k1+1) == 1 ...
                                                            && crnr_acc(k1+1) == 0)
                                                        x_mid = [x_mid;px(k1+1)];
                                                        y_mid = [y_mid;py(k1+1)];
                                                        crnr_acc(k1+1) = 1;
                                                    end;
                                                end;
                                            else
                                                for k1 = kstart:4
                                                    if (domain_inb(k1+1) == 1 ...
                                                            && crnr_acc(k1+1) == 0)
                                                        x_mid = [x_mid;px(k1+1)];
                                                        y_mid = [y_mid;py(k1+1)];
                                                        crnr_acc(k1+1) = 1;
                                                    end;
                                                end;
                                                for k1 = 1:(kend-1)
                                                    if (domain_inb(k1+1) == 1 ...
                                                            && crnr_acc(k1+1) == 0)
                                                        x_mid = [x_mid;px(k1+1)];
                                                        y_mid = [y_mid;py(k1+1)];
                                                        crnr_acc(k1+1) = 1;
                                                    end;
                                                end;
                                            end;
                                        end;
                                    
                                        if (~isempty(x_mid))
                                            bound_x = [bound_x;x_mid];
                                            bound_y = [bound_y;y_mid];
                                        end;
                                    
                                        %@@@ Adding the segment
                                    
                                        if (~isnan(out2in_xcross(min_pos)))
                                            bound_x = [bound_x;...
                                                           out2in_xcross(min_pos)];
                                            bound_y = [bound_y;...
                                                           out2in_ycross(min_pos)];
                                        end;

                                        if ((out2in(min_pos)+1) <= in2out(min_pos))                                        
                                            bound_x = [bound_x;...
                                                           bound(i).x((out2in(min_pos)+1):...
                                                                        in2out(min_pos))];
                                            bound_y = [bound_y;...
                                                           bound(i).y((out2in(min_pos)+1):...
                                                                        in2out(min_pos))];                                        
                                        else
                                            bound_x = [bound_x;...
                                                           bound(i).x((out2in(min_pos)+1):n);...
                                                           bound(i).x(2:in2out(min_pos))];
                                            bound_y = [bound_y;...
                                                           bound(i).y((out2in(min_pos)+1):n);...
                                                           bound(i).y(2:in2out(min_pos))];      
                                        end;

                                        if (~isnan(in2out_xcross(min_pos)))
                                            bound_x = [bound_x;in2out_xcross(min_pos)];
                                            bound_y = [bound_y;in2out_ycross(min_pos)];
                                        end;

                                        subseg_acc(min_pos) = 1;
                                        ending_edge = in2out_gridbox(min_pos);
                                        seg_index = min_pos;                             
                                    end;
                                    
                                end;
                                
                            end;
                            
                            %@@@ Need to close the grid;
                            
                            if (ending_edge ~= starting_edge)
                                if (ending_edge < starting_edge)
                                    for k = ending_edge:(starting_edge-1)
                                        if (crnr_acc(k+1) == 0 && ...
                                                domain_inb(k+1) == 1)
                                            bound_x = [bound_x;px(k+1)];
                                            bound_y = [bound_y;py(k+1)];
                                            crnr_acc(k+1) = 1;
                                        end;
                                    end;
                                else
                                    for k = ending_edge:4
                                        if (crnr_acc(k+1) == 0 && ...
                                                domain_inb(k+1) == 1)
                                            bound_x = [bound_x;px(k+1)];
                                            bound_y = [bound_y;py(k+1)];
                                            crnr_acc(k+1) = 1;
                                        end;
                                    end;
                                    for k =1:(starting_edge-1)
                                       if (crnr_acc(k+1) == 0 && ...
                                               domain_inb(k+1) == 1)
                                            bound_x = [bound_x;px(k+1)];
                                            bound_y = [bound_y;py(k+1)];
                                            crnr_acc(k+1) = 1;
                                        end;
                                    end;
                                end;
                            end;
                            
                            bound_x(end+1) = bound_x(1);
                            bound_y(end+1) = bound_y(1);
                            
                            %@@@ Making sure that the added points do not
                            %@@@ exceed max. defined seg length
                            
                            clear xt1 xt2 yt1 yt2 dist loc x_set y_set;
                            nsample = length(bound_x);
                            xt1 = bound_x(1:end-1);
                            yt1 = bound_y(1:end-1);
                            xt2 = bound_x(2:end);
                            yt2 = bound_y(2:end);
                            
                            dist = sqrt((xt2-xt1).^2+(yt2-yt1).^2);
                            loc = find(dist > 2*MAX_SEG_LENGTH);
                            
                            if (~isempty(loc))
                                x_set = bound_x(1:loc(1));
                                y_set = bound_y(1:loc(1));
                                nc = length(loc);
                                for k = 1:nc
                                    xp = bound_x(loc(k));
                                    yp = bound_y(loc(k));
                                    xn = bound_x(loc(k)+1);
                                    yn = bound_y(loc(k)+1);
                                    ns = floor(dist(loc(k))/MAX_SEG_LENGTH)-1;
                                    if (ns > 0)
                                        if (xp == xn)
                                            x_set = [x_set;ones(ns,1)*xp];
                                            y_set = [y_set;(yp+sign(yn-yp)...
                                                *(1:ns)'*MAX_SEG_LENGTH)];
                                        else
                                            mth = atan2(yn-yp,xn-xp);
                                            x_set = [x_set;(xp+[1:ns]'...
                                                *MAX_SEG_LENGTH*cos(mth))];
                                            y_set = [y_set;(yp+[1:ns]'...
                                                *MAX_SEG_LENGTH*sin(mth))];
                                        end;
                                    end;
                                    x_set = [x_set;xn];
                                    y_set = [y_set;yn];
                                    if k == nc
                                        if ((loc(k)+1) < nsample)
                                            x_set = [x_set;bound_x((loc(k)+2:end))];
                                            y_set = [y_set;bound_y((loc(k)+2:end))];
                                        end;
                                    else
                                        if ((loc(k)+1) < loc(k+1))
                                            x_set = [x_set;bound_x((loc(k)+2):loc(k+1))];
                                            y_set = [y_set;bound_y((loc(k)+2):loc(k+1))];
                                        end;
                                    end;
                                end;
                            else
                                x_set = bound_x;
                                y_set = bound_y;
                            end;
                            
                            %@@@ Setting up the boundary polygon
                            
                            bound_ingrid(in_coord).x = x_set;
                            bound_ingrid(in_coord).y = y_set;
                            bound_count = length(bound_ingrid(in_coord).x);
                            bound_ingrid(in_coord).n = bound_count;                                 
                            bound_ingrid(in_coord).east = max(bound_ingrid(in_coord).x);
                            bound_ingrid(in_coord).west = min(bound_ingrid(in_coord).x);
                            bound_ingrid(in_coord).north = max(bound_ingrid(in_coord).y);
                            bound_ingrid(in_coord).south = min(bound_ingrid(in_coord).y);
                            bound_ingrid(in_coord).height = bound_ingrid(in_coord).north ...
                                            - bound_ingrid(in_coord).south;
                            bound_ingrid(in_coord).width = bound_ingrid(in_coord).east ...
                                            - bound_ingrid(in_coord).west;
                            bound_ingrid(in_coord).level = lev1;
                                                       
                            in_coord=in_coord+1;    %@@@ increment boundary 
                                                    %@@@ counter                                                                       

                            crnr_acc(1) = crnr_acc(end);

                        end;         %@@@ corresponds to while loop that 
                                     %@@@ checks if all sections (subseg_acc) 
                                     %@@@ have been accounted for. 
                   
                    end;             %@@@ corresponds to if in2out_count > 0 
                                   
                end;                 %@@@ corresponds to if statement checking 
                                     %@@@ if there are boundary points inside 
                                     %@@@ the domain

            else                     %@@@ boundary lies completely inside the grid

                %@@@ initializing and adding the boundary to the list
                bound_ingrid(in_coord).x = [];
                bound_ingrid(in_coord).y = [];
                bound_ingrid(in_coord).n = 0;
                bound_ingrid(in_coord).east = 0;
                bound_ingrid(in_coord).west = 0;
                bound_ingrid(in_coord).north = 0;     
                bound_ingrid(in_coord).south = 0;
                bound_ingrid(in_coord).height = 0;
                bound_ingrid(in_coord).width = 0;
                bound_ingrid(in_coord).n = bound(i).n;
                bound_ingrid(in_coord).x = bound(i).x;
                bound_ingrid(in_coord).y = bound(i).y;
                bound_ingrid(in_coord).east = max(bound_ingrid(in_coord).x);
                bound_ingrid(in_coord).west = min(bound_ingrid(in_coord).x);
                bound_ingrid(in_coord).north = max(bound_ingrid(in_coord).y);
                bound_ingrid(in_coord).south = min(bound_ingrid(in_coord).y);
                bound_ingrid(in_coord).height = bound_ingrid(in_coord).north ...
                    - bound_ingrid(in_coord).south;
                bound_ingrid(in_coord).width = bound_ingrid(in_coord).east ...
                    - bound_ingrid(in_coord).west;
                bound_ingrid(in_coord).level = lev1;
                in_coord = in_coord+1;

            end;     %@@@ corresponds to if statement that determines if boundary 
                     %@@@ lies partially/completely in domain

        end;         %@@@ corresponds to if statement that determines if boundary  
                     %@@@ lies outside the domain

    end;             %@@@ corresponds to if statement that determines boundary type

    %@@@ counter to keep tab on the level of processing

    itmp_prev = itmp;
    itmp = floor(i/N*100);
    if (mod(itmp,5)==0 && itmp_prev ~= itmp && N > 100)
        fprintf(1,'Completed %d per cent of %d boundaries and found %d internal boundaries \n',...
            itmp,N,in_coord-1);
    end;

end;     %@@@ end of for loop that loops through all the GSHHS boundaries

Nb = in_coord-1;

if (Nb == 0)
    bound_ingrid(1) = -1;
end;

return;
