function [c1,c2,c3,c4,wdth,hgt] = compute_cellcorner(x,y,j,k,Nx,Ny)

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
%| This function determines the corners of a particular cell at the jth   |
%| row and kth column, given the 2D position matrices x and y. The        |
%| function returns 4 2 element arrays where the first element refers to  |
%| the x coordinate and the second element to the y coordinate. The       |
%| corners are defined as shown below (* indicating the position of the   |
%| cell in x,y space). The orientation of the cell will be determined by  |
%| the neighboring x,y coordinates                                        |
%|                                                                        |
%|                   c3 _____ c2                                          |
%|                     |     |                                            |
%|                     |  *  |                                            |
%|                   c4|_____|c1                                          |
%|                                                                        |
%| [c1,c2,c3,c4] = compute_cellcorner(x,y,j,k,Nx,Ny)                      | 
%|                                                                        |
%| INPUT                                                                  |
%|  x            : A 2D array specifying the longitudes of each cell      | 
%|  y            : A 2D array specifying the lattitudes of each cell      |
%|  j,k          : jth column and kth row of the 2D array (x,y)           |
%|  Nx,Ny        : Number of columns and rows of the 2D arrays            |
%|                                                                        |
%| OUTPUT                                                                 |
%|  c1,c2,c3,c4  : Corners of the cell as shown above.In each first       |
%|                 element is x coordinate and second element is y        |
%|                 coordinate                                             |
%|  wdth, hgt    : Cell width and height                                  |
%|                                                                        |
% -------------------------------------------------------------------------

%@@@ Initialize variables

c1 = [];
c2 = [];
c3 = [];
c4 = [];

x0 = x(k,j);
%@@@ Compute the corners of the cell

if ( j > 1 && j < Nx && k > 1 && k < Ny )
    
    %@@@ Internal points
    
    xt = x(k-1,j+1);
    if (abs(xt - x0) > 270)
        xt = xt - 360*sign(xt-x0);
    end;
    c1(1) = 0.5*(xt+x0);
    c1(2) = 0.5*(y(k-1,j+1) + y(k,j));
    xt = x(k+1,j+1);
    if (abs(xt - x0) > 270)
        xt = xt - 360*sign(xt-x0);
    end;
    c2(1) = 0.5*(xt+x0);
    c2(2) = 0.5*(y(k+1,j+1) + y(k,j));
    xt = x(k+1,j-1);
    if (abs(xt - x0) > 270)
        xt = xt - 360*sign(xt-x0);
    end;
    c3(1) = 0.5*(xt+x0);
    c3(2) = 0.5*(y(k+1,j-1) + y(k,j));
    xt = x(k-1,j-1);
    if (abs(xt - x0) > 270)
        xt = xt - 360*sign(xt-x0);
    end;
    c4(1) = 0.5*(xt+x0);
    c4(2) = 0.5*(y(k-1,j-1) + y(k,j));
    
    
else
    
    if ( j == 1 )
        
        %@@@ Left edge
        
        switch k
            case 1
                xt = x(k+1,j+1);
                if (abs(xt - x0) > 270)
                    xt = xt - 360*sign(xt-x0);
                end;
                c2(1) = 0.5*(xt+x0);
                c2(2) = 0.5*(y(k+1,j+1) + y(k,j));
                c4(1) = 2*x0 - c2(1);
                c4(2) = 2*y(k,j) - c2(2);
                c3(1) = x0 - (c2(2)-y(k,j));
                c3(2) = y(k,j) + (c2(1)-x(k,j));
                c1(1) = 2*x0 - c3(1);
                c1(2) = 2*y(k,j) - c3(2);
            case Ny
                xt = x(k-1,j+1);
                if (abs(xt - x0) > 270)
                    xt = xt - 360*sign(xt-x0);
                end;
                c1(1) = 0.5*(xt+x0);
                c1(2) = 0.5*(y(k-1,j+1) + y(k,j));
                c3(1) = 2*x0 - c1(1);
                c3(2) = 2*y(k,j) - c1(2);
                c2(1) = x0 - (y(k,j)-c1(2));
                c2(2) = y(k,j) + (x0-c1(1));
                c4(1) = 2*x0 - c2(1);
                c4(2) = 2*y(k,j) - c2(2);
            otherwise
                xt = x(k-1,j+1);
                if (abs(xt - x0) > 270)
                    xt = xt - 360*sign(xt-x0);
                end;
                c1(1) = 0.5*(xt+x0);
                c1(2) = 0.5*(y(k-1,j+1) + y(k,j));
                xt = x(k+1,j+1);
                if (abs(xt - x0) > 270)
                    xt = xt - 360*sign(xt-x0);
                end;
                c2(1) = 0.5*(xt+x0);
                c2(2) = 0.5*(y(k+1,j+1) + y(k,j));
                c3(1) = 2*x0 - c1(1);
                c3(2) = 2*y(k,j) - c1(2);
                c4(1) = 2*x0 - c2(1);
                c4(2) = 2*y(k,j) - c2(2);
        end;
        
    else
        
        if ( j == Nx )
            
            %@@@ Right edge
            
            switch k
                case 1
                    xt = x(k+1,j-1);
                    if (abs(xt - x0) > 270)
                        xt = xt - 360*sign(xt-x0);
                    end;
                    c3(1) = 0.5*(xt+x0); 
                    c3(2) = 0.5*(y(k+1,j-1) + y(k,j));
                    c2(1) = x0 - (c3(2)-y(k,j));
                    c2(2) = y(k,j) + (c3(1)-x0);
                    c1(1) = 2*x0 - c3(1);
                    c1(2) = 2*y(k,j) - c3(2);
                    c4(1) = 2*x0 - c2(1);
                    c4(2) = 2*y(k,j) - c2(2);
                case Ny
                    xt = x(k-1,j-1);
                    if (abs(xt - x0) > 270)
                        xt = xt - 360*sign(xt-x0);
                    end;
                    c4(1) = 0.5*(xt+x0);
                    c4(2) = 0.5*(y(k-1,j-1) + y(k,j));
                    c3(1) = x0 - (c4(2)-y(k,j));
                    c3(2) = y(k,j) + (c4(1)-x0);
                    c1(1) = 2*x0 - c3(1);
                    c1(2) = 2*y(k,j) - c3(2);
                    c2(1) = 2*x0 - c4(1);
                    c2(2) = 2*y(k,j) - c4(2);
                otherwise
                    xt = x(k+1,j-1);
                    if (abs(xt - x0) > 270)
                        xt = xt - 360*sign(xt-x0);
                    end;
                    c3(1) = 0.5*(xt+x0);
                    c3(2) = 0.5*(y(k+1,j-1) + y(k,j));
                    xt = x(k-1,j-1);
                    if (abs(xt - x0) > 270)
                        xt = xt - 360*sign(xt-x0);
                    end;
                    c4(1) = 0.5*(xt+x0);
                    c4(2) = 0.5*(y(k-1,j-1) + y(k,j));
                    c1(1) = 2*x0 - c3(1);
                    c1(2) = 2*y(k,j) - c3(2);
                    c2(1) = 2*x0 - c4(1);
                    c2(2) = 2*y(k,j) - c4(2);
            end;
            
        else
            
            if ( k == 1 )
                
                %@@@ Bottom edge
                
                xt = x(k+1,j+1);
                if (abs(xt - x0) > 270)
                    xt = xt - 360*sign(xt-x0);
                end;
                c2(1) = 0.5*(xt+x0); 
                c2(2) = 0.5*(y(k+1,j+1) + y(k,j));
                xt = x(k+1,j-1);
                if (abs(xt - x0) > 270)
                    xt = xt - 360*sign(xt-x0);
                end;
                c3(1) = 0.5*(xt+x0);
                c3(2) = 0.5*(y(k+1,j-1) + y(k,j));
                c4(1) = 2*x0 - c2(1);
                c4(2) = 2*y(k,j) - c2(2);
                c1(1) = 2*x0 - c3(1);
                c1(2) = 2*y(k,j) - c3(2);
                
            else
                
                %@@@ Top edge
                
                xt = x(k-1,j-1);
                if (abs(xt - x0) > 270)
                    xt = xt - 360*sign(xt-x0);
                end;
                c4(1) = 0.5*(xt+x0);
                c4(2) = 0.5*(y(k-1,j-1) + y(k,j));
                xt = x(k-1,j+1);
                if (abs(xt - x0) > 270)
                    xt = xt - 360*sign(xt-x0);
                end;
                c1(1) = 0.5*(xt+x0);
                c1(2) = 0.5*(y(k-1,j+1) + y(k,j));
                c2(1) = 2*x0 - c4(1);
                c2(2) = 2*y(k,j) - c4(2);
                c3(1) = 2*x0 - c1(1);
                c3(2) = 2*y(k,j) - c1(2);
                
            end;
            
        end;
        
    end;
    
end;

dx = c1(1)-c4(1);
%argd = min([1 (cosd(c1(2))*cosd(c4(2))*cosd(dx) + sind(c1(2))*sind(c4(2)))]);
%wdth = acosd(argd);
wdth = sqrt((c1(1)-c4(1))^2+(c1(2)-c4(2))^2);
%dx = c2(1)-c1(1);
%argd = min([1 (cosd(c2(2))*cosd(c1(2))*cosd(dx) + sind(c2(2))*sind(c1(2)))]);
%hgt = acosd(argd);
hgt = sqrt((c2(1)-c1(1))^2+(c2(2)-c1(2))^2);

return;
