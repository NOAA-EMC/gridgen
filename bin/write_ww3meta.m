function [messg,errno] = write_ww3meta(fname,gtype,lon,lat,varargin)

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
%| Write the meta data associated with the grids generated in this        |
%| software. This data needs to be provided as input to ww3_grid.inp when |
%| generating the mod_def files for WAVEWATCH III. Note that the paths for|
%| the actual file locations as well as the file names may be changed from|
%| what is written out in the meta file                                   |
%|                                                                        |
%| [messg,errno] = ...                                                    |
%|        write_ww3meta(fname,gtype,lon,lat,N1,N2,N3(optional),N4,N5,N6)  |
%|                                                                        |
%| INPUT                                                                  |
%|  fname       : Output file name prefix                                 | 
%|  gtype       : Grid Type. Two options                                  |
%|                  'CURV' - For curvilinear grids                        |
%|                  'RECT' - For rectilinear grids                        |
%|  lon,lat     : Longitude array (x) and lattitude array (y) of the grid |
%|                 If gtype is 'rect' these are arrays and if it is       |
%|                 'curv' then they are 2D matrices                       |
%|  N1,N2,N3    : Scaling applied to bottom bathymetry data, obstruction  |
%|                grids and coordinate (x,y) grids respectively. The      |
%|                last number is optional and needed only for curvilinear |
%|                grids.                                                  |
%|  N4,N5,N6    : Optional extensions for labeling depth, mask and obs-   |
%|                truction files (must be equal to actual files).         |
%|                                                                        |
%| OUTPUT                                                                 |
%|  messg       : Error message. Is blank if no error occurs              |
%|  errno       : Error number. Is zero for succesful write               |
% -------------------------------------------------------------------------

fid = fopen([fname,'.meta'],'w');

[messg,errno] = ferror(fid);

if (errno ~= 0)
    fprintf(1,'!!ERROR!!: %s \n',messg);
    fclose(fid);
    return;
end;

str1 = '$ Define grid -------------------------------------- $';
str2 = '$ Five records containing :';
str3 = '$  1 Type of grid, coordinate system and type of closure: GSTRG, FLAGLL,';
str4 = '$    CSTRG. Grid closure can only be applied in spherical coordinates.';
str5 = '$      GSTRG  : String indicating type of grid :';
str6 = '$               ''RECT''  : rectilinear';
str7 = '$               ''CURV''  : curvilinear';
str8 = '$      FLAGLL : Flag to indicate coordinate system :';
str9 = '$               T  : Spherical (lon/lat in degrees)';
str10= '$               F  : Cartesian (meters)';
str11= '$      CSTRG  : String indicating the type of grid index space closure :';
str12= '$               ''NONE''  : No closure is applied';
str13= '$               ''SMPL''  : Simple grid closure : Grid is periodic in the';
str14= '$                         : i-index and wraps at i=NX+1. In other words,';
str15= '$                         : (NX+1,J) => (1,J). A grid with simple closure';
str16= '$                         : may be rectilinear or curvilinear.';
str17= '$               ''TRPL''  : Tripole grid closure : Grid is periodic in the';
str18= '$                         : i-index and wraps at i=NX+1 and has closure at';
str19= '$                         : j=NY+1. In other words, (NX+1,J<=NY) => (1,J)';
str20= '$                         : and (I,NY+1) => (MOD(NX-I+1,NX)+1,NY). Tripole';
str21= '$                         : grid closure requires that NX be even. A grid';
str22= '$                         : with tripole closure must be curvilinear.';
str23= '$  2 NX, NY. As the outer grid lines are always defined as land';
str24= '$    points, the minimum size is 3x3.';

fprintf(fid,'%s\n%s\n%s\n%s\n%s\n%s\n',str1,str2,str3,str4,str5,str6);
fprintf(fid,'%s\n%s\n%s\n%s\n%s\n%s\n',str7,str8,str9,str10,str11,str12);
fprintf(fid,'%s\n%s\n%s\n%s\n%s\n%s\n',str13,str14,str15,str16,str17,str18);
fprintf(fid,'%s\n%s\n%s\n%s\n%s\n%s\n',str19,str20,str21,str22,str23,str24);

switch gtype
   case 'RECT'      
       strs1 = '$  3 Grid increments SX, SY (degr.or m) and scaling (division) factor.';
       strs2 = '$    If NX*SX = 360., latitudinal closure is applied.';
       strs3 = '$  4 Coordinates of (1,1) (degr.) and scaling (division) factor.';
       
       fprintf(fid,'%s\n%s\n%s\n',strs1,strs2,strs3);
       
    case 'CURV'
       strs1 = '$  3 Unit number of file with x-coordinate.';
       strs2 = '$    Scale factor and add offset: x <= scale_fac * x_read + add_offset.';
       strs3 = '$    IDLA, IDFM, format for formatted read, FROM and filename.';
       strs4 = '$  4 Unit number of file with y-coordinate.';
       strs5 = '$    Scale factor and add offset: y <= scale_fac * y_read + add_offset.';
       strs6 = '$    IDLA, IDFM, format for formatted read, FROM and filename.';
       
       fprintf(fid,'%s\n%s\n%s\n%s\n%s\n%s\n',strs1,strs2,strs3,strs4,...
           strs5,strs6);
    otherwise
       fprintf(1,'!!ERROR!!: Unrecognized Grid Type\n');
       fclose(fid);
       return;
end;

str25= '$  5 Limiting bottom depth (m) to discriminate between land and sea';
str26= '$    points, minimum water depth (m) as allowed in model, unit number';
str27= '$    of file with bottom depths, scale factor for bottom depths (mult.),';
str28= '$    IDLA, IDFM, format for formatted read, FROM and filename.';
str29= '$      IDLA : Layout indicator :';
str30= '$                  1   : Read line-by-line bottom to top.';
str31= '$                  2   : Like 1, single read statement.';
str32= '$                  3   : Read line-by-line top to bottom.';
str33= '$                  4   : Like 3, single read statement.';
str34= '$      IDFM : format indicator :';
str35= '$                  1   : Free format.';
str36= '$                  2   : Fixed format with above format descriptor.';
str37= '$                  3   : Unformatted.';
str38= '$      FROM : file type parameter';
str39= '$             ''UNIT'' : open file by unit number only.';
str40= '$             ''NAME'' : open file by name and assign to unit.';
str41= '$  If the Unit Numbers in above files is 10 then data is read from this file';
str42= '$';

fprintf(fid,'%s\n%s\n%s\n%s\n%s\n%s\n',str25,str26,str27,str28,str29,str30);
fprintf(fid,'%s\n%s\n%s\n%s\n%s\n%s\n',str31,str32,str33,str34,str35,str36);
fprintf(fid,'%s\n%s\n%s\n%s\n%s\n%s\n',str37,str38,str39,str40,str41,str42);
       
fprintf(fid,'   ''%s''  %s %s\n',gtype,'FLAGLL','CSTRNG');

switch gtype
    case 'RECT'
        N1 = varargin{1};
        N2 = varargin{2};
        [Ny,Nx] = size(lon);
        fprintf(fid,'%d \t %d \n',Nx,Ny);
        fprintf(fid,'%5.2f \t %5.2f \t %5.2f \n',(lon(1,2)-lon(1,1))*60,...
               (lat(2,1)-lat(1,1))*60,60);
        fprintf(fid,'%8.4f \t %8.4f \t %5.2f\n',lon(1,1),lat(1,1),1);
    case 'CURV'
        N1 = varargin{1};
        N2 = varargin{2};
        N3 = varargin{3};
        [Ny,Nx] = size(lon);
        fprintf(fid,'%d \t %d \n',Nx,Ny);
        fprintf(fid,'%d  %f  %5.2f  %d  %d %s  %s  %s  \n',20,N3,0.0,1,...
               1,'''(....)''','NAME',['''',fname,'.lon','''']);
        fprintf(fid,'%d  %f  %5.2f  %d  %d %s  %s  %s  \n',30,N3,0.0,1,...
               1,'''(....)''','NAME',['''',fname,'.lat','''']);
    otherwise
        fprintf(1,'!!ERROR!!: Unrecognized Grid Type\n');
        fclose(fid);
        return;
end;

if nargin <= 7
 ext1='.depth_ascii';
 ext2='.obstr_lev1';
 ext3='.mask';
else
 ext1=varargin{4};
 ext2=varargin{5};
 ext3=varargin{6};
end

fprintf(fid,'$ Bottom Bathymetry \n');
fprintf(fid,'%5.2f  %5.2f  %d  %f  %d  %d %s  %s  %s \n',-0.1,2.5,40,N1,1,1,...
    '''(....)''','NAME',['''',fname,ext1,'''']);
fprintf(fid,'$ Sub-grid information \n');
fprintf(fid,'%d  %f  %d  %d  %s  %s  %s  \n',50,N2,1,1,'''(....)''','NAME',...
    ['''',fname,ext2,'''']);
fprintf(fid,'$ Mask Information \n');
fprintf(fid,'%d  %d  %d  %s  %s  %s  \n',60,1,1,'''(....)''','NAME',...
    ['''',fname,ext3,'''']);

fclose(fid);

return;
