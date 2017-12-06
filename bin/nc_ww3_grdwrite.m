function nc_ww3_grdwrite(x,y,z,file,varargin);
%GRDWRITE2  Write a GMT grid file
%
% Uses built-in NetCDF capability (MATLAB R2008b or later) to 
% write a COARDS-compliant netCDF grid file
% Duplicates (some) functionality of the program grdwrite (which requires
% compilation as a mexfile-based function on each architecture) using
% Matlab 2008b (and later) built-in NetCDF functionality
% instead of GMT libraries.
%
% GRDWRITE2(X,Y,Z,'filename') will create a grid file containing the
% data in the matrix Z.  X and Y should be either vectors with
% dimensions that match the size of Z or two-component vectors
% containing the max and min values for each.
%
% See also GRDREAD2, GRDINFO2

% For more information on GMT grid file formats, see:
% http://www.soest.hawaii.edu/gmt/gmt/doc/gmt/html/GMT_Docs/node70.html
% Details on Matlab's native netCDF capabilities are at:
% http://www.mathworks.com/access/helpdesk/help/techdoc/ref/netcdf.html

% GMT (Generic Mapping Tools, <http://gmt.soest.hawaii.edu>)
% was developed by Paul Wessel and Walter H. F. Smith

% Kelsey Jordahl
% Marymount Manhattan College
% http://marymount.mmm.edu/faculty/kjordahl/software.html

% Time-stamp: <Tue Jul 19 16:28:24 EDT 2011>

% Version 1.1.2, 19-Jul-2011
% Available at MATLAB Central
% <http://www.mathworks.com/matlabcentral/fileexchange/26290-grdwrite2>

fillval=NaN;

write_mask=0;
write_obst=0;

if nargin < 4,
  help(mfilename);
  return,
elseif nargin == 5
  fillval = varargin{1};
elseif nargin == 6
  fillval = varargin{1};
  write_mask = 1;
  mask = varargin{2};
elseif nargin == 8
  fillval = varargin{1};
  write_mask = 1;
  mask = varargin{2};
  write_obst = 1;
  sx = varargin{3};
  sy = varargin{4};
end

% check for appropriate Matlab version (>=7.7)
V=regexp(version,'[ \.]','split');
if (str2num(V{1})<7) | (str2num(V{1})==7 & str2num(V{2})<7),
  ver
  error('grdread2: Requires Matlab R2008b or later!');
end

ncid = netcdf.create(file, 'NC_SHARE');
if isempty(ncid),
  return,
end

% set descriptive variables
conv='COARDS/CF-1.0';
title=file;
history='File written by MATLAB function nc_ww3_grdwrite.m';
desc=['Created ' datestr(now)];
vers='4.x';                             % is "x" OK?

% check X and Y
if (~isvector(x) | ~isvector(y)),
  error('X and Y must be vectors!');
end
if (length(x) ~= size(z,2)),    % number of columns don't match size of x
  minx=min(x); maxx=max(x);
  dx=(maxx-minx)/(size(z,2)-1);
  x=minx:dx:maxx;                       % write as a vector
end
if (length(y) ~= size(z,1)),    % number of rows don't match size of y
  miny=min(y); maxy=max(y);
  dy=(maxy-miny)/(size(z,1)-1);
  y=miny:dy:maxy;                       % write as a vector
end

% Convert all values to single and scale
scalfz = 1;
z = single(scalfz * z);

if nargin >= 6
  scalfm = 1;
  mask = single(scalfm * mask);
end

if nargin == 8
  scalfs = 1;
  sx = single(scalfs * sx);
  sy = single(scalfs * sy);
end

% match Matlab class to NetCDF data type
  nctype_ll='NC_FLOAT';

  nctype='NC_FLOAT';
  nanfill=single(fillval);
  disp(['Fill value set to ' num2str(nanfill)])

% global
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Conventions',conv);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'title',title);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'history',history);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'description',desc);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'GMT_version',vers);
% X
dimid = netcdf.defDim(ncid,'lon',length(x));
varid = netcdf.defVar(ncid,'lon',nctype_ll,dimid);
netcdf.putAtt(ncid,varid,'long_name','longitude');
netcdf.putAtt(ncid,varid,'actual_range',[min(x) max(x)]);
netcdf.endDef(ncid);
netcdf.putVar(ncid,varid,x);
% Y
netcdf.reDef(ncid);
dimid = netcdf.defDim(ncid,'lat',length(y));
varid = netcdf.defVar(ncid,'lat',nctype_ll,dimid);
netcdf.putAtt(ncid,varid,'long_name','latitude');
netcdf.putAtt(ncid,varid,'actual_range',[min(y) max(y)]);
netcdf.endDef(ncid);
netcdf.putVar(ncid,varid,y);
% Z
addoff=0.;
netcdf.reDef(ncid);
varid = netcdf.defVar(ncid,'z',nctype,[0 1]);
netcdf.putAtt(ncid,varid,'long_name','Depth');
netcdf.putAtt(ncid,varid,'_FillValue',nanfill);
netcdf.putAtt(ncid,varid,'scale_factor',scalfz);
netcdf.putAtt(ncid,varid,'add_offset',addoff);
netcdf.putAtt(ncid,varid,'actual_range',[min(z(:)) max(z(:))]);
netcdf.endDef(ncid);
netcdf.putVar(ncid,varid,z');
%
if write_mask == 1
% Mask
netcdf.reDef(ncid);
varid = netcdf.defVar(ncid,'mask',nctype,[0 1]);
netcdf.putAtt(ncid,varid,'long_name','mask');
netcdf.putAtt(ncid,varid,'_FillValue',nanfill);
netcdf.putAtt(ncid,varid,'scale_factor',scalfm);
netcdf.putAtt(ncid,varid,'add_offset',addoff);
netcdf.putAtt(ncid,varid,'actual_range',[min(mask(:)) max(mask(:))]);
netcdf.endDef(ncid);
netcdf.putVar(ncid,varid,mask');
end
%
if write_obst == 1
% Sx
netcdf.reDef(ncid);
varid = netcdf.defVar(ncid,'sx',nctype,[0 1]);
netcdf.putAtt(ncid,varid,'long_name','Subgrd obstr x');
netcdf.putAtt(ncid,varid,'_FillValue',nanfill);
netcdf.putAtt(ncid,varid,'scale_factor',scalfs);
netcdf.putAtt(ncid,varid,'add_offset',addoff);
netcdf.putAtt(ncid,varid,'actual_range',[min(sx(:)) max(sx(:))]);
netcdf.endDef(ncid);
netcdf.putVar(ncid,varid,sx');
% Sy
netcdf.reDef(ncid);
varid = netcdf.defVar(ncid,'sy',nctype,[0 1]);
netcdf.putAtt(ncid,varid,'long_name','Subgrd obstr y');
netcdf.putAtt(ncid,varid,'_FillValue',nanfill);
netcdf.putAtt(ncid,varid,'scale_factor',scalfs);
netcdf.putAtt(ncid,varid,'add_offset',addoff);
netcdf.putAtt(ncid,varid,'actual_range',[min(sy(:)) max(sy(:))]);
netcdf.endDef(ncid);
netcdf.putVar(ncid,varid,sy');
end
%
% close file
netcdf.close(ncid);

