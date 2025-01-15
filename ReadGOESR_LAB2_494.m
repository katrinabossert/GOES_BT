close all
clear all

nc_filename='OR_ABI-L1b-RadF-M6C13_G18_s20231802350215_e20231802359538_c20231802359593.nc';

ncid=netcdf.open(nc_filename,'nowrite'); 


% Get information about the contents of the file.
[numdims, numvars, numglobalatts, unlimdimID] = netcdf.inq(ncid);

disp(' '),disp(' '),disp(' ')
disp('________________________________________________________')
disp('^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~')
disp(['VARIABLES CONTAINED IN THE netCDF FILE: ' nc_filename ])
disp(' ')
for i = 0:numvars-1
    [varname, xtype, dimids, numatts] = netcdf.inqVar(ncid,i);
    disp(['--------------------< ' varname ' >---------------------'])
    flag = 0;
    for j = 0:numatts - 1
        attname1 = netcdf.inqAttName(ncid,i,j);
        attname2 = netcdf.getAtt(ncid,i,attname1);
        disp([attname1 ':  ' num2str(attname2)])
        if strmatch('add_offset',attname1)
            offset = attname2;
        end
        if strmatch('scale_factor',attname1)
            scale = attname2;
            flag = 1;
        end        
    end
    disp(' ')
    
    if flag
        eval([varname '= double(double(netcdf.getVar(ncid,i))*scale + offset);'])
    else
        eval([varname '= double(netcdf.getVar(ncid,i));'])
        
    end
end
disp('^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~')
disp('________________________________________________________')
disp(' '),disp(' ')

%%%%the commands below load latitude and longitude grid points that can be used to plot
%%%%the outlines of land masses 
load coastlines
[latcells, loncells] = polysplit(coastlat, coastlon);
numel(latcells);

%%%conversion of grid points to latitude and longitude 
arr_x=ones(length(x),1);
x_matf=x*arr_x'; % these are the same as x and y in the GOES user matrix 
y_matf=arr_x*y';
r_eq=6378137; %km radius at equator
r_pol=6356752; %km radius at pole
H=42164160; %GOES distance from earth center
lambda_0=-137; %degrees, the Satellite West longitude East for GOESR18 West

%%%Note for this lab, we will work with a subset of data. This has been
%%%predetermined by your instructor. 
st_x=3800;
stp_x=4650;
st_y=850;
stp_y=1210;
x_mat=x_matf(st_x:stp_x,st_y:stp_y);
y_mat=y_matf(st_x:stp_x,st_y:stp_y);
Rad_sub=Rad(st_x:stp_x,st_y:stp_y);

%%%follow outlined procedure for calculating latitude and longitude in the
%%%GOES user manual pages 19-23, section 5.1.2.8
aa=(sin(x_mat)).^2+((cos(x_mat)).^2).*(((cos(y_mat)).^2)+((r_eq^2)/(r_pol^2)).*sin(y_mat).^2); 
bb=-2.*H.*cos(x_mat).*cos(y_mat);
cc=H.^2-r_eq.^2;

rs=(-bb-sqrt((bb.^2-4.*aa.*cc)))./(2.*aa);
sx=rs.*cos(x_mat).*cos(y_mat);
sy=-rs.*sin(x_mat);
sz=rs.*cos(x_mat).*sin(y_mat);

lat_calc=atan(((r_eq^2)/(r_pol^2)).*sz./(sqrt((H-sx).^2+sy.^2)));
lon_calc=-137.2*pi/180-atan(sy./(H-sx)); 

lat_deg=lat_calc.*180./pi;
lon_deg=lon_calc.*180./pi;
%%%note, the above lat_deg and lon_deg are given over the same region as
%%%Rad_sub
%%%If you would like to calculate over the entire field of view, it is
%%%suggested that you take the real value of lat_deg and lon_deg for
%%%plotting, as the edges/off disk may produce imaginary results

%Problem 3
figure(1)
