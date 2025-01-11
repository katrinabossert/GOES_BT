close all
clear all

%nc_filename='/Users/katrina/Dropbox/GraduateSpacePhysicsInstrumentation/GOES_info/OR_ABI-L1b-RadF-M6C13_G18_s20231801200216_e20231801209535_c20231801209567.nc';
%nc_filename='/Users/katrina/Dropbox/GraduateSpacePhysicsInstrumentation/GOES_info/OR_ABI-L1b-RadF-M6C13_G18_s20231818050216_e20231801859536_c20231801859583.nc';
nc_filename='OR_ABI-L1b-RadF-M6C13_G18_s20231802350215_e20231802359538_c20231802359593.nc';

%nc_filename='/Users/katrina/Dropbox/GraduateSpacePhysicsInstrumentation/GOES_info/OR_ABI-L1b-RadF-M6C13_G18_s20231801850216_e20231801859536_c20231801859583.nc';
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
lon_calc=-137.2*pi/180-atan(sy./(H-sx)); %%%%***** FOR FUTURE PROJECTS, GIVE lambda_o

lat_deg=lat_calc.*180./pi;
lon_deg=lon_calc.*180./pi;

figure(1)
pcolor(Rad_sub)
shading flat
colorbar
set(gcf,'color','w')
set(gca, 'FontSize', 16)

figure(11)
pcolor((lat_deg))
shading flat 
colorbar
set(gcf,'color','w')
set(gca, 'FontSize', 16)
title('Latitude (Deg)')

figure(12)
pcolor((lon_deg))
shading flat 
colorbar 
set(gcf,'color','w')
set(gca, 'FontSize', 16)
title('Longitude (Deg)')

%%%Radiances plotted vs. lat/lon (suggestion to use real values to
%%%eliminate error associated with plotting)
%%%%plotted over disk
figure(2)
worldmap([20,60],[-120,-60])
%worldmap([-90,90],[-230,-50])
set(gcf,'color','w')
set(gca, 'FontSize', 16)
setm(gca,'FontSize', 16)
setm(gca,'MLabelParallel','south')
geoshow((lat_deg),(lon_deg),Rad_sub,'DisplayType','texturemap')
colorbar
hold on
plotm(coastlat, coastlon,'color','k','Linewidth',2)
hold off
% set(gcf,'color','w')
% set(gca, 'FontSize', 16)
% setm(gca,'FontSize', 16)
% setm(gca,'MLabelParallel','south')
title('Radiances')

%%%%%calculating temperatures 
fk2=1.38986E3;
fk1=1.073E4;
bc1=.13445;
bc2=.99955;
Temp_act=(planck_fk2./log((planck_fk1./Rad_sub)+1)-planck_bc1)./planck_bc2;

lambda=10.33e-6;
h=6.62e-34;
c=3e8;
k=1.38e-23;

figure(5)
worldmap([20,60],[-120,-60])
set(gcf,'color','w')
set(gca, 'FontSize', 16)
setm(gca,'FontSize', 16)
setm(gca,'MLabelParallel','south')
%worldmap([30,45],[-95,-75])
%worldmap([38,41],[-90,-87]) %zoomed
%geoshow(lat_deg(500:2500,2500:4200),lon_deg(500:2500,2500:4200),Temp_act(500:2500,2500:4200),'DisplayType','texturemap')
geoshow((lat_deg),(lon_deg),Temp_act,'DisplayType','texturemap')
%geoshow(lat_degn,lon_degn,Temp_act,'DisplayType','texturemap')
colorbar
hold on
plotm(coastlat, coastlon,'color','k','Linewidth',2)
hold off
caxis([200,300])
%caxis([192,200])



figure(5)
pcolor((lon_deg),lat_deg,Temp_act)
shading flat
hold on
geoshow(coastlat,coastlon,'color','k','Linewidth',2)
%plot(coastlon,coastlat,10,'k')
%scatter(coastlon,coastlat,10,'k','filled')
hold off
axis([-120,-60,20,60])
set(gcf,'color','w')
set(gca, 'FontSize', 18,'FontWeight','Bold')
title('GOES BT', 'FontSize',24)
xlabel('Longitude','FontSize',18)
ylabel('Latitude','FontSize',18)
h = colorbar
ylabel(h, 'T''[K]')
%caxis([-100,100])

%%%%find the min
Temp_arr=Temp_act(:);
lon_arr=lon_deg(:);
lat_arr=lat_deg(:);

[min_T,ind_T]=min(Temp_arr);
lon_min=lon_arr(ind_T);
lat_min=lat_arr(ind_T);