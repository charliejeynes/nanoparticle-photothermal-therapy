%%
% created by Charlie Jeynes 2019

clc
clear
close all


%% open the .nc file 
% source = '/Users/charliejeynes/Projects/arc/output/mcrt/lm_abs_dens.nc'
% source = '/Users/charliejeynes/Projects/arc/output/mcrt/mat_map_ptfe.nc'
% source = '/Users/charliejeynes/Projects/arc/output/mcrt/interfaces.nc'
% source = '/Users/charliejeynes/Projects/absorption_density.nc'
source = '/Users/charliejeynes/Projects/dia/dia/output/mcrt/absorption_density.nc'
% source = '/Users/charliejeynes/Projects/dia/dia/output/mcrt/hits.nc'
ncid = netcdf.open(source); 
data = netcdf.getVar(ncid,0); 
netcdf.close(ncid); 

%%

xslice = []
yslice = []
zslice = [32]

hss = slice(data,xslice,yslice,zslice); 
xs = get(hss,'XData');
ys = get(hss,'YData');
zs = get(hss,'ZData');
cs = get(hss,'CData');

cs_log_rotate = rot90(log10(cs), 2); 

figure, 
imagesc(cs_log_rotate(:, :));
caxis([0,8.1])% create and label the colorbar
% imagesc(cs(:, :));
% caxis([0,1e8])% create and label the colorbar
cmap = jet();
colormap(cmap);
cb=colorbar
cb.Label.String = 'absorbance density (W/m^3) (log10)';
%cb.Label.String = 'hits';
xlabel('mm');
xticklabels({'0','5','10','15','20','25','30','35', '40'});
xticks(0:8:64)
yticklabels({'0','5','10','15','20','25','30','35', '40'});
yticks(0:8:64)
ylabel('mm');
title('tumour')


% xticklabels({'0','1','2','3','4','5','6','7', '8', '9', '10'});
% xticklabels({'0','2.5','5','7.5','10','12.5','15','17.5', '20'});
% xticks(0:8:64)
% yticklabels({'0','2.5','5','7.5','10','12.5','15','17.5', '20'});
% yticks(0:8:64)
% what is the scale of the image 
scale_mm = (40 / 64) * 8
%%
64/8
%% get a 3D image
xslice = []
yslice = [32]
zslice = [32]
slice(log10(rot90(data)),xslice,yslice,zslice) 
% caxis([0,2e6])% create and label the colorbar
cmap = jet();
cb=colorbar
cb.Label.String = 'absorbance density (W/m^3) (log10)';
xlabel('mm');
% xticklabels({'0','1','2','3','4','5','6','7', '8', '9', '10'});
xticklabels({'0','2.5','5','7.5','10','12.5','15','17.5', '20'});
xticks(0:8:64)
yticklabels({'0','2.5','5','7.5','10','12.5','15','17.5', '20'});
yticks(0:8:64)
ylabel('mm');
zticklabels({'0','2.5','5','7.5','10','12.5','15','17.5', '20'});
zticks(0:8:64)
zlabel('mm');
title('tumour')

%% plot a single image 
figure, 
imagesc(data(:, : , 1));
% caxis([0,2])% create and label the colorbar
% cb.Label.String = 'absorb';
colorbar



%% get all files out of folder

% filenameS = dir('/Users/jcgj201/Documents/MATLAB/STFC innovation/vtk files to open/billion photon output'); 
filenameS = dir('/Users/jcgj201/Documents/MATLAB/STFC innovation/vtk files to open/2nd attempt output');  %point to your file

%% call the function on all files
csRotmaster = zeros(201,201,1); 
filenameLst = {}; 
for i= 4:9
    filenameShort = filenameS(i).name; 
    filenameLst{i} = filenameS(i).name
%     filenameLong = ['/Users/jcgj201/Documents/MATLAB/STFC innovation/vtk files to open/billion photon output/' filenameShort]
    filenameLong = ['/Users/jcgj201/Documents/MATLAB/STFC innovation/vtk files to open/2nd attempt output/' filenameShort]; 
    csRot = vtk_analyse(filenameShort, filenameLong); 
    csRotmaster(:,:,i) = csRot; 
end

% %% plot a single image 
% figure, 
% imagesc(csRotmaster(:, : , 5));
% caxis([1e8,1e16])% create and label the colorbar
% cb.Label.String = 'absorb';
% colorbar
% 
% %%
% 
% imtool(csRotmaster(:, : , 5)); 
%     
% %% rescale the array of slices (didnt use this) 
% 
% rescaleCSrot = rescale(csRotmaster); 
% 
% %% min and max of csRotmaster
% 
% max  = csRotmaster(:,:,5)
% 
% csRotmaster(:,:,5)
% numel(csRotmaster(:,:,5))
% 
% %% check it's getting data correctly 201
% 
% 
%     fid=fopen('half_tumour_plus_au.vtk','r');   % open the file; be sure to add the error-checking code
%     l=fgetl(fid);                % read first record
%     while ~strcmp(l,'absorb_dens 1 8120601 double')       % loop until find the magic record
%       l=fgetl(fid);
%     end
%     data=textscan(fid,'','collectoutput',1);  % read, convert cell to double array
%     fid=fclose(fid);             % done with file, go on..
%     
%     data_1 = cell2mat(data); 
%   
%  %% test using filelong
%     
%     fid=fopen(filenameLong,'r');   % open the file; be sure to add the error-checking code
%     l=fgetl(fid);                % read first record
%     while ~strcmp(l,'absorb_dens 1 1030301 double')       % loop until find the magic record
%       l=fgetl(fid);
%     end
%     data=cell2mat(textscan(fid,'','collectoutput',1));  % read, convert cell to double array
%     fid=fclose(fid);             % done with file, go on..
%     
%     dataReshape = reshape(data, 101, 101, 101); 
%     
%     a = dataReshape(50, 50, 50)
%     %% check it's getting data correctly 101
% 
% 
%     fid=fopen('base_copy.txt','r');   % open the file; be sure to add the error-checking code
%     l=fgetl(fid);                % read first record
%     while ~strcmp(l,'absorb_dens 1 1030301 double')       % loop until find the magic record
%       l=fgetl(fid);
%     end
%     data=textscan(fid,'','collectoutput',1);  % read, convert cell to double array
%     fid=fclose(fid);             % done with file, go on..
%     
%     data_1 = cell2mat(data); 
%     
%%
titles = {}
%% plot out csRotmaster (which is a 2d dataset)

csRotmasterLog = log10(csRotmaster);

for i=1:6
%    for k = 3:8
       k = i+3; 
       
       hold on
%        figure     
       subplot(2,3,i); % (mxnxp I thikn p is position)
       crop = imcrop(csRotmasterLog(:, : , k), [50 70 100 60]); 
       imagesc(crop); 
       title(filenameLst{k}); 
       cb = colorbar;  
       caxis([3,6])% create and label the colorbar
       cb.Label.String = 'energy density (log10)';
       hold off
%    end
end 

%% plot out RescaleCSRot with log(10)

csRotmasterLog = log10(csRotmaster); 

for i=1:6
%    for k = 3:8
       k = i+2; 
     
       hold on
%        figure,     
       subplot(2,3,i); % (mxnxp I thikn p is position)
       imagesc(csRotmasterLog(:, : , k)); 
       title(filenameLst{k}); 
       cb = colorbar;  
       caxis([12,16])% create and label the colorbar
       cb.Label.String = 'absorbance density';
       hold off
%    end
end 

%% Plot csRotmaster but as a line profile

csRotmasterLog = log10(csRotmaster);
imtool(csRotmasterLog(:, : , 7)); 

norm  =  csRotmasterLog(:,100,7); 
norm_plus_au  = csRotmasterLog(:,100,8); 
x = 1:201; 
figure, 
plot(x, norm, x, norm_plus_au)
legend('norm', 'norm_plus_au')
xlabel('mm')
ylabel('absorbance density (W/m^2)')
% xticklabels({'0','1','2','3','4','5','6','7', '8', '9', '10'});
% xticks(0:10:100)

%% Plot csRotmaster but as a line profile
base  =  csRotmaster(50,:,3); 
gold_e_8_Ua10  = csRotmaster(50,:,7); 
x = 1:201; 
figure, 
plot(x, base, x, gold_e_8_Ua10)
legend('base', 'gold_e_8 or Ua=10')
xlabel('mm')
ylabel('absorbance density (W/m^2)')
xticklabels({'0','1','2','3','4','5','6','7', '8', '9', '10'});
xticks(0:10:100)

%% new info from Freddy on the energy density

% energy density = J.m3.s

% dT = enDensity x time x ( watts / specific_heat x density_tissue x voxel) 
%% back of the envelop calc on dT 

% 200 mW/cm-2 for 37 J/cm-2
% 0.2 W/cm-2 
watts_p_m2 = 0.2 * 10000

totalDose = 0.2 * 185

% dT

dT = 37 / (0.0012 * 3700) % this is 37 Joules / (mass of tumour in kg [cm3] * specific heat of tissue)

%% heat up m-3 of water for 60 seconds with a 630 1W LED
time_irradiated = 1
Q = 1 * time_irradiated % 1W for 60 seconds
specific_heat_water = 4100 % kg/oC J
mass_of_water = 0.01 % 1m3 weighs 1000kg
dT = (Q) / (mass_of_water * specific_heat_water) % this is 37 Joules / (mass of tumour VOXEL in kg [cm3] * specific heat of tissue)

oneW_per_10mm3 = dT / 1e-8 % convert dT in 1m3 to in 10mm3

oneW_per_0point1mm3 = dT / 1e-10 % convert dT in 1m3 to in 0.1mm3 (same as voxel)


% watts_per_m2 = 1e15 * 1e-10 % freddy code value in m3 to 0.1mm3
% seconds_irradiated = 60
% mass_tumour_voxel = 1000 * 1e-10 % 1000 kg per m3 to 0.1mm3


%% watts m-2 in each voxel 

% each voxel is 0.1 mm3 or 1e-10 m2
watts_per_voxel = 1e15 * 1e-10 % freddy code value in m3 to 0.1mm3
seconds_irradiated = 60
mass_tumour_voxel = 1000 * 1e-10 % 1000 kg per m3 to 0.1mm3
specific_heat_tissue = 3700 % kg/oC J

dTperVoxel = (watts_per_voxel) / (mass_tumour_voxel * specific_heat_tissue) % this is 37 Joules / (mass of tumour VOXEL in kg [cm3] * specific heat of tissue)


%% 1W spead over 1m3





% 
% Heat capacity = mass x specific heat x change in temperature
% Q = mc ? T
% Q = heat capacity, J
% m = mass, g
% c = specific heat of object, J/(g-ºC)
% dT = change in temperature, ºC
% Q/mc = ? T
% dT = Q/(c * m) 

voxel = 10/100 % each voxel is 0.1 mm3 or 100 um3
meterConversion  = 0.01 / 1000
% water weighs 1g per cubic cm (i.e. 10mm3)
% mass 1L is 1 Kg so 1ml is 1g
% mass 10mm3 = 1g
% 1mm3 = 0.1g
% 0.1 mm3 = 0.01g

% the unit are in watts per m2 - need to get them to 
c = 3.7 % J/(g-ºC)
m = 0.01 % in g

dT = (1e14 * meterConversion) / (3.7 * 0.01)


%%
figure,     
subplot(2,3,1); % (mxnxp I thikn p is position)
imagesc(csRotmaster(:, : , 3)); 
subplot(2,3,2); % (mxnxp I thikn p is position)
imagesc(csRotmaster(:, : , 4)); 
subplot(2,3,3); % (mxnxp I thikn p is position)
imagesc(csRotmaster(:, : , 5)); 
subplot(2,3,4); % (mxnxp I thikn p is position)
imagesc(csRotmaster(:, : , 6));
subplot(2,3,5); % (mxnxp I thikn p is position)
imagesc(csRotmaster(:, : , 7));
subplot(2,3,6); % (mxnxp I thikn p is position)
imagesc(csRotmaster(:, : , 8));

% imshow(ind2rgb(X, map));
%%
function csRot = vtk_analyse(filenameShort, filenameLong)
    %% code from help

    fid=fopen(filenameLong,'r');   % open the file; be sure to add the error-checking code
    l=fgetl(fid);                % read first record
    while ~strcmp(l,'energy_dens 1 8120601 double')       % loop until find the magic record
      l=fgetl(fid);
    end
    data=cell2mat(textscan(fid,'','collectoutput',1));  % read, convert cell to double array
    fid=fclose(fid);             % done with file, go on..

    %% reshape the data to one long line, make 0 NaN

    absorb = reshape(data, 8120601, 1); 

    %% get x,y,z co-ordinates
    
%     x = [-0.00500000000000000,-0.00490099000000000,-0.00480198000000000,-0.00470297000000000,-0.00460396000000000,-0.00450495000000000,-0.00440594000000000,-0.00430693000000000,-0.00420792000000000,-0.00410891000000000,-0.00400990000000000,-0.00391089000000000,-0.00381188000000000,-0.00371287000000000,-0.00361386000000000,-0.00351485000000000,-0.00341584000000000,-0.00331683000000000,-0.00321782000000000,-0.00311881000000000,-0.00301980000000000,-0.00292079000000000,-0.00282178000000000,-0.00272277000000000,-0.00262376000000000,-0.00252475000000000,-0.00242574000000000,-0.00232673000000000,-0.00222772000000000,-0.00212871000000000,-0.00202970000000000,-0.00193069000000000,-0.00183168000000000,-0.00173267000000000,-0.00163366000000000,-0.00153465000000000,-0.00143564000000000,-0.00133663000000000,-0.00123762000000000,-0.00113861000000000,-0.00103960000000000,-0.000940594000000000,-0.000841584000000000,-0.000742574000000000,-0.000643564000000000,-0.000544554000000000,-0.000445545000000000,-0.000346535000000000,-0.000247525000000000,-0.000148515000000000,-4.95050000000000e-05,4.95050000000000e-05,0.000148515000000000,0.000247525000000000,0.000346535000000000,0.000445545000000000,0.000544554000000000,0.000643564000000000,0.000742574000000000,0.000841584000000000,0.000940594000000000,0.00103960000000000,0.00113861000000000,0.00123762000000000,0.00133663000000000,0.00143564000000000,0.00153465000000000,0.00163366000000000,0.00173267000000000,0.00183168000000000,0.00193069000000000,0.00202970000000000,0.00212871000000000,0.00222772000000000,0.00232673000000000,0.00242574000000000,0.00252475000000000,0.00262376000000000,0.00272277000000000,0.00282178000000000,0.00292079000000000,0.00301980000000000,0.00311881000000000,0.00321782000000000,0.00331683000000000,0.00341584000000000,0.00351485000000000,0.00361386000000000,0.00371287000000000,0.00381188000000000,0.00391089000000000,0.00400990000000000,0.00410891000000000,0.00420792000000000,0.00430693000000000,0.00440594000000000,0.00450495000000000,0.00460396000000000,0.00470297000000000,0.00480198000000000,0.00490099000000000,0.00500000000000000]; 
    x = [-0.005  -0.00495025  -0.0049005  -0.00485075  -0.004801  -0.00475124  -0.00470149  -0.00465174  -0.00460199  -0.00455224  -0.00450249  -0.00445274  -0.00440299  -0.00435323  -0.00430348  -0.00425373  -0.00420398  -0.00415423  -0.00410448  -0.00405473  -0.00400498  -0.00395522  -0.00390547  -0.00385572  -0.00380597  -0.00375622  -0.00370647  -0.00365672  -0.00360697  -0.00355721  -0.00350746  -0.00345771  -0.00340796  -0.00335821  -0.00330846  -0.00325871  -0.00320896  -0.0031592  -0.00310945  -0.0030597  -0.00300995  -0.0029602  -0.00291045  -0.0028607  -0.00281095  -0.00276119  -0.00271144  -0.00266169  -0.00261194  -0.00256219  -0.00251244  -0.00246269  -0.00241294  -0.00236318  -0.00231343  -0.00226368  -0.00221393  -0.00216418  -0.00211443  -0.00206468  -0.00201493  -0.00196517  -0.00191542  -0.00186567  -0.00181592  -0.00176617  -0.00171642  -0.00166667  -0.00161692  -0.00156716  -0.00151741  -0.00146766  -0.00141791  -0.00136816  -0.00131841  -0.00126866  -0.00121891  -0.00116915  -0.0011194  -0.00106965  -0.0010199  -0.000970149  -0.000920398  -0.000870647  -0.000820896  -0.000771144  -0.000721393  -0.000671642  -0.000621891  -0.000572139  -0.000522388  -0.000472637  -0.000422886  -0.000373134  -0.000323383  -0.000273632  -0.000223881  -0.000174129  -0.000124378  -7.46269e-05  -2.48756e-05  2.48756e-05  7.46269e-05  0.000124378  0.000174129  0.000223881  0.000273632  0.000323383  0.000373134  0.000422886  0.000472637  0.000522388  0.000572139  0.000621891  0.000671642  0.000721393  0.000771144  0.000820896  0.000870647  0.000920398  0.000970149  0.0010199  0.00106965  0.0011194  0.00116915  0.00121891  0.00126866  0.00131841  0.00136816  0.00141791  0.00146766  0.00151741  0.00156716  0.00161692  0.00166667  0.00171642  0.00176617  0.00181592  0.00186567  0.00191542  0.00196517  0.00201493  0.00206468  0.00211443  0.00216418  0.00221393  0.00226368  0.00231343  0.00236318  0.00241294  0.00246269  0.00251244  0.00256219  0.00261194  0.00266169  0.00271144  0.00276119  0.00281095  0.0028607  0.00291045  0.0029602  0.00300995  0.0030597  0.00310945  0.0031592  0.00320896  0.00325871  0.00330846  0.00335821  0.00340796  0.00345771  0.00350746  0.00355721  0.00360697  0.00365672  0.00370647  0.00375622  0.00380597  0.00385572  0.00390547  0.00395522  0.00400498  0.00405473  0.00410448  0.00415423  0.00420398  0.00425373  0.00430348  0.00435323  0.00440299  0.00445274  0.00450249  0.00455224  0.00460199  0.00465174  0.00470149  0.00475124  0.004801  0.00485075  0.0049005  0.00495025  0.005];  

    x = x(1:201); 
    y = x(1:201); 
    z = x(1:201);


    [X,Y,Z] = meshgrid(x,y,z);
    
    %% reshape x,y,z

    xR = reshape(X,8120601, 1); 
    yR = reshape(Y,8120601, 1); 
    zR = reshape(Z,8120601, 1); 
    
    %% get 2D slice 
    
    absorbCube = reshape(absorb, 201, 201, 201); 
    
    xslice = [];                               % define the cross sections to view
    yslice = [];
    zslice = 0;

    hss = slice(X,Y,Z,absorbCube,xslice, yslice, zslice);  
    xs = get(hss,'XData');
    ys = get(hss,'YData');
    zs = get(hss,'ZData');
    cs = get(hss,'CData');
    
    
%     %% Plot the 2D image
    csRot = rot90(cs); 
%     figure, 
%     imagesc(csRot)
%     cb = colorbar;  
%     % caxis([1e13,1e14])% create and label the colorbar
%     cb.Label.String = 'absorb';
% %     filenameCut = filename(80:end); 
%     title(filenameShort); 
%     
%     filenameSave = [filenameCut '.png']
%     saveas(gcf, 'filenameSave.png')

end