function [numbers,errors] = simulation_patient_error(conditions,path_in)

% function [numbers,errors,datax] = simulation_patient(conditions,path_in)

% Compared to simulation, this function adopts PCA method.
% This function performs data input, reconstruction, and prediction, 
% error calculation and reporting

%% the errors include the following:
% path_in: the folder of the data
% rms,std_dev,devp, count95_3

%% conditions include the following:
% one sample point are taken ever incrt data points, whose 
% period is approximately 0.038545s. incrt(= 4)
incrt = conditions(14);

%% data initializing ...
% read data from the data file

ModelerDataDir = path_in; % the folder of data
cd(ModelerDataDir); 
dbfolders = dir('DB*');
ndir = length(dbfolders); % the number of folders
npts=zeros(ndir,1); % the number of points in each folder
% where npts is the number of points which equals the total original data points
% divided by incrt
% define an array to accept the data
o3d=zeros(ndir,300000,4); % the array storing the original point data
aberror2=zeros(ndir,1); % the sum of the square of 3D error
ab_dev2=zeros(ndir,1);  % the sum of the square of 2D error
std_dev2=zeros(ndir,1); 
nkVx=zeros(ndir,1);
nfailx=zeros(ndir,1);
avedisx=zeros(ndir,1);
nrepx = zeros(ndir,1);

pc95_3d=zeros(ndir,1); 
pc95_2d=zeros(ndir,1); 
pc99_3d=zeros(ndir,1); 
pc99_2d=zeros(ndir,1); 
% simple statistics on the data
[o3d,npts]=read_data(ndir,incrt,path_in);
%% initializing results statistics  ...
% the number of points used
nactual=0;  % the number of points available
n_seg=0;    % the number of trajectories
% count95_3=0; % the number of time points with error greater than 3
% count95_4=0;
% the sum of the square of standard deviation

% timex=[0 0];
n_seg=0;
ndir_count = 0;
%% simulate with the data of a patient.
% k is the number of folders
% ndir is the total number of fractions for that patient.

for k=1:ndir
    ko3d=zeros(npts(k),4);
    nactual=nactual+npts(k); % the number of data points read.
             % the number of trajectories
   
 % localize the data ... 
    for ii=1:npts(k)
        ko3d(ii,1)=o3d(k,ii,1);
        ko3d(ii,2)=o3d(k,ii,2);
        ko3d(ii,3)=o3d(k,ii,3);
        ko3d(ii,4)=o3d(k,ii,4);
    end             
 % simulate with the whole data in a subfolder(corresponding to a fraction)    
     [k_segment,std_k2,pre3d_k2,pro2d_k2,threedx,twodx,parametersf] = ...
         simulation_fraction_error(npts(k),ko3d,conditions);
%  
     n_seg=n_seg + k_segment; % the number of trajectories increases by n_segment.
     if n_seg == k_segment
%         r3derr = r3derrx;
%         p3derr = p3derrx;
%         p2derr = p2derrx;
        threed = threedx;
        twod = twodx;
        parametersp = parametersf;
    else
%         r3derr = [r3derr;r3derrx];
%         p3derr = [p3derr;p3derrx];
%         p2derr = [p2derr;p2derrx];
        threed = [threed;threedx];
        twod = [twod;twodx];
        parametersp = [parametersp; parametersf];
     end    
     
    if k_segment>=1
        std_dev2(k)= std_k2;
        aberror2(k)= pre3d_k2;
        ab_dev2(k)= pro2d_k2; 
        ndir_count = ndir_count + 1;
    else
        std_dev2(k)=0;
        aberror2(k)=0;
        ab_dev2(k)=0;        
    end
end % of k...
%% reporting results
std_dev=mean(sqrt(std_dev2(1:ndir_count))); % total average standard deviation of the patient.
rms=mean(sqrt(aberror2(1:ndir_count)));     % total average rms 3D error
devp=mean(sqrt(ab_dev2(1:ndir_count)));     % total average 2D error
numbers = [ndir,nactual,n_seg];
% datax = [parametersp,r3derr,p3derr,p2derr];
datax = [parametersp,threed,twod];
errors=[std_dev,rms,devp];
%datax1: velocity
%datax2: cos2
%datax3: preerr
%datax4: prediction quality
%datax5: cost
%datax6: 3d error

fiddx=fopen(strcat(path_in,'\scalar_error'),'w+');
for jj = 1:20:size(datax,1)
    fprintf(fiddx,'%d %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f\n',...
        jj, datax(jj,1),datax(jj,2),datax(jj,3),datax(jj,4),datax(jj,5),datax(jj,6));
end
% for jj = 1:2:size(datax,1)
%     fprintf(fiddx,'%5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f\n',...
%         datax(jj,1),datax(jj,2),datax(jj,3),datax(jj,4),datax(jj,5),datax(jj,6),datax(jj,7),datax(jj,8),datax(jj,9),datax(jj,10),datax(jj,11));
% end
fclose(fiddx);  
% indicate the analysis of one patient is finished
fprintf('the number of trajectories = %u \n', n_seg); 
fprintf('the number of fractions = %u \n', ndir_count); 
