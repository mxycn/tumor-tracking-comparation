function [numbers,errors] = simulation_patient2(conditions,path_in)

% Compared to simulation, this function adopts PCA method.
% This function performs data input, reconstruction, and prediction, 
% error calculation and reporting

%% the errors include the following:
% path_in: the folder of the data
% rms,std_dev,devp, count95_3
tic()
%% conditions include the following:
% one sample point are taken ever incrt data points, whose 
% period is approximately 0.038545s. incrt(= 4)
incrt = conditions(14);

%% data initializing ...
% read data from the data file
global patient_num
ModelerDataDir = path_in; % the folder of data
cd(ModelerDataDir); 
dbfolders = dir('DB*');
ndir = length(dbfolders); % the number of folders
npts=zeros(ndir,1); % the number of points in each folder
% where npts is the number of points which equals the total original data points
% divided by incrt
% define an array to accept the data
o3d=zeros(ndir,300000,4); % the array storing the original point data
ab_dev2=zeros(ndir,1);  % the sum of the square of 2D error
std_dev2=zeros(ndir,1); 
nfail2x=zeros(ndir,1);
pc95_2d=zeros(ndir,1); 
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
    [k_segment,std_k2,pro2d_k2,nfail2,pct95_2d,pct99_2d] = ...
        simulation_fraction2(npts(k),ko3d,conditions);
%  
     n_seg=n_seg + k_segment; % the number of trajectories increases by n_segment.
   
    if k_segment>=1
        std_dev2(k)= std_k2;
        ab_dev2(k)= pro2d_k2; 
        nfail2x(k) = nfail2;
        pc95_2d(k) = pct95_2d;
        pc99_2d(k) = pct99_2d;
        ndir_count = ndir_count + 1;
    else
        pc95_2d(k)=0;
        pc99_2d(k)=0;
        std_dev2(k)=0;
        ab_dev2(k)=0;
        nfail2x(k) = 0;

    end
end % of k...
%% reporting results
std_dev=mean(sqrt(std_dev2(1:ndir_count))); % total average standard deviation of the patient.
devp=mean(sqrt(ab_dev2(1:ndir_count)));     % total average 2D error
p95_2d=mean(pc95_2d(1:ndir_count));
p99_2d=mean(pc99_2d(1:ndir_count));
numberfail2=mean(nfail2x(1:ndir_count));

errors=[std_dev,devp,p95_2d,p99_2d,numberfail2];
numbers = [ndir,nactual,n_seg];
xx=toc();
fidd=fopen('D:\Program Files\MATLAB\errorstatistics150929\error_analysis\stat_result2d_individual','a'); 
 
fprintf(fidd,' 2D beams-eye-view error: %4.2f patient %u time %4.2f \n',devp,patient_num,xx/n_seg); 
fclose(fidd);
patient_num=patient_num+1;
% indicate the analysis of one patient is finished
fprintf('the number of trajectories = %u \n', n_seg); 
fprintf('the number of fractions = %u \n', ndir_count); 
