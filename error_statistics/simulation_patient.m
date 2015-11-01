function [numbers,errors] = simulation_patient(conditions,path_in)

% function [numbers,errors,datax] = simulation_patient(conditions,path_in)

% Compared to simulation, this function adopts PCA method.
% This function performs data input, reconstruction, and prediction, 
% error calculation and reporting
global patient_num
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
rce3d = zeros(ndir,1);
nfail3x=zeros(ndir,1);
nfail2x=zeros(ndir,1);
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
    [k_segment,std_k2,rce3d_k2,pre3d_k2,pro2d_k2,nfail3,nfail2,pct95_3d,pct95_2d,pct99_3d,pct99_2d] = ...
        simulation_fraction(npts(k),ko3d,conditions);
%  
     n_seg=n_seg + k_segment; % the number of trajectories increases by n_segment.
   
    if k_segment>=1
        std_dev2(k)= std_k2;
        rce3d(k) = rce3d_k2;
        aberror2(k)= pre3d_k2;
        ab_dev2(k)= pro2d_k2; 
        nfail3x(k) = nfail3;
        nfail2x(k) = nfail2;
        pc95_3d(k) = pct95_3d;
        pc95_2d(k) = pct95_2d;
        pc99_3d(k) = pct99_3d;
        pc99_2d(k) = pct99_2d;
        ndir_count = ndir_count + 1;
    else
        pc95_3d(k)=0;
        pc95_2d(k)=0;
        pc99_3d(k)=0;
        pc99_2d(k)=0;
        std_dev2(k)=0;
        rce3d(k) = 0;
        aberror2(k)=0;
        ab_dev2(k)=0;
        nfail3x(k) = 0;
        nfail2x(k) = 0;

%         avedisx(k) = 0;
%         nrepx(k) = 0;
    end
end % of k...
%% reporting results
std_dev=mean(sqrt(std_dev2(1:ndir_count))); % total average standard deviation of the patient.
rms=mean(sqrt(aberror2(1:ndir_count)));     % total average rms 3D error
devp=mean(sqrt(ab_dev2(1:ndir_count)));     % total average 2D error
estrms = mean(sqrt(rce3d(1:ndir_count)));
p95_3d=mean(pc95_3d(1:ndir_count));
p95_2d=mean(pc95_2d(1:ndir_count));
p99_3d=mean(pc99_3d(1:ndir_count));
p99_2d=mean(pc99_2d(1:ndir_count));
% p95_3d=0;
% p95_2d=0;
% p99_3d=0;
% p99_2d=0;
numberfail3=mean(nfail3x(1:ndir_count));
numberfail2=mean(nfail2x(1:ndir_count));

errors=[std_dev,estrms,rms,devp,p95_3d,p95_2d,p99_3d,p99_2d,numberfail3,numberfail2];
numbers = [ndir,nactual,n_seg];
xx=toc();
fidd=fopen('D:\Program Files\MATLAB\errorstatistics150929\error_analysis\stat_result3d_individual','a'); 
 
fprintf(fidd,' 2D beams-eye-view error: %4.2f patient %u time %4.2f \n',devp,patient_num,xx/n_seg); 
fclose(fidd);
patient_num=patient_num+1;
% indicate the analysis of one patient is finished
fprintf('the number of trajectories = %u \n', n_seg); 
fprintf('the number of fractions = %u \n', ndir_count); 
