function patient_selection(conditions)
% this function specifies different patients, summarizes the calculation
% results, mainly the errors.

npatient=0;  % number of patients calculated
napatients = 0; % number of patients considered
numbers = zeros(46,3); % [ndir,nactual,n_seg];
errors = zeros(46,10);
% errors = zeros(46,4);
incrt = conditions(14);
durationx = conditions(10); % the time interval for one cycle of arc therapy (=72)
init_gAng = conditions(11); % initial angle(=-179)
end_gAng = conditions(12); % final angle(=+179)
angular_speed=abs(end_gAng-init_gAng)/durationx; % the rate of gantry rotation 
ang_incrt=angular_speed*incrt*0.038545; % the angle interval between data points.
points_period=floor(358.0/ang_incrt);   % the number of points in a trajectory.    
global patient_num
patient_num=1;
% By commenting and uncommenting, we can select the patient to be analyzed in the study.  
% [numbers(1,:),errors(1,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB01');
% [numbers(2,:),errors(2,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB02');
% [numbers(3,:),errors(3,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB03');
% [numbers(4,:),errors(4,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB04');
% [numbers(5,:),errors(5,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB05');
% [numbers(6,:),errors(6,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB06');
% [numbers(7,:),errors(7,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB07');
% [numbers(8,:),errors(8,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB08');
% [numbers(9,:),errors(9,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB09');
% [numbers(10,:),errors(10,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB10');
% [numbers(11,:),errors(11,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB11');
% [numbers(12,:),errors(12,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB12');
% [numbers(13,:),errors(13,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB13');
% [numbers(14,:),errors(14,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB14');
% [numbers(15,:),errors(15,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB15');
% [numbers(17,:),errors(17,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB17');
% [numbers(18,:),errors(18,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB18');
% [numbers(19,:),errors(19,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB19');
% [numbers(20,:),errors(20,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB20');
% [numbers(21,:),errors(21,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB21');
% [numbers(22,:),errors(22,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB22');
% [numbers(23,:),errors(23,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB23');
% [numbers(25,:),errors(25,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB25');
% [numbers(26,:),errors(26,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB26');
% [numbers(27,:),errors(27,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB27');
% [numbers(28,:),errors(28,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB28');
% [numbers(29,:),errors(29,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB29');
% [numbers(30,:),errors(30,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB30');
% [numbers(31,:),errors(31,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB31');
% [numbers(32,:),errors(32,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB32');
% [numbers(33,:),errors(33,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB33');
% [numbers(34,:),errors(34,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB34');
% [numbers(35,:),errors(35,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB35');
% [numbers(36,:),errors(36,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB36');
% [numbers(37,:),errors(37,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB37');
% [numbers(39,:),errors(39,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB39');
% [numbers(40,:),errors(40,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB40');
% [numbers(41,:),errors(41,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB41');
% [numbers(42,:),errors(42,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB42');
% [numbers(43,:),errors(43,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB43');
% [numbers(45,:),errors(45,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB45');
% [numbers(46,:),errors(46,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB46');
[numbers(1,:),errors(1,:)]=simulation_patient(conditions,'H:\database\DB01');
[numbers(2,:),errors(2,:)]=simulation_patient(conditions,'H:\database\DB02');
[numbers(3,:),errors(3,:)]=simulation_patient(conditions,'H:\database\DB03');
[numbers(4,:),errors(4,:)]=simulation_patient(conditions,'H:\database\DB04');
[numbers(5,:),errors(5,:)]=simulation_patient(conditions,'H:\database\DB05');
[numbers(6,:),errors(6,:)]=simulation_patient(conditions,'H:\database\DB06');
[numbers(7,:),errors(7,:)]=simulation_patient(conditions,'H:\database\DB07');
[numbers(8,:),errors(8,:)]=simulation_patient(conditions,'H:\database\DB08');
[numbers(9,:),errors(9,:)]=simulation_patient(conditions,'H:\database\DB09');
[numbers(10,:),errors(10,:)]=simulation_patient(conditions,'H:\database\DB10');
[numbers(11,:),errors(11,:)]=simulation_patient(conditions,'H:\database\DB11');
[numbers(12,:),errors(12,:)]=simulation_patient(conditions,'H:\database\DB12');
[numbers(13,:),errors(13,:)]=simulation_patient(conditions,'H:\database\DB13');
[numbers(14,:),errors(14,:)]=simulation_patient(conditions,'H:\database\DB14');
[numbers(15,:),errors(15,:)]=simulation_patient(conditions,'H:\database\DB15');
[numbers(17,:),errors(17,:)]=simulation_patient(conditions,'H:\database\DB17');
[numbers(18,:),errors(18,:)]=simulation_patient(conditions,'H:\database\DB18');
[numbers(19,:),errors(19,:)]=simulation_patient(conditions,'H:\database\DB19');
[numbers(20,:),errors(20,:)]=simulation_patient(conditions,'H:\database\DB20');
[numbers(21,:),errors(21,:)]=simulation_patient(conditions,'H:\database\DB21');
[numbers(22,:),errors(22,:)]=simulation_patient(conditions,'H:\database\DB22');
[numbers(23,:),errors(23,:)]=simulation_patient(conditions,'H:\database\DB23');
[numbers(25,:),errors(25,:)]=simulation_patient(conditions,'H:\database\DB25');
[numbers(26,:),errors(26,:)]=simulation_patient(conditions,'H:\database\DB26');
[numbers(27,:),errors(27,:)]=simulation_patient(conditions,'H:\database\DB27');
[numbers(28,:),errors(28,:)]=simulation_patient(conditions,'H:\database\DB28');
[numbers(29,:),errors(29,:)]=simulation_patient(conditions,'H:\database\DB29');
[numbers(30,:),errors(30,:)]=simulation_patient(conditions,'H:\database\DB30');
[numbers(31,:),errors(31,:)]=simulation_patient(conditions,'H:\database\DB31');
[numbers(32,:),errors(32,:)]=simulation_patient(conditions,'H:\database\DB32');
[numbers(33,:),errors(33,:)]=simulation_patient(conditions,'H:\database\DB33');
[numbers(34,:),errors(34,:)]=simulation_patient(conditions,'H:\database\DB34');
[numbers(35,:),errors(35,:)]=simulation_patient(conditions,'H:\database\DB35');
[numbers(36,:),errors(36,:)]=simulation_patient(conditions,'H:\database\DB36');
[numbers(37,:),errors(37,:)]=simulation_patient(conditions,'H:\database\DB37');
[numbers(39,:),errors(39,:)]=simulation_patient(conditions,'H:\database\DB39');
[numbers(40,:),errors(40,:)]=simulation_patient(conditions,'H:\database\DB40');
[numbers(41,:),errors(41,:)]=simulation_patient(conditions,'H:\database\DB41');
[numbers(42,:),errors(42,:)]=simulation_patient(conditions,'H:\database\DB42');
[numbers(43,:),errors(43,:)]=simulation_patient(conditions,'H:\database\DB43');
[numbers(45,:),errors(45,:)]=simulation_patient(conditions,'H:\database\DB45');
[numbers(46,:),errors(46,:)]=simulation_patient(conditions,'H:\database\DB46');
% [numbers(15,:),errors(15,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB15');
% [numbers(17,:),errors(17,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB17');
% [numbers(20,:),errors(20,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB20');
% [numbers(22,:),errors(22,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB22');
% [numbers(29,:),errors(29,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB29');
% [numbers(33,:),errors(33,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB33');
% [numbers(42,:),errors(42,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB42');

% [numbers(1,:),errors(1,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB01');
% [numbers(5,:),errors(5,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB05');
% [numbers(6,:),errors(6,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB06');
% [numbers(10,:),errors(10,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB10');
% [numbers(15,:),errors(15,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB15');
% [numbers(17,:),errors(17,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB17');
% [numbers(18,:),errors(18,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB18');
% [numbers(20,:),errors(20,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB20');
% [numbers(22,:),errors(22,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB22');
% [numbers(27,:),errors(27,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB27');
% [numbers(28,:),errors(28,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB28');
% [numbers(29,:),errors(29,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB29');
% [numbers(33,:),errors(33,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB33');
% [numbers(39,:),errors(39,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB39');
% [numbers(42,:),errors(42,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB42');
% [numbers(43,:),errors(43,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB43');
% [numbers(45,:),errors(45,:)]=simulation_patient(conditions,'e:\research\Liuwu\Database\DB45');
%% The following is to calculate the statistics for many patients
for ii = 1:46
    if numbers(ii,1) > 0
        napatients =napatients + 1;
        if numbers(ii,3) > 0
            npatient = npatient +1;
        end    
    end
end
p95_3d=0;
p95_2d=0;
p99_3d=0;        
p99_2d=0;
rmsx=0;        % patient-averaged 3D position rms error
devp=0;       % patient-averaged 2D rms error
stdd=0;  % patient-averaged standard deviation of the data points
estrmsx =0;
nfail3=0;
nfail2=0;
na_total=0;   % the number of total data point available
% totaltime=[0 0];
n_segments=0;
% sum up the errors and other statistics.
for jj=1:46
    na_total=na_total+numbers(jj,2);      
    n_segments=n_segments + numbers(jj,3);
    stdd=stdd+errors(jj,1); 
    estrmsx = estrmsx + errors(jj,2); 
    rmsx=rmsx+errors(jj,3); 
    devp=devp+errors(jj,4); 
    p95_3d=p95_3d+errors(jj,5);
    p95_2d=p95_2d+errors(jj,6);
    p99_3d=p99_3d+errors(jj,7);      
    p99_2d=p99_2d+errors(jj,8); 
    nfail3=nfail3+errors(jj,9); 
    nfail2=nfail2+errors(jj,10); 
%     num_cord=num_cord+num_plotx(jj);
%     count95_3=count95_3+count3(jj);
%     totaltime=totaltime+timex(jj,:);
end

% calculate the averages.
rmsx=rmsx/npatient; % patient-averaged
devp=devp/npatient; % patient-averaged
stdd=stdd/npatient; % patient-averaged
estrmsx=estrmsx/npatient; % patient-averaged
nfail3=nfail3/npatient; % patient-averaged
nfail2=nfail2/npatient; % patient-averaged
% 
p95_3d=p95_3d/npatient;
p95_2d=p95_2d/npatient;
p99_3d=p99_3d/npatient;
p99_2d=p99_2d/npatient; 

percentage=n_segments*points_period/na_total;
% 
fidd=fopen('D:\Program Files\MATLAB\errorstatistics150929\error_analysis\stat_result3d','a');
% 
is_pre = conditions(1);
is_block = conditions(2);
factor1 = conditions(3);
factor2 = conditions(4);
factor3 = conditions(5);
seg_err = conditions(6);
lag_time = conditions(7);
linearity = conditions(8);
amp = conditions(9); 
% init_time = conditions(13);
% num_train = conditions(24);
% marginx = conditions(25);

fprintf(fidd,'Beginning - parameters: \n  is_pre: %u, is_block: %u, factor1: %3.1f, factor2: %3.1f, factor3: %3.1f, seg_err: %3.1f,\n  lag_time: %u, linearity constraint: %4.2f, amplitude variation: %u \n',is_pre,is_block,factor1,factor2,factor3,seg_err,lag_time,linearity,amp);  

for jj=1:46
    fprintf(fidd,'No.(%d) std: %4.2f est rms: %4.2f pred rms: %4.2f devp: %4.2f; p95_3d: %4.2f p95_2d: %4.2f p99_3d: %4.2f p99_2d: %4.2f nfail3: %4.2f nfail2: %4.2f;n_seg: %u\n',...
        jj, errors(jj,1),errors(jj,2),errors(jj,3),errors(jj,4),errors(jj,5),errors(jj,6),errors(jj,7),errors(jj,8),errors(jj,9),errors(jj,10),numbers(jj,3));
end

fprintf(fidd,'results: n_segments=%u, percentage=%4.2f, std of the positions: %4.3f;\n Patient_averaged: 3D estimation error, %4.3f, 3D rms error: %4.3f, 2D rms error: %4.3f \n',n_segments,percentage,stdd,estrmsx,rmsx,devp); 
fprintf(fidd,' percentile error:  p95_3d=%4.2f, p95_2d=%4.2f, p99_3d=%4.2f, p99_2d=%4.2f, nfail3: %4.2f  nfail2: %4.2f \n end \n',p95_3d,p95_2d,p99_3d,p99_2d,nfail3,nfail2); 
fclose(fidd);  
return