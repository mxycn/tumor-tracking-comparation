function [std_dev2,rr2,numberfail2,perct95_2d,perct99_2d] = simulation_segment2(j0,so3d,conditions)
% simulate on one segment
% projecting and reconstructing the jth point.
%% initializing the parameters and arrays...
% is_pre = 1 if prediction is applied, 0 otherwise.
is_pre = conditions(1);

% is_block = 1 if block is simulated.
is_block = conditions(2);

% factor1 specifies the range of coordinate within which the points are
% selected for multivariate regression. factor1=selected range/full extent
% in that direction in the first 4 seconds.
factor1 = conditions(3);

% factor2 specifies the range of velocity within which the points are
% selected for multivariate regression. factor2=selected range/full extent
% in that direction in the first 4 seconds.
factor2 = conditions(4);

% factor3 specifies the ratio of the threshold controlling the predicted 
% changes of positions to the extent in their respective direction.
factor3 = conditions(5);

seg_err = conditions(6); % specifies the rms error of detection(=0.5mm).
lag_time = conditions(7); % the latent time
duration_ = conditions(10); % the time interval for one cycle of arc therapy (=72)
init_gAng = conditions(11); % initial angle(=-179)
end_gAng = conditions(12); % final angle(=+179)

% runlength: the number of initial points whose actual 3D coordinates
%  were measured in the first few seconds(~ 26)
init_time = conditions(13);
incrt = conditions(14);
runlength = floor(init_time/(incrt*0.0385))+1;
% one sample point are taken ever incrt data points, whose 
% period is approximately 0.038545s.incrt(= 4)
incrt = conditions(14);
sid = conditions(15); % sid: the Source-Imager distance.(=1500)
stdx = conditions(16); % stdx: the source-isocenter distance.(=1000)
r2threshold = conditions(23);
num_train = conditions(24);
marginx = conditions(25);

angular_speed=abs(end_gAng-init_gAng)/duration_; % the rate of gantry rotation 
ang_incrt=angular_speed*incrt*0.038545; % the angle interval between data points.
n_lag=round(lag_time/(incrt*38.5)); % the number of latent points
points_period=floor(358.0/ang_incrt);

r2d = zeros(points_period, 2);
o2d = zeros(points_period, 2);
r2do = zeros(points_period+n_lag, 2);
r2dop = zeros(points_period+n_lag, 2);
r2dp = zeros(points_period+n_lag, 2);
incr_2d = zeros(points_period, 2);

centerx = zeros(points_period,23);
linearity = zeros(points_period,1);
perror2d= zeros(points_period,1);
predictionquality=zeros(points_period,1);

rr2=0;
numberfail2=0;
%% Initializing
% add detection error to the points of the first few seconds, and calculate the 
% center point for the first few seconds.

%initializing the avarage
av2d=zeros(points_period,2);
av2d(1,:) = r2d(1,:);
if n_lag > 0.1
    r2dp(1:n_lag,:)=NaN;
end
jj=1;
%% start simulate the position estimation and prediction point by point
while jj <= points_period
    % project who
    % start timing
    gAng=init_gAng+ang_incrt*(jj-1);    
        [o2d(jj,:) r2d(jj,:)] = project_to_imager(seg_err,so3d,jj,gAng,sid,stdx);
        if jj>2
            av2d(jj-1,:)=(r2d(jj-2,:)+2*r2d(jj-1,:)+r2d(jj,:))/4;  
        else
            av2d(1,:) = r2d(1,:);            
        end        
            
        if jj> runlength-n_lag % start to predict with calculation
            [extentx,extentv]=near_range2d(jj,runlength,0,r2d);
            [incr_2d(jj,:),predictionquality(jj)]=predicting2(is_pre,factor1,factor2,r2threshold,num_train,n_lag,jj,extentx,extentv,av2d);    
            [r2dp(jj+n_lag,:),predq]=prediction_screening2(jj+n_lag,n_lag,r2d,incr_2d,extentv,factor3);
            perror2d(jj)=norm(r2dp(jj,:)-o2d(jj,:))*stdx/sid;
        elseif jj > n_lag % start to predict without real calculation
            r2dp(jj,:)=r2d(jj-n_lag,:); 
            perror2d(jj)=norm(r2dp(jj,:)-o2d(jj,:))*stdx/sid;
        end % predicting for the first few seconds
    rr2=rr2+perror2d(jj)*perror2d(jj);
    if perror2d(jj)>5
        numberfail2=numberfail2+1;
    end
    jj=jj+1;
end % of "while jj< = point_period."
% calculate the variation of the trajectory itself.
od=zeros(points_period,3);
% construct a new array for calculating the standard deviation
% of the trajectory.
od(:,:)=so3d(1:points_period,:);    
% sum up the square of all standard deviation.
variancex = std(od(1:points_period,1))*std(od(1:points_period,1))+...
    std(od(1:points_period,2))*std(od(1:points_period,2));
std_dev2 = variancex; 
rr2 = rr2/(points_period-n_lag);

%% plotting
% xcoord=zeros(points_period,1);
% for i=1:points_period
%     xcoord(i)=ang_incrt*(i-1)+init_gAng+180;
%     r2do(i,:) = r2d(i,:)*stdx/sid;
%     r2dop(i,:) = r2dp(i,:)*stdx/sid;
% %     norm3de(i)=sqrt(perror2d(i,1)*perror2d(i,1)+perror2d(i,2)*perror2d(i,2)+perror2d(i,3)*perror2d(i,3));
% %     normdev(i) = sqrt(displs(i,1)*displs(i,1)+displs(i,2)*displs(i,2)+displs(i,3)*displs(i,3));
% %     if normdev(i)>threshold + 2
% %         numberfail=numberfail+1;
% %     end
% %     if iskV(i)
% %         numberkV=numberkV+1;
% %     end
% %     if isrep(i)
% %         numberrep = numberrep + 1;
% %     end
% %     meandis = meandis + normdev(i);
% end
% meandis = meandis/points_period;
% 
% % std1=norm(perror2d(j0+n_lag:j0+points_period-1))/sqrt(points_period-n_lag); % average 2D error with prediction
% % max1=max(perror2d(j0+runlength+2:j0+points_period-1));
% % std2=norm(perror2dx(j0+n_lag:j0+points_period-1))/sqrt(points_period-n_lag); 
% % std3=norm(recon_error2d(j0+runlength:j0+points_period-1))/sqrt(points_period-runlength); % average 3D error of pure estimation
perct95_2d=sqrt(get_percentile_error(perror2d,0.95,points_period));
perct99_2d=sqrt(get_percentile_error(perror2d,0.99,points_period));

% % conditionx=std_dev2>3.9 & std_dev2<11.5;  
% % % conditionx= (j0==2353 | j0==4217 | j0==5615); 
% % modex = 3;
% n_segment = 1;
% j0=n_lag+1;
%         figure(n_segment), subplot(3,1,1)
%         set(gcf,'Position',[5,30,1400,800])
%         plot(xcoord(j0:points_period),o2d(j0:points_period,1),'r',... 
%             xcoord(j0:points_period),r2do(j0:points_period,1),'b',...
%             xcoord(j0:points_period),r2dop(j0:points_period,1),'c','LineWidth',1.5)
%         set(gca,'FontSize',16)
%         grid on
%         ylabel('Position(mm)'), xlim([1 359]),...
%             legend('original SI','plate SI 2','predicted SI 2');
% 
% %             title(strcat('original, estimated, and residual positions of - (j0)',...
% %             dbfolders(k).name,', j0=',int2str(j0), ', 2D rms error=', num2str(std1)),'FontSize',18)
%         
% %         figure(n_segment), subplot(4,1,3)
% %         plot(xcoord(j0:points_period),so3d(j0:points_period,3),'r',... 
% %             xcoord(j0:points_period),r3dr(j0:points_period,3),'b',...
% %             xcoord(j0:points_period),r3d(j0:points_period,3),'g','LineWidth',1.5)
% %         set(gca,'FontSize',16)
% %         grid()
% %         ylabel('Position(mm)'),xlim([1 359]),...
% %             legend('original AP','residual AP', 'checkpoint AP');
%         figure(n_segment), subplot(3,1,3)
%         plot(xcoord(j0:points_period),perror2d(j0:points_period),'k','LineWidth',1.5)
%         set(gca,'FontSize',16)
%    
%         ylabel('Error(mm)'), xlabel('angles in degree'),...
%             xlim([1 359]);
%         hold on
