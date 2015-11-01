function [std_dev2,rcerror2,aberror2,ab_dev2,numberfail3,numberfail2,perct95_3d,perct95_2d,perct99_3d,perct99_2d] = simulation_segment(j0,so3d,conditions)
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

r3d = zeros(points_period, 3);
r3dr = zeros(points_period+n_lag, 3);
r3dp = zeros(points_period+n_lag, 3);
incr_3d = zeros(points_period, 3);

kr = zeros(points_period,3); % the reconstructed direction vector 
ko = zeros(points_period,3); % the original direction vector
ym = zeros(points_period,2); % the end points of the PCA axis
checkpoint = zeros(points_period,3); % the reference of repositioning.
centerx = zeros(points_period,3);
linearity = zeros(points_period,1);
iskV = false(points_period,1);
rerror3d = zeros(points_period,3);
perror3d= zeros(points_period,3);
perror2d= zeros(points_period,2);
is_plane = zeros(points_period,1);
predictionquality=zeros(points_period,1);
Bc = zeros(points_period,3);
cost = zeros(points_period,1);

rcerror2=0;
aberror2=0;
ab_dev2=0;
numberfail3=0;
numberfail2=0;
deviat = zeros(points_period,3);
displs = zeros(points_period,3);
%% Initializing
% add detection error to the points of the first few seconds, and calculate the 
% center point for the first few seconds.

%initializing the avarage
av3d=zeros(points_period,3);
av3d(1,:) = r3d(1,:);
if n_lag > 0.1
    r3dr(1:n_lag,:)=NaN;
    r3dp(1:n_lag,:)=NaN;
end
jj=1;
%% start simulate the position estimation and prediction point by point
while jj <= points_period
    % project who
    % start timing
    gAng=init_gAng+ang_incrt*(jj-1);    
    if jj<=runlength
        iskV(jj) = true;
        r3d(jj,:)=so3d(jj,:) + seg_err * stdx/sid * randn(1,3);  
        if jj>2
            av3d(jj-1,:)=(r3d(jj-2,:)+2*r3d(jj-1,:)+r3d(jj,:))/4;  
        else
            av3d(1,:) = r3d(1,:);            
        end        
            
        if jj > n_lag
            r3dr(jj,:)=r3d(jj-n_lag,:); 
            [rerror3d(jj,:),perror3d(jj,:),perror2d(jj,:)]=errorv_calculation(jj,gAng,sid,stdx,so3d,r3d,r3dr); 
        end % predicting for the first few seconds
        if jj>runlength-n_lag
            [extentx,extentv]=near_range3d(jj,runlength,0,r3d);                     
            % calculate the difference between the point n_log
            % before and the current point, for prediction. Note the
            % input information is from the points at least n_lag
            % before.
            [incr_3d(jj,:),predictionquality(jj)]=predicting(is_pre,factor1,factor2,r2threshold,num_train,n_lag,jj,extentx,extentv,av3d);    
            
            % Add the difference to the position n_lag points before 
            % to get the predicted position of the current point.
            % If the calculated difference is too large, replace
            % it with a threshold value, which is based on the
            % amplitude of motion in the first 4 seconds. 
            [r3dp(jj+n_lag,:),predq]=prediction_screening(jj+n_lag,n_lag,r3d,incr_3d,extentv,factor3);
            predictionquality(jj) = min(predictionquality(jj),predq);
            % optimize the prediction result with the updated 
            % correlation information n_lag points before
            r3dr(jj+n_lag,:)=r3dp(jj+n_lag,:);
        end
    else  % of jj<=runlength. estimating position and prediction should be made.
%         tic();
        
        [r3dx,centerx(jj-1,:),kr(jj-1,:),ko(jj-1,:),ym(jj-1,:),linearity(jj-1),is_plane(jj-1), Bc(jj,:),cost(jj)]= ...
            proj_backproj_new(seg_err,marginx,so3d,r3d,jj,gAng,runlength,sid,stdx);        
%         xx1=toc();
%         timex(1)=timex(1)+xx1;
%         if xx1>0.010
%             fprintf('extra long estimation time:%f \n', xx1);
%         end
        if iskV(jj)  
            r3d(jj,:)=so3d(jj,:) + seg_err * stdx/sid * randn(1,3);              
        else % accepting the returned values.
            r3d(jj,:)=r3dx(:);   
        end
        
       % smoothing the position
       av3d(jj-1,:)=(r3d(jj-2,:)+2*r3d(jj-1,:)+r3d(jj,:))/4;     
%         av3d(jj-n_lag,:)=(2*r3d(jj-n_lag-1,:)+3*r3d(jj-n_lag,:))/5;
            
      %% start to predict ..
        if n_lag>0
%             tic();
            % calculate the immediate past range
            [extentx,extentv]=near_range3d(jj,runlength,0,r3d);                     
            % calculate the difference between the point n_log
            % before and the current point, for prediction. Note the
            % input information is from the points at least n_lag
            % before.
            [incr_3d(jj,:),predictionquality(jj)]=predicting(is_pre,factor1,factor2,r2threshold,num_train,n_lag,jj,extentx,extentv,av3d);    
            
            % Add the difference to the position n_lag points before 
            % to get the predicted position of the current point.
            % If the calculated difference is too large, replace
            % it with a threshold value, which is based on the
            % amplitude of motion in the first 4 seconds. 
            [r3dp(jj+n_lag,:),predq]=prediction_screening(jj+n_lag,n_lag,r3d,incr_3d,extentv,factor3);
            predictionquality(jj) = min(predictionquality(jj),predq);
            % optimize the prediction result with the updated 
            % correlation information n_lag points before 
            if jj>1+runlength-n_lag
                r3dr(jj+n_lag,:)=prediction_optimizing_new(jj+n_lag,r3dp,stdx,centerx(jj-1,:),kr(jj-1,:),linearity(jj-1)); % true realtime  
            else
                r3dr(jj+n_lag,:)=r3dp(jj+n_lag,:);
            end
            
            % the following is the simulation of blockage and the
            % update of estimated coordinates and centerpoints.
%             deviat(jj+n_lag,:) = r3dr(jj+n_lag,:) - checkpoint(jj,:);
%             displs(jj,:) = so3d(jj,:) - checkpoint(jj,:);
            [rerror3d(jj,:),perror3d(jj,:),perror2d(jj,:)]=errorv_calculation(jj,gAng,sid,stdx,so3d,r3d,r3dr);
%             [rerror3d(jj+n_lag,:) perror3d(jj+n_lag,:) perror2d(jj+n_lag)]=errorv_calculation(jj+n_lag,gAng,sid,stdx,so3d,r3d,r3dr);
                  
%             cos1 = Bc(jj,:)*kr(jj-1,:)'; % cosine of the angle between the principal axis and the treatment beam
% %         cos2 = Bc(jj,:)*(r3dr(jj,:)-r3dr(jj-1,:))'/speed; % cosine of the angle between the velocity axis and the beam
%             cos1s = cos1*cos1;
%             sqrte = sqrt(norm(perror2d(jj,:)));
% % %             sqrte = norm(perror2d(jj,:));
% % %             err_unc(1) = 0.45 + 0.41*abs(cos1) -0.66*cos1s + 0.49* sqrte;
% % %             err_unc(2) = 0.81 - 2.40*abs(cos1) +2.70*cos1s + 0.24* sqrte + 0.82*abs(cos1)*sqrte;
% % %             err_unc(3) = 0.73 - 2.51*abs(cos1) +2.41*cos1s + 0.16* sqrte + 0.9*abs(cos1)*sqrte;
% % % %             err_unc(1) = 0.44 + 0.22*abs(cos1) - 0.32*cos1s - 0.28* sqrte + 0.5*sqrte*sqrte;
% % % %             err_unc(2) = 0.48 - 1.49*abs(cos1) + 2.52*cos1s - 0.19* sqrte + 0.56*sqrte*sqrte;
% % % %             err_unc(3) = 0.41 - 1.27*abs(cos1) + 1.87*cos1s - 0.09* sqrte + 0.52*sqrte*sqrte;
%             
%             err_unc(1) = 0.41 + 0.21*abs(cos1) - 0.29*cos1s - 0.21* sqrte + 0.45*sqrte*sqrte;
%             err_unc(2) = 0.51 - 1.70*abs(cos1) + 2.79*cos1s - 0.30* sqrte + 0.67*sqrte*sqrte;
%             err_unc(3) = 0.44 - 1.24*abs(cos1) + 1.84*cos1s - 0.19* sqrte + 0.59*sqrte*sqrte;
%             
% %             if norm(err_unc) + 0.5*norm(deviat(jj+n_lag,:)) > threshold
% %                 iskV(jj+n_lag) = true;   
% %                 iskV(jj+n_lag) = true;    
% %             end
%             if iskvon && norm(err_unc) > threshold
%                 iskV(jj+n_lag) = true;  
%                 iskV(jj+n_lag+1) = true;  
%             end
%             checkpoint(jj+n_lag,:) = r3dr(jj+n_lag,:);
%             if isrep(jj)
% %             if iskV(jj)
%                 checkpoint(jj+n_lag,:) = r3dr(jj+n_lag,:);
% %                 checkpoint(jj+n_lag,:) = (r3dr(jj+n_lag,:)+r3d(jj,:))/2;
%             else
%                 checkpoint(jj+n_lag,:) = checkpoint(jj+n_lag-1,:);
% %                 checkpoint(jj+n_lag,:) = r3d(jj,:);
%             end
            
%             xx2=toc();
%             timex(2)=timex(2)+xx2;
%             if xx2>0.02
%                 fprintf('extra long prediction time:%f \n', xx2);
%             end
        else % of n_lag>0. No latency is considered
            % if the point position is measured in the first 4
            % seconds, or no latency is considered.
            % no prediction is made, use the latest point estimated  
            r3dr(jj,1:3)=r3d(jj,1:3);
            % calculate the projection of estimated position.
%             [rerror3d(jj) perror3d(jj) perror2d(jj)]=error_calculation(jj,gAng,sid,stdx,so3d,r3d,r3dr);
            [rerror3d(jj,:),perror3d(jj,:),perror2d(jj,:)]=errorv_calculation(jj,gAng,sid,stdx,so3d,r3d,r3dr);
        end % of n_lag>0.
    end % of jj<=runlength
    %         fprintf('extra long prediction time22:%f %u  \n',  r3dr(jj,1), jj);
    rr3=perror3d(jj,1)*perror3d(jj,1)+perror3d(jj,2)*perror3d(jj,2)+...
        perror3d(jj,3)*perror3d(jj,3);
    aberror2 = aberror2+rr3;
    rcerror2 = rcerror2+rerror3d(jj,1)*rerror3d(jj,1)+rerror3d(jj,2)*rerror3d(jj,2)+...
        rerror3d(jj,3)*rerror3d(jj,3);
    rr2=perror2d(jj,1)*perror2d(jj,1)+perror2d(jj,2)*perror2d(jj,2);   
    ab_dev2=ab_dev2+rr2;
    if rr3 > 25
        numberfail3=numberfail3+1;
    end
    if rr2>25
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
    std(od(1:points_period,2))*std(od(1:points_period,2))+...
    std(od(1:points_period,3))*std(od(1:points_period,3));
std_dev2 = variancex; 
rcerror2 = rcerror2/(points_period-n_lag); 
aberror2 = aberror2/(points_period-n_lag);
ab_dev2 = ab_dev2/(points_period-n_lag);

%% plotting
% xcoord=zeros(points_period,1);
% for i=1:points_period
%     xcoord(i)=ang_incrt*(i-1)+init_gAng+180;
%     norm3de(i)=sqrt(perror3d(i,1)*perror3d(i,1)+perror3d(i,2)*perror3d(i,2)+perror3d(i,3)*perror3d(i,3));
%     normdev(i) = sqrt(displs(i,1)*displs(i,1)+displs(i,2)*displs(i,2)+displs(i,3)*displs(i,3));
%     if normdev(i)>threshold + 2
%         numberfail=numberfail+1;
%     end
%     if iskV(i)
%         numberkV=numberkV+1;
%     end
% %     if isrep(i)
% %         numberrep = numberrep + 1;
% %     end
%     meandis = meandis + normdev(i);
% end
% meandis = meandis/points_period;
% 
% % std1=norm(perror2d(j0+n_lag:j0+points_period-1))/sqrt(points_period-n_lag); % average 2D error with prediction
% % max1=max(perror2d(j0+runlength+2:j0+points_period-1));
% % std2=norm(perror3dx(j0+n_lag:j0+points_period-1))/sqrt(points_period-n_lag); 
% % std3=norm(recon_error3d(j0+runlength:j0+points_period-1))/sqrt(points_period-runlength); % average 3D error of pure estimation

perct95_3d=sqrt(get_percentile_error(perror3d,0.95,points_period));
perct95_2d=sqrt(get_percentile_error(perror2d,0.95,points_period));
perct99_3d=sqrt(get_percentile_error(perror3d,0.99,points_period));
perct99_2d=sqrt(get_percentile_error(perror2d,0.99,points_period));


% perct95_3d=0;
% perct95_2d=0;
% perct99_3d=0;
% perct99_2d=0;

% % conditionx=std_dev2>3.9 & std_dev2<11.5 & numberkV > 32;  
% conditionx= (j0==2353 | j0==4217 | j0==5615); 
% modex = 3;
% n_segment = 1;
% % plotfigures(conditionx,modex,is_block,n_segment,points_period,xcoord,so3d,r3dr,r3d,sqrt(perror2d),sqrt(perror3d),sqrt(rerror3d));
% % plotfigures(1,conditionx,modex,is_block,n_segment,points_period,xcoord,so3d,deviat,r3d,normdev,norm3de,perror2d);
% plotfiguresx(j0,conditionx,modex,is_block,n_segment,points_period,xcoord,so3d,displs,checkpoint,iskV,normdev);

% modex = 2;
% n_segment = 2;
% plotfigures(1,conditionx,modex,is_block,n_segment,points_period,xcoord,so3d,displs,r3d,iskV,normdev,norm3de);
% plotfigures(conditionx,modex,is_block,n_segment,points_period,xcoord,so3d,r3dr,r3d,sqrt(perror2d),sqrt(perror3d),sqrt(rerror3d));
% plotfigures(1,conditionx,modex,is_block,n_segment,points_period,xcoord,so3d,r3dr,r3d,normdev,norm3de,perror2d);



    
        
    