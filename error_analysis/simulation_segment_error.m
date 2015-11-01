function [std_dev2,aberror2,ab_dev2,parameters,threed,twod] = simulation_segment_error(so3d,conditions)
% simulate on one segment
% projecting and reconstructing the jth point.
% to pool the information on prediction error and its potential predictors.
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
err_unc = zeros(3,1);

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

count95_3 = 0;
count95_4 = 0;

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
rerror3d = zeros(points_period,3);
perror3d= zeros(points_period,3);
perror2d= zeros(points_period,2);
is_plane = zeros(points_period,1);
predictionquality=zeros(points_period,1);
Bc = zeros(points_period,3);
cost = zeros(points_period,1);

aberror2=0;
ab_dev2=0;
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
        r3d(jj,:)=so3d(jj,:) + seg_err * stdx/sid * randn(1,3);  
        if jj>2
            av3d(jj-1,:)=(r3d(jj-2,:)+2*r3d(jj-1,:)+r3d(jj,:))/4;  
        else
            av3d(1,:) = r3d(1,:);            
        end        
            
        if jj > n_lag
            r3dr(jj,:)=r3d(jj-n_lag,:); 
            [rerror3d(jj,:) perror3d(jj,:) perror2d(jj,:)]=errorv_calculation(jj,gAng,sid,stdx,so3d,r3d,r3dr); 
        end % predicting for the first few seconds
        if jj>runlength-n_lag
            [extentx extentv]=near_range3d(jj,runlength,0,r3d);                     
            % calculate the difference between the point n_log
            % before and the current point, for prediction. Note the
            % input information is from the points at least n_lag
            % before.
            [incr_3d(jj,:) predictionquality(jj)]=predicting(is_pre,factor1,factor2,r2threshold,num_train,n_lag,jj,extentx,extentv,av3d);    
            
            % Add the difference to the position n_lag points before 
            % to get the predicted position of the current point.
            % If the calculated difference is too large, replace
            % it with a threshold value, which is based on the
            % amplitude of motion in the first 4 seconds. 
            [r3dp(jj+n_lag,:) predq]=prediction_screening(jj+n_lag,n_lag,r3d,incr_3d,extentv,factor3);
            predictionquality(jj) = min(predictionquality(jj),predq);
            % optimize the prediction result with the updated 
            % correlation information n_lag points before 
            if jj>1+runlength-n_lag
                r3dr(jj+n_lag,:)=prediction_optimizing_new(jj+n_lag,r3dp,stdx,centerx(jj-n_lag,:),kr(jj-n_lag,:),linearity(jj-n_lag)); % true realtime  
            else
                r3dr(jj+n_lag,:)=r3dp(jj+n_lag,:);
            end
        end
    else  % of jj<=runlength. estimating position and prediction should be made.
        r3dx = [0 0 0];
        tic();
        
        [r3dx,centerx(jj-1,:),kr(jj-1,:),ko(jj-1,:),ym(jj-1,:),linearity(jj-1),is_plane(jj-1), Bc(jj,:),cost(jj)]= ...
            proj_backproj_new(seg_err,marginx,so3d,r3d,jj,gAng,runlength,sid,stdx);        
        xx1=toc();
%         timex(1)=timex(1)+xx1;
        if xx1>0.010
            fprintf('extra long estimation time:%f \n', xx1);
        end
        % predicting the position the calculation is done starting from the j0+n_lag, the
        % first n_lag points are neglected.
%         if iskV(jj)  
%             r3d(jj,:)=so3d(jj,:) + seg_err * stdx/sid * randn(1,3);              
%         else % accepting the returned values.
        r3d(jj,:)=r3dx(:);   
%         end
       % smoothing the position
       av3d(jj-1,:)=(r3d(jj-2,:)+2*r3d(jj-1,:)+r3d(jj,:))/4;     
%         av3d(jj-n_lag,:)=(2*r3d(jj-n_lag-1,:)+3*r3d(jj-n_lag,:))/5;
            
      %% start to predict ..
        if n_lag>0
            tic();
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
            [rerror3d(jj,:) perror3d(jj,:) perror2d(jj,:)]=errorv_calculation(jj,gAng,sid,stdx,so3d,r3d,r3dr);
%             [rerror3d(jj+n_lag,:) perror3d(jj+n_lag,:) perror2d(jj+n_lag)]=errorv_calculation(jj+n_lag,gAng,sid,stdx,so3d,r3d,r3dr);
                  
            cos1 = Bc(jj,:)*kr(jj-1,:)'; % cosine of the angle between the principal axis and the treatment beam
%         cos2 = Bc(jj,:)*(r3dr(jj,:)-r3dr(jj-1,:))'/speed; % cosine of the angle between the velocity axis and the beam
            cos1s = cos1*cos1;
            sqrte = sqrt(norm(perror2d(jj,:)));
            xx2=toc();
%             timex(2)=timex(2)+xx2;
            if xx2>0.02
                fprintf('extra long prediction time:%f \n', xx2);
            end
        else % of n_lag>0. No latency is considered
            % if the point position is measured in the first 4
            % seconds, or no latency is considered.
            % no prediction is made, use the latest point estimated  
            r3dr(jj,1:3)=r3d(jj,1:3);
            % calculate the projection of estimated position.
%             [rerror3d(jj) perror3d(jj) perror2d(jj)]=error_calculation(jj,gAng,sid,stdx,so3d,r3d,r3dr);
            [rerror3d(jj,:) perror3d(jj,:) perror2d(jj,:)]=errorv_calculation(jj,gAng,sid,stdx,so3d,r3d,r3dr);
        end % of n_lag>0.
    end % of jj<=runlength
    %         fprintf('extra long prediction time22:%f %u  \n',  r3dr(jj,1), jj);
    aberror2=aberror2+perror3d(jj,1)*perror3d(jj,1)+perror3d(jj,2)*perror3d(jj,2)+...
        perror3d(jj,3)*perror3d(jj,3);
    ab_dev2=ab_dev2+perror2d(jj,1)*perror2d(jj,1)+perror2d(jj,2)*perror2d(jj,2);   
        if jj == 1 + runlength
%         r3derr = rerror3d(jj,:);
%         p3derr = perror3d(jj,:);      
%         p2derr = perror2d(jj,:);      
        threed = norm(rerror3d(jj,:));
        twod = norm(perror2d(jj,:));   
        speed = norm(r3dr(jj-1,:)-r3dr(jj-3,:));
        cos1 = Bc(jj,:)*kr(jj-2,:)'; % cosine of the angle between the principal axis and the treatment beam
%         cos2 = Bc(jj,:)*(r3dr(jj,:)-r3dr(jj-1,:))'/speed; % cosine of the angle between the velocity axis and the beam
        cos1s = cos1*cos1; 
%         cos2s = cos2*cos2;
        preerr1 = perror2d(jj-n_lag,1);        
        preerr2 = perror2d(jj-n_lag,2);    
        parameters = [cos1s preerr1 preerr2 speed];
    elseif jj > 1 + runlength
%          elseif jj > 1 + runlength && mod(jj,2) == 0
%         r3derr = [r3derr;rerror3d(jj,:)];
%         p3derr = [p3derr;perror3d(jj,:)];
%         p2derr = [p2derr;perror2d(jj,:)];
        threed = [threed;norm(rerror3d(jj,:))];
        twod = [twod;norm(perror3d(jj,:))];
        speed = norm(r3dr(jj-1,:)-r3dr(jj-3,:));
        cos1 = Bc(jj,:)*kr(jj-2,:)';
%         cos2 = kr(jj-1,:)*(r3dr(jj,:)-r3dr(jj-1,:))'/speed;
        cos1s = cos1*cos1;
%         cos2s = cos2*cos2;
        preerr1 = perror2d(jj-n_lag,1);        
        preerr2 = perror2d(jj-n_lag,2);   
        parameters = [parameters;[cos1s preerr1 preerr2 speed]];        
    end
    jj=jj+1;
end % of "while jj< = point_period."
% calculate the variation of the trajectory itself.
od=zeros(points_period,3);
% construct a new array for calculating the standard deviation
% of the trajectory.
od(:,:)=so3d(1:points_period,:);    
% sum up the square of all standard deviation.
variancex=std(od(1:points_period,1))*std(od(1:points_period,1))+...
    std(od(1:points_period,2))*std(od(1:points_period,2))+...
    std(od(1:points_period,3))*std(od(1:points_period,3));
std_dev2=variancex; 
aberror2 = aberror2/(points_period-n_lag);
ab_dev2 = ab_dev2/(points_period-n_lag);


    
        
    