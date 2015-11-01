function [incrx predictionquality]=predicting2(is_pre,factor1,factor2,r2threshold,num_train,n_lag,jx,extentx,extentv,av3d)
% This function returns the predicted changes in the three directions of 
% motion, incr1, incr2, and incr3, which will be added to the coordinates 
% n_lag points earlier in function workflow(...).

% is_pre = 1 if prediction is applied, 0 otherwise.

% factor1 specifies the range of coordinate within which the points are
% selected for multivariate regression. factor1=selected range/full extent
% in that direction in the first 4 seconds.

% factor2 specifies the range of velocity within which the points are
% selected for multivariate regression. factor2=selected range/full extent
% in that direction in the first 4 seconds.

% n_lag is the number of latent points, n=3 corresponds to 460ms, n=2
% cooresponds to 310ms.

% extentx gives the range of motion in the first 4 seconds.

% The following are coefficients of regression equation, for the SI, LR, AP
% directions. Two-variate regression is better than three-variate
% regression, so only three coefficient for each direction of motion is
% neccessary.

% the regression equation is
%     incrx(:)=a+b*av3d(jx-1,:)+c*(av3d(jx-3,:);
        
incrx = [0 0];
predictionquality = 1;
if is_pre == 1
    delta_y_theta_k =[0 0];
    delta_y_theta_k_l = [0 0];
    delta_y_theta_k_2l = [0 0];
    delta_y_theta_k(:)= av3d(jx-n_lag,:)-av3d(jx,:);
    delta_y_theta_k_l(:)=av3d(jx-n_lag-1,:)-av3d(jx-1,:);
    delta_y_theta_k_2l(:)=av3d(jx-1-n_lag-1,:)-av3d(jx-2,:);
    
    % The variables entering the regression, selected from the raw data,
    % the last being the dependent variable.
    v3d1=zeros(jx,2);
    v3d2=zeros(jx,2);
    v3d=zeros(jx,2);
    % The indices for different points in the regression.
    k=[0 0];    
    % Construct the independent variables and the dependent variable for the 
    % regression in the SI direction.
    % Consider only those points that are similar to the present point in terms
    % of position and velocity.
    % ii is the index for the input data, which is screened by the conditions.
    % Because the regression in a sense is an autoregression, and the input data
    % involves av3d(ii-3,1), and av3d(ii+n_lag,1), ii runs from 4 to jx-1-n_lag-3,
    % where -3 is added to eliminate the points close to the current point.
    
    % See file "details" for explanation, here av3d(ii-1,1)+av3d(ii-2,1) replaces
    % av3d(ii-1,1)
    for i2 = 1:2
        mark_vm=[0 0];
        for ii=4:jx-1-n_lag
            % the reference for the latest velocity        
            mark_vm(i2)=av3d(jx-1,i2)-av3d(jx-2,i2);  
            % av3d(ii-1,1) plays the same role in the regression as what av3d(jx-1,1)
            % plays in the output. The first condition compares the position, the 
            % second compares the velocity
            if abs(av3d(ii-1,i2)-av3d(jx-1,i2))<extentx(i2)/factor1 && ...
                    abs(av3d(ii-1,i2)-av3d(ii-2,i2)-mark_vm(i2))<extentv(i2)*factor2 
                k(i2)=k(i2)+1;
                v3d1(k(i2),i2)=av3d(ii-1,i2); % the first independent variable.
                v3d2(k(i2),i2)=av3d(ii-3,i2); % the second independent variable.
                v3d(k(i2),i2)=av3d(ii+n_lag,i2)-av3d(ii,i2);   % the dependent variable.
            end
        end    
        % two-variate regression
        if k(i2)>num_train % if more than 3 set of numbers are acquired as the input data, do regression            
            [b r2]=two_var_regress(k(i2),v3d1(1:k(i2),i2),v3d2(1:k(i2),i2),v3d(1:k(i2),i2));  
            if r2<predictionquality
                predictionquality = r2;
            end
            if r2>r2threshold
                incrx(i2)=b(1)+b(2)*av3d(jx-1,i2)+b(3)*av3d(jx-3,i2);                   
            else
                 incrx(i2)=mean(v3d(1:k(i2),i2));                  
            end            
        elseif k(i2)>0.2 % if only one, two or three sets of numbers are acuquired, use the average
            incrx(i2)=mean(v3d(1:k(i2),i2));  
            predictionquality = 0;
%             if no reference point, predict according to the shape(very rough, MULIN algorithm)
        else
            incrx(i2)=0;
            predictionquality = 0;
        end
    end
end