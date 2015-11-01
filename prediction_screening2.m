function [r2dp predictionquality]=prediction_screening2(jx,n_lag,r2d,incr_2d,extentv,factor3)
% to screen out extraordinary predictions
% input: jx-the no. of current point; n_lag- the number of lag points; r2d:
%        the 3D coordinates on which the prediction is based; extentx:
%        reference range of variation; factor3: a factor to adjust the
%        ratio of the predicted change to the reference range. incr_2d: the
%        predicted change prior to the screening.
% output: the predicted 3D coordinates.

% SI motion
r2dp=[0 0];
predictionquality = 1;
for i2 = 1:2
    maxi(i2)=extentv(i2)*factor3;
    if incr_2d(jx-n_lag,i2)>maxi(i2)
        r2dp(i2)=r2d(jx-n_lag,i2)+maxi(i2);
        if r2d(jx-n_lag,i2)-r2d(jx-2*n_lag,i2)<-extentv(i2)*0.6 %If the previous velocity is small, then threshold acceleration is exceeded
            predictionquality = 0.3;
        end
    elseif incr_2d(jx-n_lag,i2)<-maxi(i2);
        r2dp(i2)=r2d(jx-n_lag,i2)-maxi(i2);
        if r2d(jx-n_lag,i2)-r2d(jx-2*n_lag,i2)>extentv(i2)*0.6 %If the previous velocity is large, then threshold acceleration is exceeded
            predictionquality = 0.3;
        end
    else
        r2dp(i2)=r2d(jx-n_lag,i2)+incr_2d(jx-n_lag,i2);
        if abs(incr_2d(jx-n_lag,i2) + r2d(jx-n_lag,i2)-r2d(jx-2*n_lag,i2)) > extentv(i2)*factor3*factor3           
            predictionquality = 0.3;
        end
    end
end