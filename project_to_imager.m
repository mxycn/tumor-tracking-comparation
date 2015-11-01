function [o2d r2dx]=project_to_imager(seg_err,so3d,jx,gAng,sid,stdx)
% this function projects the actual points to the MV imager and 
% reconstruct them

% r3dx returns the new reconstruction result.
% is_discon signals the discontinuity in the original data.
% is_jump signals the large velocity in the actual 3D motion.

% seg_err is the detection error(=0.5mm)
% k is the no. of fraction.
% r3d is the estimated 3D positions
% jx is the no. of point to be estimated.
% delta_angle is the angle of gantry rotation during the first 4 seconds.

%%
p2d=[0 0];
% if the point data in the database is continuous(based on the time 
% information), perform estimation, flag it otherwise.
o_3d=[0 0 0];
o_3d(:)=so3d(jx,:);

%% project forward
p2d(1:3)=proj([o_3d(1) o_3d(2) o_3d(3)],gAng,sid,stdx);
o2d(1:2)=p2d(1:2);
p2d(1)=p2d(1) + seg_err * randn();
p2d(2)=p2d(2) + seg_err * randn();
r2dx(1:2) = p2d(1:2);
return