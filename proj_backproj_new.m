function [r3dx,centerx,kx,ko,ym,linearity,isplane,Bc,cost]=proj_backproj_new(seg_err,...
    marginx,so3d,r3d,jx,gAng,runlength,sid,stdx)
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

%% Initialize the parameters. 
centerx=[0 0 0];
kx=[0 0 0];
ko=[0 0 0];
ym=[0 0];
isplane=0;
linearity=1;
xfit=zeros(jx,3);
r3d0=zeros(jx,3);

% it has been compared to 1:jx-1. The following is better.
for ix=1:runlength     
    r3d0(ix,1)=r3d(ix,2);
    r3d0(ix,2)=r3d(ix,1);
    r3d0(ix,3)=stdx-r3d(ix,3);
end

p2d=[0 0 0];
wc3=[0 0 0];
prev = [0 0 0];
% if the point data in the database is continuous(based on the time 
% information), perform estimation, flag it otherwise.
o_3d=[0 0 0];
o_3d(:)=so3d(jx,:);

%% project forward
p2d(1:3)=proj([o_3d(1) o_3d(2) o_3d(3)],gAng,sid,stdx);
[coeff,score,roots] = princomp(r3d(1:jx-1,1:3)); % only earlier points can be referred to in the calculation
centerx = mean(r3d(1:jx-1,1:3),1);
xfit=repmat(centerx,jx-1,1) + score(:,1)*coeff(:,1)';
linearity = roots(1)/sum(roots);
if linearity>0.001   % a parameter can be adjusted, it is no use though.
    isplane=0;
else
    isplane=1;
end
dir_vector = coeff(:,1)';
% diro_vector = coeffo(:,1)';
kx(:)=dir_vector(:);
% ko(:)=diro_vector(:);
%% add segmentation error to the data
p2d(1)=p2d(1) + seg_err * randn();
p2d(2)=p2d(2) + seg_err * randn();
ymin=min(xfit(:,1));
ymax=max(xfit(:,1));  
ym=[ymin ymax];
Bc=[0 0 0];
marginx = marginx*(ymax-ymin);
prev = r3d(jx-1,:) - r3d(jx-2,:);
prep = r3d(jx-1,:);
%% backproject according to whether it's a plane or not.
% if isplane==0
    [wc3 B cost]=back_proj_line(centerx,kx,prev,prep,marginx,ymin,ymax,stdx,sid,gAng,p2d(1),p2d(2));
    r3dx(1)= wc3(2);
    r3dx(2)= wc3(1);
    r3dx(3)= -wc3(3)+stdx;
    Bc(1) = B(2);
    Bc(2) = B(1);
    Bc(3) = -B(3);
    Bc = Bc/norm(Bc);
% else
%     wc3=back_proj_plane(centerx,kx,r3d0,stdx,sid,gAng,p2d(1),p2d(2));
%     r3dx(1)= wc3(2);
%     r3dx(2)= wc3(1);
%     r3dx(3)= -wc3(3)+stdx;
% end
% xx=range(so3d(1:jx,1:3),1);
% range_o3d=xx(1,:); % ??
return