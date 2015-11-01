function wc3=back_proj_line(centerpoint,kp,r3dx,r,f,gAng2,u,v)
%
% % see Appendices of Liuwu for details
% Programmed by: Huagang Yan
% 07/15/2010
%% Input
% k: the slope of the tracks in x-y plane in the world coordinates,
% obtained in the first 4 seconds by regression based on the triangulation
% between kV data and MV data.
% r: the source to axis distance(SAD)
% f: the source to imager plate distance
% gAng1: the gantry angle of the last projection
% gAng2: the gantry angle of the present projection
% u1,u2: the pixel coorinate x on the portal imager 
%        left + (data coordinate y, liuwu coordinate x)
% v1,v2: the pixel coorinate z on the portal imager  
%        infer + (data coordinate x, liuwu coordiante y)
%   Anterior+: data coordinate z+, Liuwu coordinate z-
%--------------------------------------------------------------------------

O=[-r*sin(gAng2*pi/180.0),0,r*(1-cos(gAng2*pi/180.0))]';
Ap=[u,v,f]';
A=[0 0 0]';
R=[cos(gAng2*pi/180.0) 0 sin(gAng2*pi/180.0); 0 1 0; -sin(gAng2*pi/180.0) 0 cos(gAng2*pi/180.0)];

A=R*Ap+O; % the 3D coordinate of the current projection point

B=A-O;
is_useyz=0;
P=[0 0 0]';
wc3=[0 0 0]';

% in the following, we are going to seek a point on the current projection
% line that is of the shortest distance extended from the the last
% projection line in the plane that satisfies y= k x.-- according to liuwu's convention
% intermediate step1
% in liuwu's coordinate, x=k1(1)*y+k1(2);z=k2(1)*y+k2(2);z=k3(1)*x+k3(2);
% and the new form is:
% the plane equation is: (x-centerpoint(2))*b+(y-centerpoint(1))*a+(z-std+centerpoint(3)*(-c)=0;

cos_theta=kp*B/norm(B)/norm(kp);
is_tooparallel=0;
if cos_theta*cos_theta<0.8
    is_tooparallel=1;
end

center_p=[centerpoint(2) centerpoint(1) r-centerpoint(3)]';
norm_p=[kp(2) kp(1) -kp(3)]';
ref_l=O;
direc_l=B;
prepoint=[r3dx(size(r3dx,1),2) r3dx(size(r3dx,1),1) r-r3dx(size(r3dx,1),3)];

point_intersection=inter_line_plane(center_p,norm_p,ref_l,direc_l);

% if is_tooparallel==0
%     P=point_intersection;
% else
    [point_nearest nearestpoint]=find_nearest(r3dx,ref_l,direc_l);
    % find a point in the line that is closest to the the immediate
    % previous point.
    point_adjacent=find_closest(prepoint,ref_l,direc_l); 
% adjust whether the intersection point is far beyond the previous points
% define a farther point to the center
%     nearestpoint(1)=(nearestpoint(1)-centerpoint(2))*1.1+centerpoint(2);
%     nearestpoint(2)=(nearestpoint(2)-centerpoint(1))*1.1+centerpoint(1);
%     nearestpoint(3)=(nearestpoint(3)-r+centerpoint(3))*1.1+r-centerpoint(3);
%     cosbeta=(centerpoint(2)-nearestpoint(1))*(point_intersection(1)-nearestpoint(1))+...
%         (centerpoint(1)-nearestpoint(2))*(point_intersection(2)-nearestpoint(2))+...
%         (r-centerpoint(3)-nearestpoint(3))*(point_intersection(3)-nearestpoint(3));
%     if cosbeta<0 % yes, it is far beyond the range, meaning the error is too large.
%        P=point_nearest;
        P=point_nearest;
%        P=point_adjacent;
%        P=(point_nearest+point_adjacent)/2;
%     else
%        P=(point_nearest+point_adjacent)/2*(1-abs(cos_theta))+abs(cos_theta)*point_intersection;
%     end
% end

wc3(1:3)=P;
return;