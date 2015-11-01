t = 0:pi/20:4*pi;
x = sin(t);
y = 2*sin(t+pi/2) + 1/3*sin(2*t);
z = 3*sin(t+0.5) + 1/2*sin(2*t+0.5);
figure(1)
% plot(x,y)
plot3(x,y,z)
y = 2*sin(t+pi) + 1/3*sin(2*t+0.5);
z = 3*sin(t) + 1/2*sin(2*t+2.5);
figure(2)
% plot(x,y)
plot3(x,y,z)
y = 2*sin(t) + 1/3*sin(2*t+pi/2);
z = 3*sin(t+0.5) + 1/2*sin(2*t+1.5);
figure(3)
% plot(x,y)
plot3(x,y,z)
y = 2*sin(t) + 1/3*sin(2*t+pi);
z = 3*sin(t+0.5) + 1/2*sin(2*t+0.5);
figure(4)
% plot(x,y)
plot3(x,y,z)
y = 2*sin(t+pi/2) + 1/3*sin(2*t+pi/2);
z = 3*sin(t) + 1/2*sin(2*t+2.5);
figure(5)
% plot(x,y)
plot3(x,y,z)
y = 2*sin(t+pi/2) + 1/3*sin(2*t+pi);
z = 3*sin(t+0.5) + 1/2*sin(2*t+1.5);
figure(6)
% plot(x,y)
plot3(x,y,z)