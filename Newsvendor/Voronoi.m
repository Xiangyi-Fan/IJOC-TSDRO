clc
clear all
n= 10;
x=10*rand(1,n);
y=10*rand(1,n);
h=voronoi(x,y);
[vx,vy] =voronoi(x,y);
[v,c] = voronoin([x(:) y(:)]);
close all
plot(x,y,'r+',vx,vy,'b-'); 
rectangle('Position',[0,0,12,12]);
axis equal
figure(2)
for j=1:length(vx(1,:))
    line([vx(1,j) vx(2,j)],[vy(1,j) vy(2,j)])
    rectangle('Position',[0,0,12,12]);
end
axis equal