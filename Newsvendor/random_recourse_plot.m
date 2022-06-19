seed = 42;
rng(seed);

X = [3.96;6.25];


size1 = 5;
size2 = 2;
c = 0;
d = 5;
data = c+d*rand(size1,size2);
data = data';


% set color by distance order
N = size(data,2);
dis = zeros(N,1);
for i = 1:N
    dis(i) = norm(data(:,i)-[0;0]);
end
[out,idx] = sort(dis);

color = cell(N,1);
cn = 1;
step = 0.8/N;
for i = 1:N
    color{idx(i)} = [0,cn,0];
    cn = cn -step;
end    


figure()
[VX,VY] = voronoi(data(1,:),data(2,:));
h = plot(VX,VY,'-k',data(1,:),data(2,:),'.r');
rectangle('Position',[0,0,5.5,6]);
axis equal

xlim([0,5.5])
ylim([0,6])
[V,R] = voronoin(data');


% Assign labels to the points X.
nump = size(data',1);
plabels = arrayfun(@(n) {sprintf('X%d', n)}, (1:nump)');
hold on
Hpl = text(data(1,:), data(2,:)+0.2, plabels, 'color', 'r', ...
      'FontWeight', 'bold', 'HorizontalAlignment',...
      'center', 'BackgroundColor', 'none');
  
% Assign labels to the Voronoi vertices V.
% By convention the first vertex is at infinity.
numv = size(V,1);
vlabels = arrayfun(@(n) {sprintf('V%d', n)}, (2:numv)');
hold on
Hpl = text(V(2:end,1), V(2:end,2)+.2, vlabels, ...
      'FontWeight', 'bold', 'HorizontalAlignment',...
      'center', 'BackgroundColor', 'none');
hold off

VI = [];
i = 1;
hold on
for j=1:length(VX(1,:))
     m = (VY(2,j)-VY(1,j))/(VX(2,j)-VX(1,j));
     yEqn = @(x)m*(x-VX(1,j))+VY(1,j);
     for x=[0 5.5]
        % check for straddling of x
        if (VX(1,j)<=x && VX(2,j)>=x) || (VX(1,j)>=x && VX(2,j)<=x)
            y = yEqn(x);
            % if y is between 0 and 10, then point of intersection is on 
            % line x=0 or x=10
            if y>=0 && y<=6
                VI = [VI;[x,y]];
                i = i + 1;
                plot(x,y,'r*');
                % skip to next edge
                continue;
            end
        end
    end
    xEqn = @(y)(y-VY(1,j)+m*VX(1,j))/m;
    % for equations y=0 and y=10
    for y=[0 6]
     % check for straddling of y
        if (VY(1,j)<=y && VY(2,j)>=y) || (VY(1,j)>=y && VY(2,j)<=y)
            x = xEqn(y);
            % if x is between 0 and 10, then point of intersection is on 
            % line y=0 or y=10
            if x>=0 && x<=5.5
                VI = [VI;[x,y]];
                i = i + 1;
                plot(x,y,'g*');
                % skip to next edge
                continue;
            end
        end
    end
end
hold off


vxi = V';
Z = max(X-vxi,0);
vi = VI';
Z0 = max(X-vi,0);
vb = [0,0,5.5,5.5;0,6,6,0];
Zb = max(X-vb,0);


X = cell(size(data,2),1);
Y = cell(size(data,2),1);
Z1 = cell(size(data,2),1);
Z2 = cell(size(data,2),1);
for j = 1:size(data,2)
    x = [];
    y = [];
    z1 = [];
    z2 = [];
    for i = R{j}
        if i ~= 1
            if (V(i,1) >= 0) && (V(i,2) >= 0)
                x = [x;V(i,1)];
                y = [y;V(i,2)];
                z1 = [z1;Z(1,i)];
                z2 = [z2;Z(2,i)];
            end
        end
    end
    X{j} = x;
    Y{j} = y;
    Z1{j} = z1;
    Z2{j} = z2;
end

figure
hold on
p1 = [vb(1,1);vi(1,1);vi(1,4);X{1}];
p2 = [vb(2,1);vi(2,1);vi(2,4);Y{1}];
p3 = [Zb(2,1);Z0(2,1);Z0(2,4);Z2{1}];
p = cell(size(p1,1),1);
for i = 1:size(p1,1)
    p{i} = [p1(i),p2(i),p3(i)];
end
perm = perms(p);  
for i = 1:size(perm,1)
    p1 = zeros(size(perm,2),1);
    p2 = zeros(size(perm,2),1);
    p3 = zeros(size(perm,2),1);
    for j = 1:size(perm,2)
        point = perm{i,j};
        p1(j) = point(1);
        p2(j) = point(2);
        p3(j) = point(3);
    end
    h = fill3(p1,p2,p3,color{1});
    set(h,'EdgeColor','none')
end
% alpha(0.3);
grid on;
xlabel("\xi_1")
ylabel("\xi_2")
zlabel("x_2({\bf \xi})")

hold on
p1 = [X{2};vi(1,4);vi(1,2)];
p2 = [Y{2};vi(2,4);vi(2,2)];
p3 = [Z2{2};Z0(2,4);Z0(2,2)];
p = cell(size(p1,1),1);
for i = 1:size(p1,1)
    p{i} = [p1(i),p2(i),p3(i)];
end
perm = perms(p);  
for i = 1:size(perm,1)
    p1 = zeros(size(perm,2),1);
    p2 = zeros(size(perm,2),1);
    p3 = zeros(size(perm,2),1);
    for j = 1:size(perm,2)
        point = perm{i,j};
        p1(j) = point(1);
        p2(j) = point(2);
        p3(j) = point(3);
    end
    h = fill3(p1,p2,p3,color{2});
    set(h,'EdgeColor','none')
end

hold on
p1 = [X{3};vi(1,1);vi(1,3)];
p2 = [Y{3};vi(2,1);vi(2,3)];
p3 = [Z2{3};Z0(2,1);Z0(2,3)];
p = cell(size(p1,1),1);
for i = 1:size(p1,1)
    p{i} = [p1(i),p2(i),p3(i)];
end
perm = perms(p);  
for i = 1:size(perm,1)
    p1 = zeros(size(perm,2),1);
    p2 = zeros(size(perm,2),1);
    p3 = zeros(size(perm,2),1);
    for j = 1:size(perm,2)
        point = perm{i,j};
        p1(j) = point(1);
        p2(j) = point(2);
        p3(j) = point(3);
    end
    h = fill3(p1,p2,p3,color{3});
    set(h,'EdgeColor','none')
end

hold on
p1 = [X{4};vi(1,4);vi(1,3);vb(1,1);vb(1,4)];
p2 = [Y{4};vi(2,4);vi(2,3);vb(2,1);vb(2,4)];
p3 = [Z2{4};Z0(2,4);Z0(2,3);Zb(2,1);Zb(2,4)];
p = cell(size(p1,1),1);
for i = 1:size(p1,1)
    p{i} = [p1(i),p2(i),p3(i)];
end
perm = perms(p);  
for i = 1:size(perm,1)
    p1 = zeros(size(perm,2),1);
    p2 = zeros(size(perm,2),1);
    p3 = zeros(size(perm,2),1);
    for j = 1:size(perm,2)
        point = perm{i,j};
        p1(j) = point(1);
        p2(j) = point(2);
        p3(j) = point(3);
    end
    h = fill3(p1,p2,p3,color{4});
    set(h,'EdgeColor','none')
end

hold on
p1 = [X{5};vi(1,2);vi(1,5);vb(1,3)];
p2 = [Y{5};vi(2,2);vi(2,5);vb(2,3)];
p3 = [Z2{5};Z0(2,2);Z0(2,5);Zb(2,3)];
p = cell(size(p1,1),1);
for i = 1:size(p1,1)
    p{i} = [p1(i),p2(i),p3(i)];
end
perm = perms(p);  
for i = 1:size(perm,1)
    p1 = zeros(size(perm,2),1);
    p2 = zeros(size(perm,2),1);
    p3 = zeros(size(perm,2),1);
    for j = 1:size(perm,2)
        point = perm{i,j};
        p1(j) = point(1);
        p2(j) = point(2);
        p3(j) = point(3);
    end
    h = fill3(p1,p2,p3,color{5});
    set(h,'EdgeColor','none')
end
