pointA = [0,0,0];
pointB = [-10,-20,10];
pointC = [10,20,10];
points=[pointA' pointB' pointC']; % using the data given in the question
fill3(points(1,:),points(2,:),points(3,:),0.5,'LineStyle','none');
grid on
alpha(0.3) 