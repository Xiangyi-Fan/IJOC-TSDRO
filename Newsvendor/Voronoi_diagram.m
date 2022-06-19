seed = 3;
rng(seed);

Npoints = 10;
data1 = generate_data_fix(5);
data2 = generate_data_fix(5);
data3 = generate_data_fix(15);
data4 = generate_data_fix(25);
data5 = generate_data_fix(50);
data6 = generate_data_fix(9900);
data = [data1,data2,data3,data4,data5,data6];
 
voronoi(data(1,:),data(2,:))