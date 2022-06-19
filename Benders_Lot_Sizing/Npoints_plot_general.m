
Npoints = [160];
num = size(Npoints,2);
n = size(Npoints,2);


avg_T_benders_list = zeros(1,n);
avg_T_max_sub_list = zeros(1,n);
avg_num_iter_list = zeros(1,n);

i = 0;
while i < n 
    i = i+1;
    N = Npoints(i);
    K = N;
    [avg_T_benders, avg_T_max_sub, avg_num_iter] = outperformance_LS_general(N,K);    
    avg_T_benders_list(i) = avg_T_benders;  
    avg_T_max_sub_list(i) = avg_T_max_sub;
    avg_num_iter_list(i) = avg_num_iter;
end

% plot of out of sample performance
font_size = 28;
alpha = 0.2;
color = [0.4940, 0.1840, 0.5560;
         0.4660, 0.6740, 0.1880];


% plot of running time
fig1 = figure(1);
hold on
plot_bender_T = plot(Npoints(1:num),avg_T_benders_list(1:num),'linewidth', 3,'color', color(1,:)); 

grid on
set(gca, 'FontSize', font_size - 6);
xlabel('Number of In-sample Points', 'Interpreter', 'latex', 'FontSize', font_size);
ylabel('Average Runtime(sec)', 'Interpreter', 'latex', 'FontSize', font_size);
lgd = legend(plot_bender_T, 'Benders C0','Location', 'northeast');
set(lgd,'Interpreter','latex', 'FontSize', font_size-6);
