
Npoints = [10,20,40,80,160];
num = size(Npoints,2);
n = size(Npoints,2);
optimal_cost = 33565;

avg_exp_benders_list = zeros(1,n);
avg_exp_saa_list =  zeros(1,n);
avg_exp_benders_list_neg10 = zeros(1,n);
avg_exp_benders_list_pos90 = zeros(1,n);
avg_exp_saa_list_neg10 = zeros(1,n);
avg_exp_saa_list_pos90 = zeros(1,n);
avg_T_benders_list = zeros(1,n);
avg_T_saa_list = zeros(1,n);
avg_T_max_sub_list = zeros(1,n);
avg_num_iter_list = zeros(1,n);

i = 0;
while i < n 
    i = i+1;
    N = Npoints(i)
    K = N;
    [avg_exp_benders,avg_exp_saa,benders10,benders90,saa10,saa90,avg_T_benders,avg_T_saa,avg_T_max_sub, avg_num_iter] = outperformance_CV_FLP(N,K);

    
    avg_exp_benders_list(i) = avg_exp_benders;
    avg_exp_saa_list(i) = avg_exp_saa;
    avg_exp_benders_list_neg10(i) = benders10;
    avg_exp_benders_list_pos90(i) = benders90;
    avg_exp_saa_list_neg10(i) = saa10;
    avg_exp_saa_list_pos90(i) = saa90;
    
    avg_T_benders_list(i) = avg_T_benders;
    avg_T_saa_list(i) = avg_T_saa;
    
    avg_T_max_sub_list(i) = avg_T_max_sub;
    avg_num_iter_list(i) = avg_num_iter;   
end

% plot of out of sample performance
font_size = 24;
alpha = 0.2;
color = [0.4940, 0.1840, 0.5560;
         0.4660, 0.6740, 0.1880];

fig1 = figure(1);
hold on

plot_benders = plot_with_shade(Npoints(1:num), avg_exp_benders_list(1:num), avg_exp_benders_list_neg10(1:num),avg_exp_benders_list_pos90(1:num), alpha, color(1,:),'o');
plot_saa = plot_with_shade(Npoints(1:num), avg_exp_saa_list(1:num), avg_exp_saa_list_neg10(1:num),avg_exp_saa_list_pos90(1:num), alpha, color(2,:),'|');
plot_optimal = plot(Npoints(1:num), ones(1,num)*optimal_cost,'linewidth', 3, 'color','r','marker','.','markersize',12);


grid on
set(gca, 'FontSize', font_size - 6);
xlabel('Number of In-sample Points N', 'Interpreter', 'latex', 'FontSize', font_size);
ylabel('Out of Sample Cost', 'Interpreter', 'latex', 'FontSize', font_size);
lgd = legend([plot_benders plot_saa plot_optimal], 'Benders C0','SAA','optimal','Location', 'northeast');
set(lgd,'Interpreter','latex', 'FontSize', font_size-6);

% plot of running time
fig2 = figure(2);
hold on
plot_bender_T = plot(Npoints(1:num),avg_T_benders_list(1:num),'linewidth', 3,'color', color(1,:)); 
plot_saa_T = plot(Npoints(1:num),avg_T_saa_list(1:num),'linewidth', 3,'color', color(2,:)); 

grid on
set(gca, 'FontSize', font_size - 6);
xlabel('Number of In-sample Points', 'Interpreter', 'latex', 'FontSize', font_size);
ylabel('Average Runtime(sec)', 'Interpreter', 'latex', 'FontSize', font_size);
lgd = legend([plot_bender_T plot_saa_T], 'Benders C0','SAA', 'Location', 'northeast');
set(lgd,'Interpreter','latex', 'FontSize', font_size-6);