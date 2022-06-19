Npoints = [10,20,40,80,160];
num = size(Npoints,2);
n = size(Npoints,2); 
optimal_cost = 10697;

avg_cop_list = zeros(1,n);
avg_sp_list = zeros(1,n);
avg_saa_list =  zeros(1,n);
avg_wass_list =  zeros(1,n);
avg_cop_list_neg0 = zeros(1,n);
avg_cop_list_pos100 = zeros(1,n);
avg_sp_list_neg0 = zeros(1,n);
avg_sp_list_pos100 = zeros(1,n);
avg_saa_list_neg0 = zeros(1,n);
avg_saa_list_pos100 = zeros(1,n);
avg_wass_list_neg0 = zeros(1,n);
avg_wass_list_pos100 = zeros(1,n);
avg_T_cop_list = zeros(1,n);
avg_T_sp_list = zeros(1,n);
avg_T_saa_list = zeros(1,n);
avg_T_wass_list = zeros(1,n);
feasibility_num_list = zeros(1,n);

seed = 1;
rng(seed);

trial = 20;
trial_new = trial*1;
N = 5; % number of stocks
data0_in = cell(trial_new,1);
data_cons0 = LS_generate_data_uniform(20);
for i = 1:trial_new
    data0_in{i,1} = LS_generate_data(20);
end
disp("finish generating data")


epsilon_cop_list = zeros(trial_new,1);
epsilon_sp_list = zeros(trial_new,1);
i = 0;
while i < trial_new
    i = i+1;
    epsilon_cop_list(i) = LS_epsilon_cop_value_general(data_cons0,data0_in{i,1});
    epsilon_sp_list(i) = LS_epsilon_sp_value_general(data_cons0,data0_in{i,1});
end
epsilon_cop_fix = mode(epsilon_cop_list)
epsilon_sp_fix = mode(epsilon_sp_list)

% cross validation for gamma
epsilon_cop_p_list = zeros(trial_new,1);
epsilon_sp_p_list = zeros(trial_new,1);


for i = 1:trial_new
    epsilon_cop_p_list(i) = LS_epsilon_cop_p_value_general(data_cons0,data0_in{i,1},epsilon_cop_fix);
    epsilon_sp_p_list(i) = LS_epsilon_cop_p_value_general(data_cons0,data0_in{i,1},epsilon_sp_fix);
end
epsilon_cop_p_fix = mode(epsilon_cop_p_list)
epsilon_sp_p_fix = mode(epsilon_sp_p_list)


i = 0;
% cont = true;
while i < n
    i = i+1;
    N = Npoints(i)
    K = N;
    [avg_cop,avg_sp,avg_saa,avg_wass,cop0,cop100,sp0,sp100,saa0,saa100,wass0,wass100,feasibility_num,avg_T_cop,avg_T_sp,avg_T_saa,avg_T_wass] = outperformance_LS_general(N,epsilon_cop_fix,epsilon_sp_fix,epsilon_cop_p_fix,epsilon_sp_p_fix);
    
    avg_cop_list(i) = avg_cop;
    avg_sp_list(i) = avg_sp;
    avg_saa_list(i) = avg_saa;
    avg_wass_list(i) = avg_wass;
    avg_cop_list_neg0(i) = cop0;
    avg_cop_list_pos100(i) = cop100;
    avg_sp_list_neg0(i) = sp0;
    avg_sp_list_pos100(i) = sp100;
    avg_saa_list_neg0(i) = saa0;
    avg_saa_list_pos100(i) = saa100;
    avg_wass_list_neg0(i) = wass0;
    avg_wass_list_pos100(i) = wass100; 
    
    avg_T_cop_list(i) = avg_T_cop;
    avg_T_sp_list(i) = avg_T_sp;
    avg_T_saa_list(i) = avg_T_saa;
    avg_T_wass_list(i) = avg_T_wass;
    
    feasibility_num_list(i) = feasibility_num;
end

% plot of out of sample performance
font_size = 24;
alpha = 0.2;
color = [0.9290, 0.6940, 0.1250;
         0.4940, 0.1840, 0.5560;
         0.4660, 0.6740, 0.1880;
         0.3010, 0.7450, 0.9330];

fig1 = figure(1);
hold on

plot_cop = plot_with_shade(Npoints(1:num), avg_cop_list(1:num), avg_cop_list_neg0(1:num),avg_cop_list_pos100(1:num), alpha, color(1,:),'^');
plot_sp = plot_with_shade(Npoints(1:num), avg_sp_list(1:num), avg_sp_list_neg0(1:num),avg_sp_list_pos100(1:num), alpha, color(2,:),'*');
plot_saa = plot_with_shade(Npoints(1:num), avg_saa_list(1:num), avg_saa_list_neg0(1:num),avg_saa_list_pos100(1:num), alpha, color(4,:),'|');
plot_optimal = plot(Npoints(1:num), ones(1,num)*optimal_cost,'linewidth', 3,'color', 'r','marker','.','markersize',12);

grid on
set(gca, 'FontSize', font_size - 6);
xlabel('Number of In-sample Points', 'Interpreter', 'latex', 'FontSize', font_size);
ylabel('Out of Sample Cost', 'Interpreter', 'latex', 'FontSize', font_size);
lgd = legend([plot_cop plot_sp plot_saa], 'C1','C0','Wasserstein','SAA','Location', 'northeast');
set(lgd,'Interpreter','latex', 'FontSize', font_size-6);
