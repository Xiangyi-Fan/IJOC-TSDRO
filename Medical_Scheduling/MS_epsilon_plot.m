epsilons = [0.1,0.5,1.5,3,6,10];
Npoints = 10;
K = Npoints;

num = size(epsilons,2);
n = size(epsilons,2);  
avg_cvar_cop_list = zeros(1,n);
avg_cvar_sp_list = zeros(1,n);
avg_cvar_saa_list =  zeros(1,n);
avg_cvar_cop_list_neg10 = zeros(1,n);
avg_cvar_cop_list_pos90 = zeros(1,n);
avg_cvar_sp_list_neg10 = zeros(1,n);
avg_cvar_sp_list_pos90 = zeros(1,n);
avg_cvar_saa_list_neg10 = zeros(1,n);
avg_cvar_saa_list_pos90 = zeros(1,n);

for i = 1:n
    epsilon = epsilons(i)*K
    [avg_cvar_cop,avg_cvar_sp,avg_cvar_saa,cop10,cop90,sp10,sp90,saa10,saa90] = outperformance_MS(epsilon,Npoints,K);
    avg_cvar_cop_list(i) = avg_cvar_cop;
    avg_cvar_sp_list(i) = avg_cvar_sp;
    avg_cvar_saa_list(i) = avg_cvar_saa;
    avg_cvar_cop_list_neg10(i) = cop10;
    avg_cvar_cop_list_pos90(i) = cop90;
    avg_cvar_sp_list_neg10(i) = sp10;
    avg_cvar_sp_list_pos90(i) = sp90;
    avg_cvar_saa_list_neg10(i) = saa10;
    avg_cvar_saa_list_pos90(i) = saa90;  
end


% plot of out of sample performance
font_size = 24;
alpha = 0.2;
color = [0.9290, 0.6940, 0.1250;
         0.4940, 0.1840, 0.5560;
         0.4660, 0.6740, 0.1880];

fig1 = figure(1);
hold on

plot_cop = plot_with_shade(epsilons(1:num), avg_cvar_cop_list(1:num), avg_cvar_cop_list_neg10(1:num),avg_cvar_cop_list_pos90(1:num), alpha, color(1,:),'^');
plot_sp = plot_with_shade(epsilons(1:num), avg_cvar_sp_list(1:num), avg_cvar_sp_list_neg10(1:num),avg_cvar_sp_list_pos90(1:num), alpha, color(2,:),'o');
plot_saa = plot_with_shade(epsilons(1:num), avg_cvar_saa_list(1:num), avg_cvar_saa_list_neg10(1:num),avg_cvar_saa_list_pos90(1:num), alpha, color(3,:),'|');

grid on
set(gca, 'FontSize', font_size - 6);
xlabel('Robustness Parameter $\epsilon_k$', 'Interpreter', 'latex', 'FontSize', font_size);
ylabel('Out of Sample Cost', 'Interpreter', 'latex', 'FontSize', font_size);
lgd = legend([plot_cop plot_sp plot_saa], 'C1','C0','SAA', 'Location', 'northeast');
set(lgd,'Interpreter','latex', 'FontSize', font_size-6);

