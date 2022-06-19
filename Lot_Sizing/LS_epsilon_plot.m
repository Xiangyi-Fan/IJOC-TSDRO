
epsilons = [5,10,20];
n = size(epsilons,2);  
avg_dro_list = zeros(1,n);
avg_saa_list =  zeros(1,n);
avg_dro_list_neg20 = zeros(1,n);
avg_dro_list_pos80 = zeros(1,n);
avg_saa_list_neg20 = zeros(1,n);
avg_saa_list_pos80 = zeros(1,n);

for i = 1:n
    epsilon = epsilons(i)
    [avg_dro,avg_saa,dro20,dro80,saa20,saa80] = outperformance_LS(epsilon);
    avg_dro_list(i) = avg_dro;
    avg_saa_list(i) = avg_saa;
    avg_dro_list_neg20(i) = dro20;
    avg_dro_list_pos80(i) = dro80;
    avg_saa_list_neg20(i) = saa20;
    avg_saa_list_pos80(i) = saa80;  
end

figure()
xlabel("epsilon");
ylabel("Out of Sample Cost)($)");
hold on
plot_dro=plot(epsilons,avg_dro_list,'b'); 
plot_saa=plot(epsilons,avg_saa_list,'r'); 
legend([plot_dro plot_saa],{'DRO','SAA'});    