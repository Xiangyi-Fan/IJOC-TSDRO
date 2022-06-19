function [epsilon_k] =  LS_PLD_epsilon_thm(rho1, data, cons_data, k)
[P_k] = LS_PLD_partitions(cons_data,k);
[R_k] =  LS_PLD_radius(P_k, cons_data, k);
[I_k] =  LS_PLD_nums(P_k, data);
epsilon_k = R_k*R_k/sqrt(I_k)*(2+sqrt(2*log(1/rho1)));

end