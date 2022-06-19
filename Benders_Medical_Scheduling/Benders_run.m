Npoints_list = [10,20,40,80,160];
N = size(Npoints_list,2);
runtime = zeros(N,1);
epsilon0 = 20;

for i = 1:N
    Npoints = Npoints_list(i);
    K = Npoints;
    data = MS_generate_data(Npoints);
    optimal = Benders(data,epsilon0,K);
    runtime(i) = optimal.runtime;
    text = ['finish Benders with',num2str(Npoints),'in_sample data'];
    disp(text)
end

runtime