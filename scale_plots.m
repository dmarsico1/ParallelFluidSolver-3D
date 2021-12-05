proc_vec=[2 4 6 8 10 12 16 20];
time_vec=[62.43 31.8 21.75 15.72 12.71 11.03 8.63 6.97];
time_vec_ideal=[62.43 63.43/2 62.43/4 62.43/5 62.43/8 62.43/10];
plot(log(proc_vec),log(time_vec),'r*');
legend('Actual','Ideal')
xlabel('Log(Number of Processors)')
ylabel('Log(CPU time)')
