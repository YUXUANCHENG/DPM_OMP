clc
clear
close all


folder = "C:\Users\Yuxuan Cheng\source\repos\distritution\Debug\negative1\";

file = folder + "eigen_energy.txt";
energy_list = csvread(file);


scatter(log10(1:100),log10(energy_list))
xlabel("log(r)");
ylabel("log(E)");