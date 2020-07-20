clc
clear
close all


folder = "C:\Users\Yuxuan Cheng\source\repos\distritution\Debug\negative8\";

matrix_eigen = {};
contact_number = [];
contact_number_a = [];
contact_number_r = [];

for t_index = 1:41
    
    extend = "_" + int2str(t_index) +".txt";
    contact_file = folder + "contact" + extend;
    contact = csvread(contact_file);
    mean_contact = mean(sum(contact>0,2),'all');
    mean_contact_r = mean(sum(contact==1,2),'all');
    mean_contact_a = mean(sum(contact==2,2),'all');   
    contact_number = [contact_number, mean_contact];
    contact_number_r = [contact_number_r, mean_contact_r];
    contact_number_a = [contact_number_a, mean_contact_a];
end

figure(1); hold on
scatter(-23:0.5:-3, flip(contact_number)/max(contact_number));
scatter(-23:0.5:-3, flip(contact_number_a)/max(contact_number));
scatter(-23:0.5:-3, flip(contact_number_r)/max(contact_number));
%scatter(-3:0.025:-2, flip(contact_number)/contact_number(1));
xlabel("log(T)");
ylabel("contact number");
legend("total contact number", "attractive region", "repulsive region")
















