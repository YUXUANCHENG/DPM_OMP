clc
clear
close all


%folder = "C:\Users\Yuxuan Cheng\OneDrive\Windows_Version\remove3\";
%folder = "C:\Users\Yuxuan Cheng\source\repos\distritution\Debug\remove6\";
folder = "C:\Users\Yuxuan Cheng\source\repos\distritution\Debug\negative4\";

matrix_eigen = {};
contact_number = [];

for t_index = 1:41
    
    extend = "_" + int2str(t_index) +".txt";
    contact_file = folder + "contact" + extend;
    contact = csvread(contact_file);
    mean_contact = mean(sum(contact,2),'all');
    
    contact_number = [contact_number, mean_contact];
end

scatter(-23:0.5:-3, flip(contact_number)/max(contact_number));
%scatter(-3:0.025:-2, flip(contact_number)/contact_number(1));
xlabel("log(T)");
ylabel("contact number");
