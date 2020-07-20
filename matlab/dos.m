clc
clear
close all



%folder = "C:\Users\Yuxuan Cheng\source\repos\distritution\Debug\decompress14\";
%folder = "C:\Users\Yuxuan Cheng\source\repos\distritution\Debug\64\";
folder = "C:\Users\Yuxuan Cheng\source\repos\distritution\Debug\162\";

matrix_eigen = [];
phi_e = [];
contact_r = [];

for t_index = 1:200
    try
    extend = "_" + int2str(t_index) +".txt";
    eigen_file = folder + "dynamical_file_remove" + extend;
    eigen = csvread(eigen_file);
    phi_file = folder + "final_phi" + extend;
    phi_v = csvread(phi_file);
    contact_file = folder + "contact" + extend;
    contact = csvread(contact_file);
    
    catch
        continue
    end
    
    if (phi_v<0.38)
        continue
    end
    matrix_eigen = [matrix_eigen, eigen'];
    phi_e = [phi_v, phi_e];
    contact_r = [sum(contact==1,2)/sum(contact>0,2), contact_r];
    %contact_r = [sum(contact==1,2), contact_r];
end

matrix_eigen(matrix_eigen< 10^-7) = 0;
matrix_eigen = matrix_eigen(matrix_eigen > 0);
figure(1);
[counts, binCenters] = hist(sqrt(matrix_eigen),30);
counts = counts / sum(counts);
plot(binCenters, counts,'linewidth',1.75);
xlabel("w");
ylabel("D(w)");

figure(2);hold on;
scatter(phi_e,contact_r);
P = polyfit(phi_e,contact_r,1);
yfit = P(1)*phi_e+P(2);
plot(phi_e,yfit);
xlabel("phi");
ylabel("repulsive contact");
