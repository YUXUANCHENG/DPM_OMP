clc
clear
close all


%folder = "C:\Users\Yuxuan Cheng\source\repos\distritution\Debug\decompress16\";
folder = "C:\Users\Yuxuan Cheng\source\repos\distritution\Debug\over3\";
t_index = 1;

extend = "_" + int2str(t_index) +".txt";
contact_file = folder + "contact" + extend;
%contact = csvread(contact_file);
fre_file = folder + "dynamical_file" + extend;
fre = dlmread(fre_file);
%[eigenVec_remove, eigenV_remove] = eig(Hmatrix_remove, Mmatrix_remove);
%[eigenValues_remove,ind_remove] = sort(diag(eigenV_remove));
%eigenVectors_remove = eigenVec_remove(:,ind_remove);

%fre(fre< 10^-7) = 0;
fre = sqrt(fre);

contact = csvread(contact_file);

%coordinate_file = folder + "coordinate_file_remove" + extend;
coordinate_file = folder + "coordinate_file" + extend;

%radius_file = folder + "radius_remove" + extend;
radius_file = folder + "radius" + extend;

v_file =  folder + "v_file_remove" + extend;
dt_file = folder + "dt" + extend;
p_and_phi_file = folder + "pressure" + extend;

coordinate = csvread(coordinate_file);
radius = csvread(radius_file);
p_and_phi = csvread(p_and_phi_file);
pressure = p_and_phi(:,1);
phi = p_and_phi(:,2);

%N=5;
N=16;
frames= size(coordinate,1)/N ;

for ii = 1 : round(frames/10):frames
    start_point = 1 + N * ( ii - 1 );
    end_point = N * ii;
    plot_particles_2d(1,[1,1],radius(start_point:end_point),coordinate(start_point:end_point,1),coordinate(start_point:end_point,2))
    pause(0.2);
end

% vobj = VideoWriter('decompress.mp4','MPEG-4');
% vobj.FrameRate = 100;
% open(vobj);
% for i = 1 : frames
%     start_point = 1 + N * ( i - 1 );
%     end_point = N * i;
%     plot_particles_2d(1,[1,1],radius(start_point:end_point),coordinate(start_point:end_point,1),coordinate(start_point:end_point,2))
%     frame = getframe(gcf) ;
%     writeVideo(vobj, frame);
% end
% 
% close(vobj);

figure(2);
plot(phi,pressure);
set(gca, 'XDir','reverse')
xlabel("phi");
ylabel("pressure");

% figure(3);
% plot(phi,fre(3:(2*N):end-2*N));
% xlabel("phi");
% ylabel("w3");
% 
% figure(4);
% plot(phi,fre(4:(2*N):end-2*N));
% xlabel("phi");
% ylabel("w3");


% index = find(phi<0.714 & phi>0.712);
% index1 = (index-1)*10+3;
% figure(5);
% plot(pressure(index),fre(index1));
% xlabel("pressure");
% ylabel("w3");

% [c ,index] = min(abs(phi-0.5382));
% ii = index;
% start_point = 1 + N * ( ii - 1 );
% end_point = N * ii;
% plot_particles_2d(3,[1,1],radius(start_point:end_point),coordinate(start_point:end_point,1),coordinate(start_point:end_point,2))
% 
% 
% 
% [c ,index] = min(abs(phi-0.3882));
% ii = index;
% start_point = 1 + N * ( ii - 1 );
% end_point = N * ii;
% plot_particles_2d(4,[1,1],radius(start_point:end_point),coordinate(start_point:end_point,1),coordinate(start_point:end_point,2))
% 
% [c ,index] = min(abs(phi-0.4972));
% ii = index;
% start_point = 1 + N * ( ii - 1 );
% end_point = N * ii;
% plot_particles_2d(5,[1,1],radius(start_point:end_point),coordinate(start_point:end_point,1),coordinate(start_point:end_point,2))
% 

q_pressure = 10.^(-(1:7));
q_frequency = [0.0362841,0.00348675,0.000337248,2.236569657e-05,2.691877618e-06,2.047437992e-07,9.295732479e-09];

figure(5);
plot(log10(q_pressure) ,log10(sqrt(q_frequency)));
xlabel("log(pressure)");
ylabel("log(w3)");


mean_contact = sum(contact>0,2);
mean_contact_r = sum(contact==1,2);
mean_contact_a = sum(contact==2,2);   
iii = size(mean_contact,1);
% contact_number = [contact_number, mean_contact];
% contact_number_r = [contact_number_r, mean_contact_r];
% contact_number_a = [contact_number_a, mean_contact_a];




figure(6); hold on
plot(phi, mean_contact);
%plot(1:iii, flip(mean_contact_r));
plot(phi, mean_contact_a);
set(gca, 'XDir','reverse')
xlabel("phi");
ylabel("contact number");

figure(7);
per = (mean_contact_a./mean_contact);
plot(phi, per);
set(gca, 'XDir','reverse')
xlabel("phi");
ylabel("ratio");



