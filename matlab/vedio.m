clc
clear
close all


folder = "C:\Users\Yuxuan Cheng\source\repos\distritution\Debug\remove5\";
t_index = 5;

extend = "_" + int2str(t_index) +".txt";
H_matrix_file_remove = folder + "H_matrix_remove" + extend;
Hmatrix_remove = dlmread(H_matrix_file_remove);
M_matrix_file_remove = folder + "M_matrix_remove" + extend;
Mmatrix_remove = dlmread(M_matrix_file_remove);
[eigenVec_remove, eigenV_remove] = eig(Hmatrix_remove, Mmatrix_remove);
[eigenValues_remove,ind_remove] = sort(diag(eigenV_remove));
eigenVectors_remove = eigenVec_remove(:,ind_remove);

eigenValues_remove(eigenValues_remove< 10^-7) = 0;
fre = sqrt(eigenValues_remove);

coordinate_file = folder + "coordinate_file_remove" + extend;
radius_file = folder + "radius_remove" + extend;
v_file =  folder + "v_file_remove" + extend;
dt_file = folder + "dt" + extend;

coordinate = csvread(coordinate_file);
radius = csvread(radius_file);

N=8;
frames= size(coordinate,1)/N ;

% for ii = 1 : round(frames/100):frames
%     start_point = 1 + N * ( ii - 1 );
%     end_point = N * ii;
%     plot_particles_2d(1,[1,1],radius(start_point:end_point),coordinate(start_point:end_point,1),coordinate(start_point:end_point,2))
%     %pause(0.2);
% end

vobj = VideoWriter('vibration1.mp4','MPEG-4');
vobj.FrameRate = 100;
open(vobj);
for i = 1 :  round(frames/500):frames
    start_point = 1 + N * ( i - 1 );
    end_point = N * i;
    plot_particles_2d(1,[1,1],radius(start_point:end_point),coordinate(start_point:end_point,1),coordinate(start_point:end_point,2))
    frame = getframe(gcf) ;
    writeVideo(vobj, frame);
end

close(vobj);










