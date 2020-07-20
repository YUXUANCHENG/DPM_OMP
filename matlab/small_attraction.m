clear;
clc;
close all;

%coordinate_file = '/Users/yc757/Downloads/test/test/DerivedData/test/Build/Products/Debug/coordinate_file4.txt';
%radius_file = '/Users/yc757/Downloads/test/test/DerivedData/test/Build/Products/Debug/radius4.txt';
%contact_file = '/Users/yc757/Downloads/test/test/DerivedData/test/Build/Products/Debug/contact4.txt';

coordinate_file = "C:\Users\Yuxuan Cheng\OneDrive\Windows_Version\movie_5/coordinate_file_5.txt";
radius_file = "C:\Users\Yuxuan Cheng\OneDrive\Windows_Version\movie_5/radius_5.txt";
contact_file = "C:\Users\Yuxuan Cheng\OneDrive\Windows_Version\movie_5/contact_5.txt";

coordinate = csvread(coordinate_file);
radius = csvread(radius_file);
contact = csvread(contact_file);

N=64;
frames= size(coordinate,1)/N ;


if (N * (N-1) /2 ~= size(contact,2))
    print("stop");
end

size_evolution = [];

for i = 1 : size(contact,1)
    A = transfer_A(N, contact(i,:));
    [nComponents,sizes,members] = networkComponents(A);
    size_evolution = [size_evolution, sizes(1)];
end
figure(2)
scatter(1:size(size_evolution,2),size_evolution)



% vobj = VideoWriter('attraction4.mp4','MPEG-4');
% vobj.FrameRate = 20;
% open(vobj);
% for i = 1 : frames
%     start_point = 1 + N * ( i - 1 );
%     end_point = N * i;
%     plot_particles_2d(1,[1,1],radius(start_point:end_point),coordinate(start_point:end_point,1),coordinate(start_point:end_point,2))
% 
%     frame = getframe(gcf) ;
%     writeVideo(vobj, frame);
% end
% 
% close(vobj);










