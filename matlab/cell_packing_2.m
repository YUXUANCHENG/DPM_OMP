clc
clear
close all


%folder = "C:\Users\Yuxuan Cheng\source\repos\cell\Debug\New folder5\";
%folder = "C:\Users\Yuxuan Cheng\source\repos\cell_packing\cell_packing\";
%folder = "C:\Users\Yuxuan Cheng\source\repos\cell_packing\Debug\activity6\";
folder = "C:\Users\Yuxuan Cheng\source\repos\cells\forked-cells\forked-cells\";
%folder = "D:\cells\bash\yc\";

all_mean_cal_A = [];
contact_file = folder + "contact.txt";
packing_contact = dlmread(contact_file);
coordinate_file = folder + "jam.txt";
coordinate = csvread(coordinate_file);
%length_file = folder + "length" + extend;
length_file = folder + "length.txt";
lengthscale = csvread(length_file);
cal_A_file = folder + "calA.txt";
cal_A = csvread(cal_A_file);
N=sum(lengthscale(1:end-2),'all');
frames= size(coordinate,1)/N ;
Ncell = size(lengthscale,1)-2;

i = frames;
start_point = 1 + N * ( i - 1 );
end_point = N * i;
plot_particles_2d(1,[lengthscale(end-1),lengthscale(end)],...
coordinate(start_point:end_point,3)/2,coordinate(start_point:end_point,1),coordinate(start_point:end_point,2))

xpos_at_frame = coordinate(start_point:end_point,1);
ypos_at_frame = coordinate(start_point:end_point,2);
[vAll,cAll] = draw_voronoi(xpos_at_frame, ypos_at_frame, Ncell, lengthscale, 3);
jammed_voronoi_con_net = voronoi_contact(cAll, Ncell);

v0_file = folder + "v0.txt";
v0 = csvread(v0_file);
order_per = [];

for t_index = 0:9
    try
    extend = "_" + int2str(t_index) +".txt";

    coordinate_file = folder + "jam" + extend;
    coordinate = csvread(coordinate_file);
    %length_file = folder + "length" + extend;
    length_file = folder + "length.txt";
    lengthscale = csvread(length_file);
    cal_A_file = folder + "calA" + extend;
    cal_A = csvread(cal_A_file);
    contact_file = folder + "contact" + extend;
    contact = dlmread(contact_file);
    
    catch
        continue
    end
    
    contact = contact(end,:);
    difference = sum(contact ~= packing_contact,'all');
%     if difference > 0
%         disp(t_index)
%         disp(difference)
%     end
    
    N=sum(lengthscale(1:end-2),'all');
    frames= size(coordinate,1)/N ;
    Ncell = size(lengthscale,1)-2;

%     if t_index == 9
%         vobj = VideoWriter('active_cell3.mp4','MPEG-4');
%         vobj.FrameRate = 20;
%         open(vobj);
%             for i = 1 :  round(frames/100):frames
%                 start_point = 1 + N * ( i - 1 );
%                 end_point = N * i;
%                 plot_particles_2d(2,[lengthscale(end-1),lengthscale(end)],...
%                     coordinate(start_point:end_point,3)/2,coordinate(start_point:end_point,1),coordinate(start_point:end_point,2))
%                 frame = getframe(gcf) ;
%                 writeVideo(vobj, frame);
%             end
%         
%         close(vobj);
%     end
    
    for i = 1 :  round(frames/5):frames
        start_point = 1 + N * ( i - 1 );
        end_point = N * i;
        plot_particles_2d(2,[lengthscale(end-1),lengthscale(end)],...
            coordinate(start_point:end_point,3)/2,coordinate(start_point:end_point,1),coordinate(start_point:end_point,2))
    %     frame = getframe(gcf) ;
    %     writeVideo(vobj, frame);
    end
   
    
    i = frames;
    start_point = 1 + N * ( i - 1 );
    end_point = N * i;
    plot_particles_2d((t_index+1)*10,[lengthscale(end-1),lengthscale(end)],...
    coordinate(start_point:end_point,3)/2,coordinate(start_point:end_point,1),coordinate(start_point:end_point,2))
    xpos_at_frame = coordinate(start_point:end_point,1);
    ypos_at_frame = coordinate(start_point:end_point,2);
    [vAll,cAll] = draw_voronoi(xpos_at_frame, ypos_at_frame, Ncell, lengthscale, (t_index+1)*10+1);
    voronoi_con_net = voronoi_contact(cAll, Ncell);
    is_equal = isequal(voronoi_con_net,jammed_voronoi_con_net);
    if is_equal == 0
        disp({'voronoi_con_net is different for ',t_index})
    end
    order = cal_order(coordinate(:,1:2),N)/(0.005*100)/v0(t_index+1,1);
    order_per = [order_per, order];    

end

order_per = reshape(order_per, 10, []);
% figure(3);
% plot(v0,all_mean_cal_A)
% xlabel("v0");
% ylabel("CalA");

function [vAll,cAll] = draw_voronoi(xpos_at_frame, ypos_at_frame, Ncell, lengthscale, fig)
    xcomp = zeros(1,Ncell);
    ycomp = zeros(1,Ncell);
    for ci = 1:Ncell
        index = sum(lengthscale(1:ci),'all');
        start_point_last = 1 + index - lengthscale(ci);
        end_point_last = index;
        cx_tmp = mean(xpos_at_frame(start_point_last:end_point_last),'all');
        cy_tmp = mean(ypos_at_frame(start_point_last:end_point_last),'all');
%         cx_tmp = cx_tmp - lengthscale(end-1) * round(cx_tmp/lengthscale(end-1));
%         cy_tmp = cy_tmp - lengthscale(end) * round(cy_tmp/lengthscale(end));
%         if cx_tmp <0
%             cx_tmp = cx_tmp + lengthscale(end-1);
%         end
%         if cy_tmp <0
%             cy_tmp = cy_tmp + lengthscale(end);
%         end
        xcomp(ci) = mod(cx_tmp,lengthscale(end-1));
        ycomp(ci) = mod(cy_tmp,lengthscale(end));
    end
    voronoiX = [];
    voronoiY = [];
    for i = -1:1
        for j = -1:1
            voronoiX = [voronoiX,xcomp + lengthscale(end-1) * i];
            voronoiY = [voronoiY,ycomp + lengthscale(end)* j];
        end
    end
    text_range = (1 + 4 * Ncell) : (5 * Ncell);
    figure(fig); hold on
    voronoi(voronoiX,voronoiY)
    plabels = arrayfun(@(n) {sprintf('C%d', mod(n-1,Ncell)+1)}, text_range');
    text(voronoiX(text_range), voronoiY(text_range), plabels, 'FontWeight', ...
          'bold', 'HorizontalAlignment','center', ...
          'BackgroundColor', 'none')
    xlim([0,lengthscale(end-1)]);
    ylim([0,lengthscale(end)]);
    axis square
    hold off
    dt = delaunayTriangulation(voronoiX',voronoiY');
    [vAll,cAll] = voronoiDiagram(dt);
end

function voronoi_con_net = voronoi_contact(cAll, Ncell)
    voronoi_con_net = {};
    for i = (1 + 4 * Ncell) : (5 * Ncell)
        mod_i = mod(i-1,Ncell)+1;
        neighbers = [];
        for element =  cAll{i}
            for j = 1:Ncell*9
                mod_j = mod(j-1,Ncell)+1;
                element_of_other_cell = cAll{j};
                if (sum(element_of_other_cell==element)>0 && mod_j > mod_i)
                    neighbers = [neighbers, mod_j];
                end
            end
        end
        voronoi_con_net{mod_i} = sort(unique(neighbers));
    end
end

function order = cal_order(coordinate,N)

    vel = coordinate(1:end - N,:) - coordinate(N+1:end,:);
    
    order = mean(vel,1);
    order = norm(order);

end




