clc
clear
close all

%folder = "C:\Users\Yuxuan Cheng\source\repos\cells\forked-cells\forked-cells\";
%folder = "C:\Users\Yuxuan Cheng\source\repos\cell_packing\Debug\activity6\";
%folder = "D:\project\cells4\";
folder = "D:\cells21\";
%folder = "D:\project\cells2_32\";

all_mean_cal_A = [];
contact_file = folder + "contact.txt";
packing_contact = dlmread(contact_file);
packing_contact = packing_contact(end,:);
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

% for i = 1 :frames
%     start_point = 1 + N * ( i - 1 );
%     end_point = N * i;
%     plot_particles_2d(2,[lengthscale(end-1),lengthscale(end)],...
%         coordinate(start_point:end_point,3)/2,coordinate(start_point:end_point,1),coordinate(start_point:end_point,2))
% %         frame = getframe(gcf) ;
% %         writeVideo(vobj, frame);
% end

% vobj = VideoWriter('bring_to_jam.mp4','MPEG-4');
% vobj.FrameRate = 2;
% open(vobj);
% for i = 2 :frames
%     start_point = 1 + N * ( i - 1 );
%     end_point = N * i;
%     plot_particles_2d(2,[lengthscale(end-1),lengthscale(end)],...
%         coordinate(start_point:end_point,3)/2,coordinate(start_point:end_point,1),coordinate(start_point:end_point,2))
%     frame = getframe(gcf) ;
%     writeVideo(vobj, frame);
% end
% close(vobj);

i = frames;
start_point = 1 + N * ( i - 1 );
end_point = N * i;
plot_particles_2d(1,[lengthscale(end-1),lengthscale(end)],...
coordinate(start_point:end_point,3)/2,coordinate(start_point:end_point,1),coordinate(start_point:end_point,2))
 
xpos_at_frame = coordinate(start_point:end_point,1);
ypos_at_frame = coordinate(start_point:end_point,2);
[vAll,cAll,xcomp_j,ycomp_j] = draw_voronoi(xpos_at_frame, ypos_at_frame, Ncell, lengthscale, 3);
jammed_voronoi_con_net = voronoi_contact(cAll, Ncell);

v0_file = folder + "v0.txt";
v0 = csvread(v0_file);
order_per = [];
ifjammed = [];
msd = {};


for t_index_i = 0:9
    %close all
    for t_index_j = 0:9
    try
    extend = "_" + int2str(t_index_i) + int2str(t_index_j) +".txt";

    coordinate_file = folder + "jam" + extend;
    coordinate = csvread(coordinate_file);
    %length_file = folder + "length" + extend;
    length_file = folder + "length.txt";
    lengthscale = csvread(length_file);
    cal_A_file = folder + "calA" + extend;
    cal_A = csvread(cal_A_file);
    contact_file = folder + "contact" + extend;
    contact = dlmread(contact_file);
    v_file = folder + "v" + extend;
    vel = csvread(v_file);
    
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

%     vobj = VideoWriter('low_kl.mp4','MPEG-4');
%     vobj.FrameRate = 20;
%     open(vobj);
%         for i = 1 :  round(frames/100):frames
%             start_point = 1 + N * ( i - 1 );
%             end_point = N * i;
%             plot_particles_2d(2,[lengthscale(end-1),lengthscale(end)],...
%                 coordinate(start_point:end_point,3)/2,coordinate(start_point:end_point,1),coordinate(start_point:end_point,2))
%             frame = getframe(gcf) ;
%             writeVideo(vobj, frame);
%         end
% 
%     close(vobj);

%    
%     for i = 1 : round(frames/10):frames
%         start_point = 1 + N * ( i - 1 );
%         end_point = N * i;
%         plot_particles_2d(2,[lengthscale(end-1),lengthscale(end)],...
%             coordinate(start_point:end_point,3)/2,coordinate(start_point:end_point,1),coordinate(start_point:end_point,2))
%     end
    
    i = frames;
    start_point = 1 + N * ( i - 1 );
    end_point = N * i;
    plot_particles_2d((t_index_j+1)*10,[lengthscale(end-1),lengthscale(end)],...
        coordinate(start_point:end_point,3)/2,coordinate(start_point:end_point,1),coordinate(start_point:end_point,2))
    xpos_at_frame = coordinate(start_point:end_point,1);
    ypos_at_frame = coordinate(start_point:end_point,2);
    [vAll,cAll,xcomp,ycomp] = draw_voronoi(xpos_at_frame, ypos_at_frame, Ncell, lengthscale, (t_index_j+1)*10+1);
    voronoi_con_net = voronoi_contact(cAll, Ncell);
   
    
    v_d_index = t_index_i * 10 + t_index_j + 1;
    
%     order = cal_order1(vel)/v0(v_d_index,1);
%     order_per = [order_per, order];
    
    
%     n_particle = density_fluctuation(N, coordinate, frames, Ncell, lengthscale);
%     N_std = std(n_particle);
%     mean_N = mean(n_particle);
%     disp(N_std/sqrt(mean_N))
    
    [MSD,deltaT] = cal_msd(N, coordinate, frames, Ncell, lengthscale);
    % open figure window
    figure((t_index_j+1)*10+2), clf, hold on, box on;
    % plot curve, add units to axes, etc
    plot(deltaT, MSD,'color','red','linewidth',3);
    xlabel('time');ylabel('MSD');
    %P = polyfit(log10(deltaT(round(5*frames/6): end)), log10(MSD(round(5*frames/6): end))', 1);
    P = polyfit(log10(deltaT(1: end)), log10(MSD(1: end))', 1);
    yfit = P(1)*log10(deltaT)+P(2);
    plot(deltaT,10.^(yfit),'r-.');
    theString = sprintf('slope = %.3f ', P(1));
    text(50, 0.01, theString, 'FontSize', 24);
    ax = gca;
    ax.FontSize = 22;
    ax.XScale = "log";
    ax.YScale = "log";
    
    msd{t_index_i+1,t_index_j+1} = MSD;
    
    %is_equal = isequal(voronoi_con_net,jammed_voronoi_con_net);
    %is_equal = compare_contact(voronoi_con_net,jammed_voronoi_con_net,Ncell);
    %is_equal = compare_contact1(voronoi_con_net,jammed_voronoi_con_net,xcomp_j,ycomp_j,lengthscale,Ncell);
    is_equal = (max(MSD,[],'all')<sqrt(lengthscale(end)*lengthscale(end-1)/(Ncell*3.14)));
    
    ifjammed = [ifjammed, is_equal];
    
    if is_equal == 0
        disp({'voronoi_con_net is different for ',t_index_i,t_index_j})
    end
    
    cal_A = reshape(cal_A,Ncell,[]);
    mean_cal_A = mean(cal_A,'all');
%     
    all_mean_cal_A = [all_mean_cal_A, mean_cal_A];
    end
end

% order_per = reshape(order_per, 10, []);

% figure(5);
% heatmap(order_per)
% xlabel("v0");
% ylabel("noise");
% title("order parameter");

figure(6);
ifjammed = reshape(ifjammed, 10, []);
heatmap(ifjammed,'CellLabelColor','none')
xlabel("v0");
ylabel("Kb");


% figure(7);
% plot(order_per(:,1))
% ylabel("order");
% xlabel("noise");

figure(8), clf, hold on, box on;
% plot curve, add units to axes, etc
plot(deltaT, msd{1,1},'color','red','linewidth',3);
plot(deltaT, msd{1,2},'color','blue','linewidth',3);
%P = polyfit(log10(deltaT(round(5*frames/6): end)), log10(MSD(round(5*frames/6): end))', 1);
MSD1=msd{1,1};
P = polyfit(log10(deltaT(1: end)), log10(MSD1(1: end))', 1);
yfit = P(1)*log10(deltaT)+P(2);
plot(deltaT,10.^(yfit),'r-.');
legend({'Liquid like','Solid like','fit'});
xlabel('time');ylabel('MSD');
theString = sprintf('slope = %.3f ', P(1));
text(5, 1, theString, 'FontSize', 24);
ax = gca;
ax.FontSize = 22;
ax.XScale = "log";
ax.YScale = "log";
    
% v0_file = folder + "v0.txt";
% v0 = csvread(v0_file);
% figure(3);
% plot(v0,all_mean_cal_A)
% xlabel("v0");
% ylabel("CalA");


function [vAll,cAll,xcomp,ycomp] = draw_voronoi(xpos_at_frame, ypos_at_frame, Ncell, lengthscale, fig)
    xcomp = zeros(1,Ncell);
    ycomp = zeros(1,Ncell);
    for ci = 1:Ncell
        index = sum(lengthscale(1:ci),'all');
        start_point_last = 1 + index - lengthscale(ci);
        end_point_last = index;
        cx_tmp = mean(xpos_at_frame(start_point_last:end_point_last),'all');
        cy_tmp = mean(ypos_at_frame(start_point_last:end_point_last),'all');
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
    clf;
    voronoi(voronoiX,voronoiY)
    plabels = arrayfun(@(n) {sprintf('C%d', mod(n-1,Ncell)+1)}, text_range');
    text(voronoiX(text_range), voronoiY(text_range), plabels, 'FontWeight', ...
          'bold', 'HorizontalAlignment','center', ...
          'BackgroundColor', 'none')
    xlim([0,lengthscale(end-1)]);
    ylim([0,lengthscale(end)]);
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
                if (sum(element_of_other_cell==element)>0 && mod_j ~= mod_i)
                    neighbers = [neighbers, mod_j];
                end
            end
        end
        voronoi_con_net{mod_i} = sort(unique(neighbers));
    end
end

function order = cal_order(coordinate,N)

    vel = coordinate(1:end - N,:) - coordinate(N+1:end,:);
    velx = reshape(vel(:,1) , N, []);
    vely = reshape(vel(:,2) , N, []);
    vel_norm = sqrt(velx.^2+ vely.^2);
    velx = mean(velx./vel_norm,1);
    vely = mean(vely./vel_norm,1);
    
%     velx = mean(velx,1);
%     vely = mean(vely,1);
    
    order = sqrt(velx.^2+ vely.^2);
    order = mean(order,'all');
      
end

function order = cal_order1(vel)

    order = sqrt(vel(:,1).^2+ vel(:,2).^2);
    order = mean(order,'all');
      
end




function isequal = compare_contact(voronoi_con_net,jammed_voronoi_con_net,Ncell)
    isequal = 1;
    for i = 1:Ncell
        if sum(ismember(jammed_voronoi_con_net{i},voronoi_con_net{i}),'all') < 3
            isequal = 0;
        end
    end
end

function isequal = compare_contact1(voronoi_con_net,jammed_voronoi_con_net,xcomp,ycomp,lengthscale,Ncell)
    isequal = 1;
    for i = 1:Ncell
        changedcontact = find(~ismember(voronoi_con_net{i},jammed_voronoi_con_net{i}));
        for index = changedcontact
            dis1 = cell_distance(i,voronoi_con_net{i}(index),xcomp,ycomp,lengthscale);
            dis2 = mean(cell_distance(i,jammed_voronoi_con_net{i},xcomp,ycomp,lengthscale),'all');
            if dis1> 2 * dis2
                isequal = 0;
            end
        end
    end
end

function distance = cell_distance(i,j,xcomp,ycomp,lengthscale)
    ci_x = xcomp(i);
    ci_y = ycomp(i);
    cj_x = xcomp(j);
    cj_y = ycomp(j);
    dx = ci_x - cj_x;
    dy = ci_y - cj_y;
    dx = dx - lengthscale(end-1)*round(dx./lengthscale(end-1));
    dy = dy - lengthscale(end)*round(dy./lengthscale(end));
    distance = sqrt(dx.^2+dy.^2);
end

function n_particle = density_fluctuation(N, coordinate, frames, Ncell, lengthscale)
    xcomp=zeros(frames,Ncell);
    ycomp=zeros(frames,Ncell);
    for i = 1 :frames
        start_point = 1 + N * ( i - 1 );
        end_point = N * i;
        [xcomp_t,ycomp_t]=cal_c_pos(coordinate(start_point:end_point,1),coordinate(start_point:end_point,2), Ncell, lengthscale);
        xcomp(i,:)= xcomp_t;
        ycomp(i,:)= ycomp_t;
    end
    
    xcomp = mod(xcomp,lengthscale(end-1));
    ycomp = mod(ycomp,lengthscale(end));
    n_particle = xcomp<3*lengthscale(end-1)/6 &  xcomp>lengthscale(end-1)/6 & ycomp<3*lengthscale(end)/6 & ycomp>lengthscale(end)/6;
    n_particle = sum(n_particle,2);
    
end



function [xcomp,ycomp]=cal_c_pos(xpos_at_frame, ypos_at_frame, Ncell, lengthscale)
    xcomp = zeros(1,Ncell);
    ycomp = zeros(1,Ncell);
    for ci = 1:Ncell
        index = sum(lengthscale(1:ci),'all');
        start_point_last = 1 + index - lengthscale(ci);
        end_point_last = index;
        cx_tmp = mean(xpos_at_frame(start_point_last:end_point_last),'all');
        cy_tmp = mean(ypos_at_frame(start_point_last:end_point_last),'all');
        xcomp(ci) = cx_tmp;
        ycomp(ci) = cy_tmp;
    end
end


function [MSD,deltaT] = cal_msd(N, coordinate, frames, Ncell, lengthscale)


    xcomp=zeros(frames,Ncell);
    ycomp=zeros(frames,Ncell);
    for i = 1 :frames
        start_point = 1 + N * ( i - 1 );
        end_point = N * i;
        [xcomp_t,ycomp_t]=cal_c_pos(coordinate(start_point:end_point,1),coordinate(start_point:end_point,2), Ncell, lengthscale);
        xcomp(i,:)= xcomp_t;
        ycomp(i,:)= ycomp_t;
    end
    
    
    % create MSD array (y-axis of MSD plot)
    MSD = zeros(round(9*frames/10),1);
    NT = length(MSD);
    % loop over the different possible time windows, calculate MSD for each
    % time window size
    for ii = 1:NT

        % calculate x displacements, separated by ii indices
        dx = xcomp(1+ii:end,:) - xcomp(1:end-ii,:);

        % calculate y displacements similarly
        dy = ycomp(1+ii:end,:) - ycomp(1:end-ii,:);

        % take mean over all displacements
        dispMean = mean(dx.^2 + dy.^2,'all');

        % store in MSD array
        MSD(ii) = dispMean;
    end

    % create deltaT array, using a for loop or vectorization
    deltaT = 1:NT;


end
















