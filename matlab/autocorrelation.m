clc
clear
close all

folder = "C:\Users\Yuxuan Cheng\OneDrive\Windows_Version\remove1\";
extend = "_1.txt";

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
v = csvread(v_file);
dt = csvread(dt_file);

N=8;
frames= size(coordinate,1)/N ;

for ii = 1 : round(frames/50):frames
    start_point = 1 + N * ( ii - 1 );
    end_point = N * ii;
    plot_particles_2d(1,[1,1],radius(start_point:end_point),coordinate(start_point:end_point,1),coordinate(start_point:end_point,2))
    
end

vx = v(:,1);
vy = v(:,2);
vx = reshape(vx,[N,frames]);
vy = reshape(vy,[N,frames]);
v = sqrt(vx.^2 +vy.^2);

coordinate_x = coordinate(:,1);
coordinate_y = coordinate(:,2);
coordinate_x = reshape(coordinate_x,[N,frames]);
coordinate_y = reshape(coordinate_y,[N,frames]);


autocorre = cal_autocorrelation(frames, vx, vy, v, 0);

% wv = [];
% %w_range = 0.001:0.0005:0.6; 
% w_range = 0.001:0.1:250; 
% for w = w_range 
%     auto = autocorre .* exp(sqrt(-1) * 0.01 * 20 * w * (1:frames-1))';
%     wv = [wv, sum(auto,'all');];
% end
% figure(2)
% plot(w_range,wv)

wv = [];

%w_range = 0.001:0.1:250; 
w_range = 1:0.001:30; 
for w = w_range 
    %auto = autocorre .* cos(4.321028*10^-5 * 100 * w * (1:frames-1))';
    auto = autocorre .* cos(dt * 100 * w * (1:frames-1))';
    %auto = autocorre .* cos(0.0003848 * 10 * w * (1:frames-1))';
    wv = [wv, sum(auto,'all');];
end
figure(2); hold on;
plot(w_range,log10(abs(wv)))
for i = 3:length(fre)
    xline(fre(i),'--r','color','b');
end


v_matrix = zeros(2*N);
for ii = 1: 2*N
    index1 = ceil(ii/2);
    
    if mod(ii,2)==1
        v1 = vx(index1,:);
    else
        v1 = vy(index1,:);
    end
    
    for jj = 1: 2*N
        
        index2 = ceil(jj/2);
    
        if mod(jj,2)==1
            v2 = vx(index2,:);
        else
            v2 = vy(index2,:);
        end
        
        v_matrix(ii,jj) = mean(v1.*v2,"all");
          
        
    end
end

c_matrix = zeros(2*N);
for ii = 1: 2*N
    
    index1 = ceil(ii/2);
    
    if mod(ii,2)==1
        brac1 = coordinate_x(index1,:)-coordinate_x(index1,1);
    else
        brac1 = coordinate_y(index1,:)-coordinate_y(index1,1);
    end
    
    for jj = 1: 2*N
        
        index2 = ceil(jj/2);
    
        if mod(jj,2)==1
            brac2 = coordinate_x(index2,:)-coordinate_x(index2,1);
        else
            brac2 = coordinate_y(index2,:)-coordinate_y(index2,1);
        end
        
        c_matrix(ii,jj) = mean(brac1.*brac2,"all");
        
    end
end

[eigenVec, eigenV] = eig(v_matrix/c_matrix);
[eigenValues,ind] = sort(diag(eigenV));
eigenVectors = eigenVec(:,ind);

figure(3); hold on
scatter((1:(2*N-2))/(2*N-2),sqrt(abs(eigenValues(3:end))'),20);
scatter((1:(2*N-2))/(2*N-2),fre(3:end),100,'s');
xlabel("i/2N-2");
ylabel("w");
legend("displacement", "Hessian")

function autocorrelation = cal_autocorrelation(NT, vx, vy, v, skip)
% loop over the different possible time windows, calculate MSD for each
% time window size

    % create MSD array (y-axis of MSD plot)
    % time = round(logspace(0,log10(NT-skip-1),100));
    
    autocorrelation = zeros(NT - 1,1);
    
    for ii = 1: NT-1
        % calculate x displacements, separated by ii indices
        dotv = (vx(:,1+skip+ii:end) .* vx(:,1+skip:end-ii)) + ...
            (vy(:,1+skip+ii:end) .* vy(:,1+skip:end-ii));

        % take mean over all displacements
        corre = mean(dotv, 'all')/mean(v(:,1:(NT - ii)).^2,'all');
       
        % store in MSD array
        autocorrelation(ii) = corre;
              
    end
    

end








