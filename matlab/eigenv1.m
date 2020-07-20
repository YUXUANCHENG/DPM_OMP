clc
clear
close all


folder = "C:\Users\Yuxuan Cheng\source\repos\distritution\Debug\negative4\";

t_index = 1;

extend = "_" + int2str(t_index) +".txt";
contact_file = folder + "contact" + extend;
contact = csvread(contact_file);
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

for ii = 1 : round(frames/10):frames
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


%     autocorre = cal_autocorrelation(frames, vx, vy, v, 0);
% 
%     wv = [];
% 
%     w_range = 1:0.001:30; 
%     for w = w_range 
%         auto = autocorre .* cos(dt * 100 * w * (1:frames-1))';
%         wv = [wv, sum(auto,'all');];
%     end
%     figure(2); hold on;
%     plot(w_range,log10(abs(wv)))
%     for i = 3:length(fre)
%         xline(fre(i),'--r','color','b');
%     end


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

[eigenVec1, eigenV1] = eig(v_matrix*pinv(c_matrix));
[eigenValues1,ind1] = sort(diag(eigenV1));
eigenVectors1 = eigenVec(:,ind1);

figure(4); hold on
scatter((1:(2*N-2))/(2*N-2),sqrt(abs(eigenValues1(3:end))'),20);
scatter((1:(2*N-2))/(2*N-2),fre(3:end),100,'s');
xlabel("i/2N-2");
ylabel("w");
legend("displacement", "Hessian")



len = size(eigenVectors,1);
x_tran = zeros(1,len);
x_tran(1:2:end)=1;
x_tran = x_tran'/norm(x_tran);

y_tran = zeros(1,len);
y_tran(2:2:end)=1;
y_tran = y_tran'/norm(y_tran);

e={};

for i = 1:size(eigenVectors_remove,1)
    e{i} = eigenVectors(:,i)/norm(eigenVectors(:,i));
    % ceil(find(abs(round(e{i}))>0.5)/2)
end

e{1} = x_tran;
e{2} = y_tran;

e = GramSchmidt(e);

%coordinate = csvread(coordinate_file);
%radius = csvread(radius_file);


N=8;
frames= size(coordinate,1)/N ;


if (N * (N-1) /2 ~= size(contact,2))
    print("stop");
end

for i = 1 : 1
    start_point = 1 + N * ( i - 1 );
    end_point = N * i;
    plot_particles_2d(1,[1,1],radius(start_point:end_point),coordinate(start_point:end_point,1),coordinate(start_point:end_point,2))
end

ee=e{3};
quiver(coordinate(start_point:end_point,1),coordinate(start_point:end_point,2),ee(1:2:end),ee(2:2:end),'color','r','linewidth',3);



function v = GramSchmidt(v)

k = size(v,2);
assert(k>=2,'The input matrix must include more than one vector.');
for ii = 1:k
    v{ii} = v{ii} / norm(v{ii});
    for jj = ii+1:k
        v{jj} = v{jj} - proj(v{jj},v{ii});
        v{jj} = v{jj} / norm(v{jj});
    end
end
end

function w = proj(u,v)
    % This function projects vector v on vector u
    w = (u'*v).*v;
end































