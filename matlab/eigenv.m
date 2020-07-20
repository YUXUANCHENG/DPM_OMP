clc
clear
close all


folder = "C:\Users\Yuxuan Cheng\source\repos\distritution\Debug\negative1\";
%folder = "C:\Users\Yuxuan Cheng\source\repos\distritution\Debug\decompress7\";
t_index = 1;
extend = "_" + int2str(t_index) +".txt";
contact_file = folder + "contact" + extend;
contact = csvread(contact_file);
H_matrix_file_remove = folder + "H_matrix_remove" + extend;
Hmatrix_remove = dlmread(H_matrix_file_remove);
M_matrix_file_remove = folder + "M_matrix_remove" + extend;
Mmatrix_remove = dlmread(M_matrix_file_remove);
H_matrix_file = folder + "H_matrix" + extend;
Hmatrix = dlmread(H_matrix_file);
M_matrix_file = folder + "M_matrix" + extend;
Mmatrix = dlmread(M_matrix_file);

%coordinate_file = folder + "coordinate_file_remove" + extend;
coordinate_file = folder + "coordinate_file" + extend;
%radius_file = folder + "radius_remove" + extend;
radius_file = folder + "radius" + extend;

[eigenVec, eigenV] = eig(Hmatrix, Mmatrix);
[eigenValues,ind] = sort(diag(eigenV));
eigenVectors = eigenVec(:,ind);

[eigenVec_remove, eigenV_remove] = eig(Hmatrix_remove, Mmatrix_remove);
[eigenValues_remove,ind_remove] = sort(diag(eigenV_remove));
eigenVectors_remove = eigenVec_remove(:,ind_remove);


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

for i = 1:5
    ceil(find(abs(round(e{i}))>0.5)/2)
end

coordinate = csvread(coordinate_file);
radius = csvread(radius_file);


N=5;
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


% eigenValues_remove(eigenValues_remove< 10^-7) = 0;
% figure(2);
% [counts, binCenters] = hist(sqrt(eigenValues_remove),10);
% counts = counts / sum(counts);
% plot(binCenters, counts,'linewidth',1.75);
% xlabel("w");
% ylabel("D(w)");
% figure(3);
% scatter((1:2*N)/(2*N),sqrt(eigenValues_remove));
% xlabel("i/2N");
% ylabel("w");
% 
% 
% pariticipation = [];
% for i = 1:size(e,2)
%     ee=e{i};
%     ratio = sum(ee(1:2:end).^2+ee(2:2:end).^2)^2/((size(e,2)/2)*sum((ee(1:2:end).^2+ee(2:2:end).^2).^2));
%     pariticipation = [pariticipation,ratio];
% end
% figure(4);
% scatter((1:size(e,2))/(size(e,2)),pariticipation);
% xlabel("i/2N");
% ylabel("participation ratio");


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


% successful for i = 2




dlmwrite('eigenv.txt',ee);










