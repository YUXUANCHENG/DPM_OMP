clc
clear
close all

H_matrix_file = "C:\Users\Yuxuan Cheng\OneDrive\Windows_Version\remove\H_matrix.txt";
Hmatrix = dlmread(H_matrix_file);

M_matrix_file = "C:\Users\Yuxuan Cheng\OneDrive\Windows_Version\remove\M_matrix.txt";
Mmatrix = dlmread(M_matrix_file);

H_matrix_file_remove = "C:\Users\Yuxuan Cheng\OneDrive\Windows_Version\remove\H_matrix_remove.txt";
Hmatrix_remove = dlmread(H_matrix_file_remove);

M_matrix_file_remove = "C:\Users\Yuxuan Cheng\OneDrive\Windows_Version\remove\M_matrix_remove.txt";
Mmatrix_remove = dlmread(M_matrix_file_remove);

coordinate_file = "C:\Users\Yuxuan Cheng\OneDrive\Windows_Version\remove\coordinate_file_1.txt";
radius_file = "C:\Users\Yuxuan Cheng\OneDrive\Windows_Version\remove\radius_1.txt";
contact_file = "C:\Users\Yuxuan Cheng\OneDrive\Windows_Version\remove\contact_1.txt";

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

for i = 1:5
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
contact = csvread(contact_file);

N=64;
frames= size(coordinate,1)/N ;


if (N * (N-1) /2 ~= size(contact,2))
    print("stop");
end

for i = 1 : frames
    start_point = 1 + N * ( i - 1 );
    end_point = N * i;
    plot_particles_2d(1,[1,1],radius(start_point:end_point),coordinate(start_point:end_point,1),coordinate(start_point:end_point,2))
end

ee=e{4};
quiver(coordinate(:,1),coordinate(:,2),ee(1:2:end),ee(2:2:end),'color','r','linewidth',3);


eigenValues_remove(eigenValues_remove< 10^-7) = 0;
figure(2);
[counts, binCenters] = hist(sqrt(eigenValues_remove),10);
counts = counts / sum(counts);
plot(binCenters, counts,'linewidth',1.75);
xlabel("w");
ylabel("D(w)");
figure(3);
scatter((1:2*N-8)/(2*N-8),sqrt(eigenValues_remove));
xlabel("i/2N");
ylabel("w");

e_r={};

for i = 1:size(eigenVectors_remove,1)
    e_r{i} = eigenVectors_remove(:,i)/norm(eigenVectors_remove(:,i));
end

len = size(eigenVectors_remove,1);
x_tran = zeros(1,len);
x_tran(1:2:end)=1;
x_tran = x_tran'/norm(x_tran);

y_tran = zeros(1,len);
y_tran(2:2:end)=1;
y_tran = y_tran'/norm(y_tran);
e_r{1} = x_tran;
e_r{2} = y_tran;

e_r = GramSchmidt(e_r);

pariticipation = [];
for i = 1:size(e_r,2)
    ee=e_r{i};
    ratio = sum(ee(1:2:end).^2+ee(2:2:end).^2)^2/((size(e_r,2)/2)*sum((ee(1:2:end).^2+ee(2:2:end).^2).^2));
    pariticipation = [pariticipation,ratio];
end
figure(4);
scatter((1:size(e_r,2))/(size(e_r,2)),pariticipation);
xlabel("i/2N");
ylabel("participation ratio");


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















