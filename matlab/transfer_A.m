function A = transfer_A(N, contact)

A = zeros(N);
for i = 1: N-1
    for j = (i+1):N
        A(i,j) = contact((i-1)*N - (i*(i+1))/2 + j);
    end
end

end