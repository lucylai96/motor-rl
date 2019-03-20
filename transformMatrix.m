function tapNum =transformMatrix(bins)
%this function assigns a bunch of binned taps to a linear tap number

N = max(bins(:));
tapMat = zeros(N);
for i =1:N
    for j = 1:N
        tapMat(i,j) = (i-1)*N+j;
    end
end

for i =1:size(bins,1)
     tapNum(i) =  tapMat(bins(i,1),bins(i,2));
end



end