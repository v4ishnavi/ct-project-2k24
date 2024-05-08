function threshold= iterative_threshold(X_, nc)
%here threshold is the output array 
%X is the distribution 
%nc -> number of clusters: here 4 
X = zeros(length(X_), nc);
X(:, 1) = X_;
threshold = zeros(1, nc);
[~, C] = kmeans(X, nc, 'Replicates', 5);
disp(C)
for p = 1:nc
threshold(p) = C(p,1);
end

end