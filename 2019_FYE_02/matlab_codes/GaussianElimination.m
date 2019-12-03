function X = GaussianElimination(A,b)
% perform Gaussian Elimination to solve Ax=b problem


% Jianhong Chen
% 09-21-2019


% useful parameters
[m,~] = size(A);

% perform Gaussian elimination to obtain upper trianglar matrix A
for j = 1:m-1
    for i = j+1:m
        a_ratio = A(i,j)/A(j,j);
        A(i, :) = A(i, :) - A(j, :) * a_ratio;
        b(i) = b(i) - b(j) * a_ratio;
        
    end
end
%

% perform back-substitution to compute X-vector (soln)
% calculate the solution of last row
X(m) = b(m)/A(m, m);
for i = m-1:-1:1 % row iteration start from the bottom
  sum = 0;
  for j = i+1:m % column iteration
    %sum up all the values from the solved terms
    sum(:) = sum(:) + A(i,j) * X(j);
  end 
  X(i) = (b(i) - sum(:))/A(i,i);
end

X = reshape(X,m,1);

end











