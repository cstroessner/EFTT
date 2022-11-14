function Ak = MatrixUnfolding(A,k,rowpoints,columnpoints)
% A function handle: m * d -> m * 1
% m number of points, d dimension

% MatrixUnfolding: ( m * 1, n * (d-1) ) -> m * n
m = size(rowpoints,1);
n = size(columnpoints,1);
d = size(columnpoints,2) + 1;
Ak = zeros(m,n);

temp = repmat(columnpoints,m,1);
temp = A([temp(:,1:k-1),reshape(repmat(rowpoints,1,n)',[],1),temp(:,k:(d-1))]);
Ak = reshape(temp,n,m)';

%Ak = bsxfun(@(row,col) A([col(1,1:k-1),row,col(1,k:(d-1))]),rowpoints,columnpoints);


% if m <= n
%     for i = 1:m
%         Ak(i,:) = A([columnpoints(:,1:k-1),repmat(rowpoints(i),n,1),columnpoints(:,k:(d-1))]);
%     end
% else
%     for i = 1:n
%         temp = repmat(columnpoints(i,:),m,1);
%         Ak(:,i) = A([temp(:,1:k-1),rowpoints,temp(:,k:(d-1))]);
%     end
% end

end