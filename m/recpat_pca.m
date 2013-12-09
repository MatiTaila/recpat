function [y, V] = recpat_pca(x,thr)

% x	     = x - repmat(mean(x,2),1,size(x,2));
sigma    = cov(x');
[V, D]   = eig(sigma);
D		 = diag(D); 
[~, ind] = sort(D,'descend');
V		 = V(:,ind);
D		 = D(ind);

if nargin > 1
	ind  = find(cumsum(D)>=thr*sum(D),1,'first');
	V    = V(:,1:ind);
end

y = V'*x;