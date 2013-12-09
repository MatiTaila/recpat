function [Y,S,V,sigmaK,D] = dr_diffusionMap (X,sigmaK,alpha,newDim,options)
%
% function Y = dr_diffusionMap (X,sigmaK,alpha,newDim)
%
% Name: Diffusion Map
%
% Category: Dimensionality Reduction
%
% Description: diffusionKernel function
%
% Input:
%   X: data as dim0 x N matrix (dim0 = dimensionality, N = #points)
%   sigmaK: exponent in the rotation-invariant kernel, if sigmaK==-1, it's
%   estimated from the data.
%     - sigmaK == -1, sigma is estimated from the data
%     - sigmaK == 0 a local scaling is computed (similar to Zelnik-Perona
%     Local Scaling).
%   alpha: parameter of the family of anisotropic diffusion processes
%   newDim: max embedding dimensionality
%   options: various options 
%   * labels: label of each point (if the data is supervised), 1 x N.
%   * distance: kind of distance 
%   options: valid (and necesary arguments are)
%     - distance: which distance use to compute the distance between nodes
%     ('L1', 'L2', 'SDF')
%     - LocalScaling_NN: when sigmaK==0, this is the nearest neighbor to
%     look for compute the local scaling.
%
% Output:
%   Y: embedding as newDim x N matrix
%
% Author: federico lecumberry <fefo at fing dot edu dot uy>. Original code
% written by R. Coifman & S. Lafon; adapted by 
%
% Version: $Id: dr_diffusionMap.m 204 2008-06-16 14:08:08Z fefo $
%

%%% Number of points
N = size(X,2);

%%% If there are some options
if (nargin == 5)
  if (isfield(options,'labels'))
    %%% If there are labels
    m1 = repmat(options.labels,N,1);
    beta = 1;
    gamma = 1;
    sameClass = double(m1 == m1');
    differentClass = 1 - sameClass;
    distanceWeights = beta*sameClass + gamma*differentClass;
  else
    distanceWeights = ones(N,N);
  end  
else
  distanceWeights = ones(N,N);
  options.distance = 'L2';
end

%%% Compute distance between points, default distance is L2
switch options.distance
  case 'L1'
    D = L1Distance(X,X,1);
  case 'L2'
    D = L2Distance(X,X,1);
  case 'SDF'
    D = SDFDistance(X,X,240,120,1);
end

% figure(4500)
% imagesc(D)

%%% Affect the distance matrix with some external/a priori criteria
D = D.*distanceWeights; % afectar D con los pesos propuestos

%%% With sigmaK==-1, sigma is estimated from the data
if (sigmaK == -1)
  D1 = sort(D); % ordenar las distancias entre los puntos
  Dmin = max(D1(3,:)); % max. distancia al segundo vecino más cercano
  Dmax = max(D1(7,:)); % max. distancia al quinto vecino más cercano
  pres = 0.3; %
  sigma_min = Dmin/sqrt(log(1/pres));
  sigma_max = Dmax/sqrt(log(1/pres));
  sigmaK = 0.5*(sigma_min + sigma_max);
  %sigmaK = mean(mean(D1(3:5,:)));
  disp(['Sigma: ' num2str(sigmaK)])
  K = exp(-(D/sigmaK).^2);
elseif (sigmaK == 0)
  %%% If sigmaK==0 a local scaling is computed (similar to Zelnik-Perona
  %%% Local Scaling). K is the rotation-invariant kernel
  K = LocalScaling(D,options.LocalScaling_NN);
else
  K = exp(-(D/sigmaK).^2);
end

p = sum(K);
p = p(:);
K1 = K./((p*p').^alpha);
v = sqrt(sum(K1));
v = v(:);
A = K1./(v*v');

if (sigmaK >= 0.5)
  %%% Sparse version
  thre = 1e-7;
  M = max(max(A));
  A = sparse(A .* double(A > thre*M) );
  [U,S,V] = svds(A,newDim+1);
  U = U./(U(:,1)*ones(1,newDim+1));
else
  %%% Full version
  [U,S,V] = svd(A,0); 
  U = U./(U(:,1)*ones(1,size(U,1)));
end
Y = U(:,2:newDim+1);
Y = Y';

% figure(343897623)
% plot(diag(S),'r.-')
% title('Eigenvalues')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function AN = LocalScaling(D,K)
%%% Sort the distance matrix to find the distance to the K/2 neighbor
Dsort = sort(D);

%%% Take the local_dist as the median over the distances from each point to
%%% its neighbors
local_dist = Dsort(floor((K+1)/2)+1,:)';
N = local_dist*local_dist';
DN = (D.^2) ./ N;            %% Locally scaled distance matrix
AN = exp(-DN);               %% Locally scaled affinity matrix
AN = AN .* ~eye(size(D,1));  %% Zero out the diagonal
