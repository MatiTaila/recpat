function d = L1Distance(a,b,df)
% function d = L1Distance(a,b,df)
%
% Name: L1Distance
%
% Category: Vector Distance
%
% Description: Compute the L1 distance between all the D-dimensional vetor
%
% Input:
%    a: (DxM) matrix 
%    b: (DxN) matrix
%    df: 1, force diagonals to be zero; 0 (default), do not force
%
% Output:
%    d - (MxN) L1 distances between vectors in A and B
%
% Author: federico lecumberry <fefo at fing dot edu dot uy>
%
% Version: $Id: L1Distance.m 70 2007-08-22 22:00:48Z fefo $
% 

%%% Check enough elements
if (nargin < 2)
  error('Not enough input arguments');
end

%%% Force 0 in the diagonal
if (nargin < 3)
  df = 0; % By default, do not force 0 on the diagonal
end

%%% Check the dimensionality
if (size(a,1) ~= size(b,1))
  error('A and B should be of same dimensionality');
end

%%% Check for real elements
if ~(isreal(a)&&isreal(b))
  disp('Warning: running with imaginary numbers.  Results may be off.'); 
end

%%% Run the for loop for the matrix with less elements. For each column
%%% compute the L1 distance with all the column in the other matrix.
if (size(b,2) < size(a,2))
  for i = 1:size(b,2)
    tb = repmat(b(:,i), [1 size(a,2)]);
    d(:,i) = sum(abs(a - tb),1)'; % Sum along the columns
  end
else
  for i = 1:size(a,2)
    ta = repmat(a(:,i), [1 size(b,2)]);
    d(i,:) = sum(abs(b - ta),1); % Sum along the columns
  end
end

%%% Make sure result is all real
d = real(d); 

%%% Force zero on the diagonal? 
if (df == 1)
  d = d.*(1-eye(size(d)));
end

