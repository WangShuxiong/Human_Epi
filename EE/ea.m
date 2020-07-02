% [W,s] = ea(X,K[,k]) Gaussian entropic affinities
%
% This computes the matrix W of Gaussian entropic affinities (EAs) and the
% corresponding bandwidth values s for a dataset X and a desired perplexity
% K. By default, W is a full matrix. Use the optional argument k (with k > K)
% to restrict the EAs to the k nearest neighbors of each point, in which case
% W will be sparse. Note that, if using k neighbors, the maximum perplexity
% achievable is K=k, in the limit where s=Inf (which results in affinities
% equal to 1/k for each neighbor), so k < K makes no sense. In practice, K
% should not approach k too much.
%
% W is normalized so each row sums 1, thus it is a stochastic matrix. It is
% the random-walk matrix commonly used in machine learning, but where each
% row (data point) has its own bandwidth in the Gaussian kernel.
%
% When using sparse affinities, the runtime is mostly due to computing the
% nearest neighbors in nnsqdist.m, which is O(D.N²), rather than to computing
% the EAs, which is O(N). If you have precomputed the nearest neighbors, you
% should pass them as input argument to ea.m. See also nnsqdist.m about using
% selection rather than sorting.
%
% In:
%   X: NxD matrix of N row D-dimensional vectors.
%   K: scalar in (0,N), the perplexity.
%   k: either a scalar k > K, the number of neighbors; or a cell array {D2,nn}
%      containing two N x k matrices, where D2(n,i) is the square distance to
%      the ith nearest neighbor and nn(n,i) is the index of the ith nearest
%      neighbor. Default: N-1.
% Out:
%   W: NxN matrix of EAs (random-walk matrix).
%   s: Nx1 vector of bandwidth values.
%
% Any non-mandatory argument can be given the value [] to force it to take
% its default value.

% Copyright (c) 2013 by Max Vladymyrov and Miguel A. Carreira-Perpinan

function [W,s] = ea(X,K,k)

N = size(X,1);
% ---------- Argument defaults ----------
if ~exist('k','var') || isempty(k) k = N-1; end;
if iscell(k) 
  D2 = k{1}; nn = k{2}; k = size(D2,2);
else
  [D2,nn] = nnsqdist(X,k);		% Square distances to k nn
  %[D2,nn] = nnsqdist(X,k,'mink');	% Usually faster, see nnsqdist.m
end
% ---------- End of "argument defaults" ----------

b = zeros(N,1); Wp = zeros(N,k); logK = log(K);
[B,D2] = eabounds(logK,D2);		% Log-beta bounds
[~,p] = sort(D2(:,ceil(K)));		% Point order: distance to Kth nn
j = p(1); b0 = mean(B(j,:)); p=[p;0];	% Initialization
for i=1:N				% Compute log-beta & EAs for each point
  [b(j),Wp(j,:)] = eabeta(D2(j,:),b0,logK,B(j,:));
  b0 = b(j); j = p(i+1);		% Next point
end
W = sparse(repmat((1:N)',1,k),nn,Wp,N,N); if k>=N-1 W = full(W); end;
if nargout==2 s=1./sqrt(2*exp(b)); end	% Bandwidths from log-beta values


% [B,D2] = eabounds(logK,D2) Gaussian EAs: bounds
%
% Computes simple bounds (in constant time) of the beta values given the
% distances of each point to its k nearest neighbors, contained in matrix D2:
%   D2(n,i) = squared distance from point X(n,:) to its ith nearest neighbor.
%
% We return D2 because the bounds' formula needs strictly increasing neighbor
% distances, so in case of tied distances we perturb them a little. Tied
% distances can happen with quantized data values, e.g. if the dataset is an
% image that has areas with constant color.
%
% In:
%   logK: scalar, log of the perplexity. 
%   D2: N x k matrix of sorted square distances to the k nearest neighbors.
% Out:
%   B: Nx2 matrix of log-beta bounds for each data point.
%   D2: same as D2 with a tiny perturbation of the distance to the first
%      nearest neighbor if K or more nearest neighbors are at exactly the
%      same distance.

% Copyright (c) 2013 by Max Vladymyrov and Miguel A. Carreira-Perpinan

function [B,D2] = eabounds(logK,D2)

N = size(D2,2);	% We call "N" the number of nearest neighbors, as in the paper
logN = log(N); logNK = logN-logK;

delta2 = D2(:,2)-D2(:,1);
% Ensure delta2 >= eps
ind = find(delta2<eps); i=3;
flag = 1;
while ~isempty(ind)
  if i>exp(logK) && flag
    % Permute ever so slightly the distance to the first neighbor for points
    % that have K or more closest neighbors at the same distance. Practically 
    % this happens only when K is very small.
    D2(ind,1) = D2(ind,1)*0.99; flag = 0;
  end
  delta2(ind) = D2(ind,i)-D2(ind,1); ind = find(delta2<eps); i = i+1;
end

deltaN = D2(:,N)-D2(:,1);

% Compute p1(N,logK)
if logK > log(sqrt(2*N))
  p1 = 3/4;
else
  p1 = 1/4;
  for i=1:100 e = -p1*log(p1/N)-logK; g = -log(p1/N)+1; p1 = p1-e/g; end
  p1 = 1-p1/2;
end

bU1 = (2*log(p1*(N-1)/(1-p1)))./delta2;
bL1 = (2*logNK/(1-1/N))./deltaN;
bL2 = (2*sqrt(logNK))./sqrt(D2(:,N).^2-D2(:,1).^2);
B = log([max(bL1,bL2) bU1]);


% [b,W] = eabeta(d2,b0,logK,B) Gaussian EAs: beta and affinities
%
% Computes the values of beta and the corresponding Gaussian affinities for
% one point. It does root finding with Newton's method embedded in a bisection
% loop to ensure global convergence.
%
% In:
%   d2: 1 x k vector of square distances to the k nearest neighbors.
%   b0: initial value of beta.
%   logK: log of the perplexity K. 
%   B: 1x2 vector of lower and upper bounds on log(beta).
% Out:
%   b: log(beta) value.
%   W: 1 x k vector of EAs computed with this beta. 

% Copyright (c) 2013 by Max Vladymyrov and Miguel A. Carreira-Perpinan

function [b,W] = eabeta(d2,b0,logK,B)

% No need to change these ever, probably:
maxit = 20;		% Max. no. iterations before a bisection
tol = 1e-10;		% Convergence tolerance to stop iterating

if (b0<B(1) || b0>B(2)) b=(B(1)+B(2))/2; else b=b0; end

i = 1;			% maxit counter
% Inside the loop, pbm is a flag to detect any numerical problems (zero
% gradient, infinite function value, etc.).
while 1
  bE = exp(b); pbm = 0;
  
  % Compute the function value: m0, m1v, m1 are O(N)
  ed2 = exp(-d2*bE); m0 = sum(ed2);
  if m0<realmin		% Numerical error
    e = -logK; pbm = 1;
  else
    m1v = ed2.*d2/m0; m1 = sum(m1v); e = bE*m1 + log(m0) - logK;
  end
  
  if abs(e) < tol break; end
  
  % Very narrow bounds, no need to iterate. This can happen if K is very small.
  if B(2)-B(1) < 10*eps break; end
  
  % Update the bounds
  if (e<0 && b<=B(2))
    B(2) = b;
  elseif (e>0 && b>=B(1))
    B(1) = b;
  end
  
  pbm = pbm || isinf(e) || e<-logK || e>log(length(d2))-logK;
  
  if ~pbm
    if i==maxit		% Exceeded maxit, bisection step
      b = (B(1)+B(2))/2; i=1; continue;
    end
    % Compute the gradient: m2 is O(N)
    eg2 = bE^2; m2 = m1v*d2'; m12 = m1^2-m2; g = eg2*m12;
    if g==0 pbm=1; end
  end
  
  % If there was a problem with the function, do bisection with old bounds.
  % If there was a problem with the gradient, do bisection with new bounds.
  if pbm 
    % If there is a numerical problem on both ends of the bounds, return 
    % current value. Practically this happens only when K is very small.
    esqd1 = exp(-d2*exp(B(1))); esqd2 = exp(-d2*exp(B(2)));
    if sum(esqd1+esqd2) < 2*sqrt(realmin) break; end
    b = (B(1)+B(2))/2; i=1; continue; 
  end
  
  % Newton step ok, update bounds
  p = -e/g; b = b + p;
  
  if (b<B(1) || b>B(2))	% Out of bounds, bisection step
    b = (B(1)+B(2))/2; i = 0;
  end
  i=i+1;
end

W = ed2/m0;		% Affinities


% [D2,nn] = nnsqdist(X,k[,method]) Nearest-neighbor squared distances
%
% This computes the k nearest neighbors of each data point in X and its
% squared Euclidean distances:
% - D2(n,i) = squared distance from point X(n,:) to its ith nearest neighbor.
% - nn(n,i) = index of the ith nearest neighbor of X(n,:).
%
% Finding the nearest neighbors is done by default by sorting N distances for
% each point. If k<<N, a faster way is to use selection in O(N+k.logk) rather
% than sorting in O(N.logN). To use selection, install this package and set
% method='mink':
%   Min/Max selection
%   http://www.mathworks.com/matlabcentral/fileexchange/23576-minmax-selection
% Its mink() function selects the k smallest elements in the array with a
% partial quicksort in O(N+k.logk).
% Note that computing the distances between pairs of points is O(D.N²) anyway,
% so the improvement of selection is larger the smaller D is.
%
% In:
%   X: NxD matrix of N row D-dimensional vectors.
%   k: number of nearest neighbors.
%   method: find nearest neighbors with sorting ('sort') or selection ('mink').
%      Default: 'sort'.
% Out:
%   D2: N x k matrix of sorted square distances to the k nearest neighbors.
%   nn: N x k matrix of indices of the corresponding nearest neighbors.
%
% Any non-mandatory argument can be given the value [] to force it to take
% its default value.

% Copyright (c) 2013 by Max Vladymyrov and Miguel A. Carreira-Perpinan

function [D2,nn] = nnsqdist(X,k,method)

% ---------- Argument defaults ----------
if ~exist('method','var') || isempty(method) method='sort'; end
% ---------- End of "argument defaults" ----------

N = size(X,1); k = min(N-1,k);
% Block size (see below), select as large as your memory allows
mem = 1; B = floor((mem*1024^3)/(4*N*8)/2);	% This will fit in mem GB RAM

% Process by blocks to save memory
i1 = 1; i2 = min(N,B);
Xt = X'; X2 = X.^2; x2 = sum(X2,2)'; D2 = zeros(N,k); nn = D2;

while i1 <= N
  % This computes squared distances in a fast, vectorized way but can have
  % cancellation error for points closer than sqrt(eps).
  sd1 = max(bsxfun(@plus,sum(X2(i1:i2,:),2),...
                   bsxfun(@minus,x2,2*X(i1:i2,:)*Xt)),0);
  if method == 'sort'
    [sd1,ind] = sort(sd1,2);
  else
    [sd1,ind] = mink(sd1,k+1,2);
  end
  D2(i1:i2,:) = sd1(:,2:(k+1)); nn(i1:i2,:) = ind(:,2:(k+1));
  i1 = i1 + B; i2 = min(N,i1+B-1);
end

