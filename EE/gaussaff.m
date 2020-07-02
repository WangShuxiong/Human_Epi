% [W,nn] = gaussaff(X,G[,s,a,sqd]) Gaussian affinity matrix
%
% Given a graph with nodes at X, it assigns to edge i->j a weight
% Wij = exp(-|Xi-Xj|²/(2.si.sj)), or 0 if such edge doesn't exist.
% It can use the following Euclidean-distance neighbourhood graphs:
% - e-ball (symmetric): i~j if dist(i,j) <= e.
% - k-nearest-neighbour:
%   . nonsymmetric: i->j if j is among the K nearest neighbours of i
%   . symmetric: i~j if j is among the K nearest neighbours of i or vice versa
%   . mutual: i~j if j is among the K nearest neighbours of i and vice versa.
% If a is given, W is normalised as DD^(-a).W.DD^(-a) with DD=diag(sum(W,2)).
% W is sparse unless the graph is fully connected (large K or e).
%
% In:
%   X: NxD data set of row vectors.
%   G: cell array, the graph type:
%      {'k',K}: symmetric k-nn graph with K neighbours;
%      {'K',K}: nonsymmetric k-nn graph with K neighbours;
%      {'m',K}: mutual k-nn graph with K neighbours;
%      {'e',e}: e-ball graph with distance e.
%   s: scalar, the Gaussian width; or a Nx1 array of widths (one per point).
%      Default: Inf (ie, binary edges).
%   a: scalar in [0,1], to normalise W as in diffusion maps. Default: 0.
%   sqd: NxN nonsparse matrix of X's squared distances. Default: compute them.
% Out:
%   W: the NxN weight matrix, possibly sparse.
%   nn: Nx(K+1) array of indices so that nn(n,k+1) points to the kth nearest
%      neighbour of X(n,:). If the graph is fully connected or of type 'e'
%      we return nn=[].
%
% Any non-mandatory argument can be given the value [] to force it to take
% its default value.

% Copyright (c) 2011 by Miguel A. Carreira-Perpinan

function [W,nn] = gaussaff(X,G,s,a,sqd)

% ---------- Argument defaults ----------
if nargout>1 nn = []; end;
if ~exist('s','var') | isempty(s) s = Inf; end;
if ~exist('sqd','var') sqd = []; end;
if ~isscalar(s) invs = sparse(1./s(:)); end;
[N,D] = size(X);
if (s==0) | (ismember(G{1},{'k','K','m','e'}) & G{2}<=0)
  W = speye(N,N); if nargout>1 nn = (1:N)'; end; return;
end
% ---------- End of "argument defaults" ----------
% Block size (see below), select as large as your memory allows
mem = 1; B = floor((mem*1024^3)/(4*N*8)/2);	% This will fit in mem GB RAM

if isinf(G{2}) | (ismember(G{1},{'k','K','m'}) & (G{2}>=N-1)) % Full graph
  if isinf(s)
    W = ones(N,N);
  else
    if isempty(sqd)
      if isscalar(s)
        W = exp(-sqdist(X)/(2*s^2));
      else
        W = exp(-(diag(invs)*sqdist(X)*diag(invs))/2);
      end
    else
      if isscalar(s)
        W = exp(-sqd/(2*s^2));
      else
        W = exp(-(diag(invs)*sqd*diag(invs))/2);
      end
    end
  end
else	% Sparse graph (k-nn or e-ball)
  K = G{2}; e = G{2};
  
  % Find neighbours and distances
  if isempty(sqd)
    % Process by blocks to save memory
    i1 = 1; i2 = min(N,B); ind1 = []; ind2 = []; sd = [];
    X2 = X.^2; x2 = sum(X2,2)';
    while i1 <= N
      sd1 = max(bsxfun(@plus,sum(X2(i1:i2,:),2),...
                       bsxfun(@minus,x2,2*X(i1:i2,:)*X')),0);
      if G{1}=='e'
        % sd1 = sqdist(X(i1:i2,:),X);
        % Stupid Matlab screws this for B=1
        [Ind1,Ind2] = find(sd1 <= e^2);
        ind1 = [ind1;i1-1+Ind1]; ind2 = [ind2;Ind2];
        sd = [sd;sd1(sub2ind(size(sd1),Ind1,Ind2))];
      else
        % [sd1,ind] = sort(sqdist(X(i1:i2,:),X),2);
        [sd1,ind] = sort(sd1,2);
        ind1 = [ind1;repmat((i1:i2)',1,K+1)]; ind2 = [ind2;ind(:,1:(K+1))];
        sd = [sd;sd1(:,1:(K+1))];
      end
      i1 = i1 + B; i2 = min(N,i1+B-1);
    end
  else
    if G{1}=='e'
      [ind1,ind2] = find(sqd <= e^2); sd = sqd(sub2ind([N N],ind1,ind2));
    else      
      [sd,ind] = sort(sqd,2);
      ind1 = repmat((1:N)',1,K+1); ind2 = ind(:,1:(K+1)); sd = sd(:,1:(K+1));
    end
  end
  if (nargout>1) & (G{1}~='e') nn = ind2; end;
  
  % Correct neighbours depending on type of k-nn graph
  if G{1}~='e'
    ind1 = ind1(:); ind2 = ind2(:); sd = sd(:);
    if G{1}~='K'	% 'k' and 'm' are symmetric graphs
      A = sparse(ind1,ind2,sd,N,N); [I,J] = find(A<A');	% Nonsymmetric edges
      ind = find(A(sub2ind([N N],I,J))==0); I = I(ind); J = J(ind);
      if G{1}=='k'	% symmetric k-nn: A OR A' -> add missing edges [I J]
        ind1 = [ind1;I]; ind2 = [ind2;J]; sd = [sd;A(sub2ind([N N],J,I))];
      else		% mutual k-nn: A AND A' -> remove surplus edges [J I]
        ind = setdiff(sub2ind([N N],ind1,ind2),sub2ind([N N],J,I));
        [ind1,ind2] = ind2sub([N N],ind); sd = A(ind);
      end
    end
  end
  
  % Compute Gaussian affinity for each edge
  if isinf(s)
    W = sparse(ind1,ind2,ones(length(ind1),1),N,N);
  else
    if isscalar(s)
      W = sparse(ind1,ind2,exp(-sd/(2*s^2)),N,N);
    else
      W = sparse(ind1,ind2,exp(-(invs(ind1).*sd.*invs(ind2))/2),N,N);
    end
    if G{1}~='K' W = (W+W')/2; end;	% Ensure symmetry
  end
  
end

if exist('a','var') & ~isempty(a)
  DD = diag(sparse(sum(W,2).^(-a))); W = DD*W*DD;
end;

