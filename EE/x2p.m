function [pp,beta] = x2p(xx,kk,tol)
%[pp,beta] = x2p(xx,kk,tol)
%
% convert raw input data xx into stochastic neighbourhood
% probabilities pp, using effective number of neighbours kk
% around each point
% (does a binary search in temperature to achieve this, tol is
% the convergence criterion for that search)
%
% Based on the original version of Sam Roweis.
%
% History:
%   05/27/12, Max Vladymyrov: vectorize square distnce computation.

more off;

if(nargin<3) tol=1e-4; end %tol=1e-4*log(kk)

N=size(xx,2);
pp=zeros(N,N);
beta = zeros(N,1);
lkk=log(kk);

sumxx = sum(xx .^ 2, 1);
D = bsxfun(@plus, sumxx', bsxfun(@plus, sumxx, -2*xx'*xx));

for nn=1:N
  thisdd=D(nn, [1:nn-1 nn+1:end]);
  
  betamin=-Inf; betamax=Inf;
  beta(nn)=1; %???
  [H,thispp]=Hbeta(thisdd,beta(nn));
  hdiff = H-lkk;
  while(abs(hdiff)>tol)
    if(hdiff>0)
      betamin=beta(nn);
      if(isinf(betamax))
        beta(nn)=beta(nn)*2;
      else
        beta(nn)=(beta(nn)+betamax)/2;
      end
    else
      betamax=beta(nn);
      if(isinf(betamin))
        beta(nn)=beta(nn)/2;
      else
        beta(nn)=(beta(nn)+betamin)/2;
      end
    end
    [H,thispp] = Hbeta(thisdd,beta(nn));
    hdiff = H-lkk;
  end
  
  %fprintf(1,'Point %d/%d beta=%g\r',nn,N,beta(nn));
  pp(nn,[1:(nn-1),(nn+1):end])=thispp;
  
end

function [H,ppn] = Hbeta(ddn,beta)
%  ff = exp(-ddn*beta);
shift = log(realmax)/2 - max(-ddn*beta);
ff = exp(-ddn*beta + shift);
zz = sum(ff);
ppn=ff/zz;
%  H = log(zz)+beta*sum(ddn.*ff)/zz;
H = log(zz)-shift+beta*sum(ddn.*ppn);
