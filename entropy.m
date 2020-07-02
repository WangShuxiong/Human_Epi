function z = entropy(P)
% P: kxn probability matrix


[~,n] = size(P);
z = zeros(n,1);

% for i = 1:n
%     z(i) = wentropy(P(:,i),'shannon'); % 'log energy','shannon'
% end

for i = 1:n
    zzv = zeros(size(P(:,i)));
    zz = find(P(:,i) == 0);
    zzv(zz) = 1;
    z(i) = -sum(P(:,i).*log2(P(:,i) +zzv ));
end
