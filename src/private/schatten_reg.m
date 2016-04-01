function [V2,V3] = schatten_reg(NU2,NU3,lambda)
V2 = zeros(size(NU2));
V3 = zeros(size(NU3));
for j = 1:size(NU2,2),
        [U,S,V] = svd([NU2(:,j),NU3(:,j)],'econ');
        aux = U*max(abs(S)-lambda(j),0)*V';
        V2(:,j) = aux(:,1);
        V3(:,j) = aux(:,2);
end
        