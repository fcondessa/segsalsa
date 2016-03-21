% patchbased Jacobian
function [outx,outy] = patch_jacobian(Z,nr,nc,no_classes,Lk,gamma)
[X,Y] = meshgrid(-Lk:Lk);
K = exp(-(X.^2+Y.^2)/gamma);
L = length(K(:));
outx = zeros(nr,nc,no_classes*L);
outy = zeros(nr,nc,no_classes*L);
[prex,prey] = gradient(reshape(Z',[nr,nc,no_classes]));
for j = 1:L;
    for i = 1:no_classes,
        outx(:,:,(j-1)*no_classes + i) = K(j) * circshift(prex(:,:,i),[X(j),Y(j)]);
        outy(:,:,(j-1)*no_classes + i) = K(j) * circshift(prey(:,:,i),[X(j),Y(j)]);
    end
end
