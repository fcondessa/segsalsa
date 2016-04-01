% overimposeimage
% h1 = overimposeimage(image,labelmap,scrible,alpha_val,colo)
%   image is a imshowable image (uint8 if needed)
%   labelmap is a map of labels (1 through K)
%   scrible is a map of scrible locations (-1 -> no scrible; 0 -> scrible
%   of class 1; ...; K-1 -> scrible of class K
%   alpha_val is the value of the nonscrible transparency
%   colo is colormap of the labelmap
% FJCC 2015.12.03
function h1 = overimposeimage(image,labelmap,scrible,alpha_val,colo)
alpha_map = alpha_val * ones(size(image,1),size(image,2));
enlarge_scrible = 2;

% IF scrible goes from -1 to K-1
scrib_mask = imdilate(scrible ~= -1,ones(enlarge_scrible));
labelmap = labelmap .*(1-scrib_mask) + (1+imdilate(scrible,ones(enlarge_scrible))).*scrib_mask;
% IF scrible goes from 0 to K
%scrib_mask = imdilate(scrible ~= 0,ones(enlarge_scrible));
%labelmap = labelmap .*(1-scrib_mask) + (imdilate(scrible,ones(enlarge_scrible))).*scrib_mask;

alpha_map(scrib_mask) = 1;
h1 = imshow(image);
hold on;
h= imshow(labelmap,[]);
colormap(colo)
hold off;
set(h,'AlphaData',alpha_map)
