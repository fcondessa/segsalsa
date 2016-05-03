%multi_segs_show


colors_segms = [0  0 1; 0 1 1; 1 0 1; 1 1 0];

szimg = size(image);

figure(119);
imshow(image);
for lev = 1:size(segmentations,3);


hold on;
segm_col = cat(3,colors_segms(lev,1)*ones(szimg(1),szimg(2)),...
    colors_segms(lev,2)*ones(szimg(1),szimg(2)),...
    colors_segms(lev,3)*ones(szimg(1),szimg(2)));
h= imshow(segm_col);
[Xa,Ya] = gradient(double(segmentations(:,:,lev)));
segm_map = 0.2*imdilate(Xa | Ya,ones(4)) + 0.6*(Xa | Ya);
hold off;
set(h,'AlphaData',segm_map)

end