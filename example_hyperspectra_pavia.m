% fcondessa
% script showing the usage of Graph-SegSALSA on hyperspectral data
addpath('data')
addpath('src')
%% load image
% Indian Pine scene
% 10 samples per class with LORSAL-MRL classifier
load('hyperspectral_pavia.mat')
% false color composition of the image
false_color_image;
% class probabilities
probabilities;
% unsupervised segmentations of the HE data
segmentations;
% image of gradients
weight_image;


szimg = size(false_color_image);
[X1,Y1] = gradient(double(mask));
delim = X1|Y1;
white = cat(3,ones(szimg(1),szimg(2)),ones(szimg(1),szimg(2)),ones(szimg(1),szimg(2)));

%% local clusters
% parameter definition
VIS_TAG = 0;
iterations = 20;
mu = 5;
%% GTV only
Z0 = segsalsa(probabilities,mu,iterations,...
    'GTV',1,'tau_gtv',[.2 .5 1 2],'clusters',segmentations,'VIS',VIS_TAG );
[c0,d0] = max(Z0,[],3);
%
figure(123)
imshow(d0,[0 size(probabilities,3)+1]);colormap(colormap_hsv)
hold on;
h= imshow(white);
hold off;
set(h,'AlphaData',delim)
disp('GTV')
mean(d0(mask(:)) == ground_truth(mask(:)))
%title(['GTV only. Acc = ' num2str(mean(d0(mask(:)) == ground_truth(mask(:))))]);
%% VTV only
Z0 = segsalsa(probabilities,mu,iterations,...
    'VTV',1,'tau_vtv',5,'weight_image_vtv',weight_image,'VIS',VIS_TAG );
[c0,d0] = max(Z0,[],3);

figure(124)
imshow(d0,[0 size(probabilities,3)+1]);colormap(colormap_hsv)
hold on;
white = cat(3,ones(szimg(1),szimg(2)),ones(szimg(1),szimg(2)),ones(szimg(1),szimg(2)));
h= imshow(white);
hold off;
set(h,'AlphaData',delim)
disp('VTV')
mean(d0(mask(:)) == ground_truth(mask(:)))
%% GTV and VTV
Z0 = segsalsa(probabilities,mu,iterations,...
    'VTV',1,'tau_vtv',5,'weight_image_vtv',weight_image,...
    'GTV',1,'tau_gtv',[.2 .5 1 2],'clusters',segmentations,'VIS',VIS_TAG);
[c0,d0] = max(Z0,[],3);

figure(125)
imshow(d0,[0 size(probabilities,3)+1]);colormap(colormap_hsv)
hold on;
white = cat(3,ones(szimg(1),szimg(2)),ones(szimg(1),szimg(2)),ones(szimg(1),szimg(2)));
h= imshow(white);
hold off;
set(h,'AlphaData',delim)
disp('GTV + VTV')
mean(d0(mask(:)) == ground_truth(mask(:)))

%% STR only
Z0 = segsalsa(probabilities,mu,iterations,...
    'STR',1,'tau_str',.75,'weight_image_str',weight_image,'gamma_str',1,'window_str',1,'VIS',VIS_TAG );
[c0,d0] = max(Z0,[],3);

figure(126)
imshow(d0,[0 size(probabilities,3)+1]);colormap(colormap_hsv)
hold on;
white = cat(3,ones(szimg(1),szimg(2)),ones(szimg(1),szimg(2)),ones(szimg(1),szimg(2)));
h= imshow(white);
hold off;
set(h,'AlphaData',delim)
disp('STR')
mean(d0(mask(:)) == ground_truth(mask(:)))


%% STR and VTV
Z0 = segsalsa(probabilities,mu,iterations,...
    'STR',1,'tau_str',.75,'weight_image_str',weight_image,'gamma_str',1,'window_str',1,...
    'VTV',1,'tau_vtv',5','weight_image_vtv',weight_image,'VIS',VIS_TAG );
[c0,d0] = max(Z0,[],3);

figure(127)
imshow(d0,[0 size(probabilities,3)+1]);colormap(colormap_hsv)
hold on;
white = cat(3,ones(szimg(1),szimg(2)),ones(szimg(1),szimg(2)),ones(szimg(1),szimg(2)));
h= imshow(white);
hold off;
set(h,'AlphaData',delim)
disp('STR VTV')
mean(d0(mask(:)) == ground_truth(mask(:)))

%% STR and GTV
Z0 = segsalsa(probabilities,mu,iterations,...
    'STR',1,'tau_str',.75,'weight_image_str',weight_image,'gamma_str',1,'window_str',1,...
    'GTV',1,'tau_gtv',[.2 .5 1 2],'clusters',segmentations,'VIS',VIS_TAG );
[c0,d0] = max(Z0,[],3);

figure(128)
imshow(d0,[0 size(probabilities,3)+1]);colormap(colormap_hsv)
hold on;
white = cat(3,ones(szimg(1),szimg(2)),ones(szimg(1),szimg(2)),ones(szimg(1),szimg(2)));
h= imshow(white);
hold off;
set(h,'AlphaData',delim)
disp('STR GTV')
mean(d0(mask(:)) == ground_truth(mask(:)))
