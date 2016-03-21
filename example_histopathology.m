% fcondessa
% 2016.03.20
% script showing the usage of Graph-SegSALSA on histopathology data
addpath('data')
addpath('src')
%% load image
load('histopathology_data.mat')
% struct containing HE data
imageData;
% class probabilities
probabilities;
% unsupervised segmentations of the HE data
segmentations;
% image of gradients
weight_image;
%%
% subsampling of the original image to account for the patch-based
% probabilities
step = 4;
img = imageData.originalImage(1:step:end,1:step:end,:);
[~,~,ground_truth] = unique(imageData.roiMask(1:step:end,1:step:end,:));
ground_truth = reshape(ground_truth,size(imageData.roiMask(1:step:end,1:step:end,:)));

%% local clusters
% parameter definition
iterations = 5;
mu = 10;
%% GTV only
Z0 = segsalsa(probabilities,mu,iterations,...
    'GTV',1,'tau_gtv',[1 2 5 10],'clusters',segmentations,'VIS',0);
[c0,d0] = max(Z0,[],3);

figure(123)
imshow(d0,[]);colormap('jet')
title(['GTV only. Acc = ' num2str(mean(d0(:) == ground_truth(:)))]);
%% VTV only
Z0 = segsalsa(probabilities,mu,iterations,...
    'VTV',1,'tau_vtv',10,'weight_image_vtv',weight_image,'VIS',0);
[c0,d0] = max(Z0,[],3);

figure(124)
imshow(d0,[]);colormap('jet')
title(['VTV only. Acc = ' num2str(mean(d0(:) == ground_truth(:)))]);
%% GTV and VTV
Z0 = segsalsa(probabilities,mu,iterations,...
    'VTV',1,'tau_vtv',10/2,'weight_image_vtv',weight_image,...
    'GTV',1,'tau_gtv',[1 2 5 10]/2,'clusters',segmentations,'VIS',0);
[c0,d0] = max(Z0,[],3);

figure(125)
imshow(d0,[]);colormap('jet')
title(['GTV and VTV. Acc = ' num2str(mean(d0(:) == ground_truth(:)))]);

%% STR only
Z0 = segsalsa(probabilities,mu,iterations,...
    'STR',1,'tau_str',1,'weight_image_str',weight_image,'gamma_str',1,'window_str',2,'VIS',0);
[c0,d0] = max(Z0,[],3);

figure(126)
imshow(d0,[]);colormap('jet')
title(['STR only. Acc = ' num2str(mean(d0(:) == ground_truth(:)))]);


%% STR and VTV
Z0 = segsalsa(probabilities,mu,iterations,...
    'STR',1,'tau_str',1,'weight_image_gtv',weight_image,'gamma_str',1,'window_str',2,...
    'VTV',1,'tau_vtv',1','VIS',0);
[c0,d0] = max(Z0,[],3);

figure(127)
imshow(d0,[]);colormap('jet')
title(['STR and VTV Acc = ' num2str(mean(d0(:) == ground_truth(:)))]);

%% STR and GTV
Z0 = segsalsa(probabilities,mu,iterations,...
    'STR',1,'tau_str',1,'weight_image_str',weight_image,'gamma_str',1,'window_str',2,...
    'GTV',1,'tau_gtv',[1 2 5 10]/2,'clusters',segmentations,'VIS',0);
[c0,d0] = max(Z0,[],3);

figure(128)
imshow(d0,[]);colormap('jet')
title(['STR and GTV Acc = ' num2str(mean(d0(:) == ground_truth(:)))]);
