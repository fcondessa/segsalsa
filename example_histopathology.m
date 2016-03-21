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
iterations = 50;
mu = 10;
%% GTV only
Z0 = segsalsa(probabilities,mu,iterations,...
    'GTV',1,'tau_gtv',[1 2 5 10],'clusters',segmentations);
[c0,d0] = max(Z0,[],3);

figure(123)
imshow(d0,[]);colormap('jet')
title(['GTV only. Acc = ' num2str(mean(d0(:) == ground_truth(:)))]);
%% VTV only
Z0 = segsalsa(probabilities,mu,iterations,...
    'VTV',1,'tau_vtv',10,'weight_image',weight_image);
[c0,d0] = max(Z0,[],3);

figure(124)
imshow(d0,[]);colormap('jet')
title(['GTV only. Acc = ' num2str(mean(d0(:) == ground_truth(:)))]);
%% GTV and VTV
Z0 = segsalsa(probabilities,mu,iterations,...
    'VTV',1,'tau_vtv',10/2,'weight_image',weight_image,...
    'GTV',1,'tau_gtv',[1 2 5 10]/2,'clusters',segmentations);
[c0,d0] = max(Z0,[],3);

figure(125)
imshow(d0,[]);colormap('jet')
title(['GTV only. Acc = ' num2str(mean(d0(:) == ground_truth(:)))]);