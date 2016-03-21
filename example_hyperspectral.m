% fcondessa
% script showing the usage of Graph-SegSALSA on hyperspectral data
addpath('data')
addpath('src')
%% load image
% Indian Pine scene
% 10 samples per class with LORSAL-MRL classifier
load('hyperspectral_data.mat')
% false color composition of the image
false_color_image;
% class probabilities
probabilities;
% unsupervised segmentations of the HE data
segmentations;
% image of gradients
weight_image;


%% local clusters
% parameter definition
VIS_TAG = 0;
iterations = 10;
mu = 10;
%% GTV only
Z0 = segsalsa(probabilities,mu,iterations,...
    'GTV',1,'tau_gtv',[1 2 5 10],'clusters',segmentations,'VIS',VIS_TAG );
[c0,d0] = max(Z0,[],3);

figure(123)
imshow(d0,[]);colormap('jet')
title(['GTV only. Acc = ' num2str(mean(d0(mask(:)) == ground_truth(mask(:))))]);
%% VTV only
Z0 = segsalsa(probabilities,mu,iterations,...
    'VTV',1,'tau_vtv',10,'weight_image_vtv',weight_image,'VIS',VIS_TAG );
[c0,d0] = max(Z0,[],3);

figure(124)
imshow(d0,[]);colormap('jet')
title(['VTV only. Acc = ' num2str(mean(d0(mask(:)) == ground_truth(mask(:))))]);
%% GTV and VTV
Z0 = segsalsa(probabilities,mu,iterations,...
    'VTV',1,'tau_vtv',10/2,'weight_image_vtv',weight_image,...
    'GTV',1,'tau_gtv',[1 2 5 10]/2,'clusters',segmentations,'VIS',1);
[c0,d0] = max(Z0,[],3);

figure(125)
imshow(d0,[]);colormap('jet')
title(['GTV and VTV. Acc = ' num2str(mean(d0(mask(:)) == ground_truth(mask(:))))]);

%% STR only
Z0 = segsalsa(probabilities,mu,iterations,...
    'STR',1,'tau_str',0.1,'weight_image_str',weight_image,'gamma_str',1,'window_str',1,'VIS',VIS_TAG );
[c0,d0] = max(Z0,[],3);

figure(126)
imshow(d0,[]);colormap('jet')
title(['STR only. Acc = ' num2str(mean(d0(mask(:)) == ground_truth(mask(:))))]);


%% STR and VTV
Z0 = segsalsa(probabilities,mu,iterations,...
    'STR',1,'tau_str',.1,'weight_image_gtv',weight_image,'gamma_str',1,'window_str',1,...
    'VTV',1,'tau_vtv',1','VIS',VIS_TAG );
[c0,d0] = max(Z0,[],3);

figure(127)
imshow(d0,[]);colormap('jet')
title(['STR and VTV Acc = ' num2str(mean(d0(mask(:)) == ground_truth(mask(:))))]);

%% STR and GTV
Z0 = segsalsa(probabilities,mu,iterations,...
    'STR',1,'tau_str',.1,'weight_image_str',weight_image,'gamma_str',1,'window_str',1,...
    'GTV',1,'tau_gtv',[1 2 5 10]/2,'clusters',segmentations,'VIS',VIS_TAG );
[c0,d0] = max(Z0,[],3);

figure(128)
imshow(d0,[]);colormap('jet')
title(['STR and GTV Acc = ' num2str(mean(d0(:) == ground_truth(:)))]);
