% fcondessa
% script showing the usage of Graph-SegSALSA on hyperspectral data
addpath('data')
addpath('src')
addpath('ext/SLIC_mex')
%addpath('/Users/filipe/work/code/segm_data');
%%
%load Graz benchmark
%load('graz_data120.mat');
load('graz_data2.mat');
%load('graz_data92.mat');
%%

[clusters50r, ~] = slicmex(image,100,20);
[clusters100r, ~] = slicmex(image,200,20);
[clusters200r, ~] = slicmex(image,400,20);
%[clusters500r, ~] = slicmex(image,1000,2);

segmentations = cat(3,clusters50r,clusters100r,clusters200r)+1;
multi_segs_show
%% load image


% false color composition of the image
false_color_image = image;
% mask
mask = ones(size(false_color_image,1),size(false_color_image,2));
% class probabilities
probabilities  = exp(-data);
% unsupervised segmentations of the HE data
%segmentations = cat(3,clusters50r,clusters100r,clusters200r,clusters500r)+1;

% image of gradients
[x1,x2] = gradient(double(image)/255);
weight_image = exp(-sum((x1.^2 + x2.^2),3)*5);

ground_truth = labels;

figure(120);
colo = colormap('jet');
overimposeimage(image,ground_truth,seeds,0.5,colo);


[c,d] = max(probabilities,[],3);
figure(121);
colo = colormap('jet');
overimposeimage(image,d,seeds,0.5,colo);
%% local clusters
% parameter definition
VIS_TAG = 1;
iterations = 20;
mu = 5;
%% GTV only
Z0 = segsalsa(probabilities,mu,iterations,...
    'GTV',1,'tau_gtv',[5 5 5],'clusters',segmentations,'VIS',VIS_TAG );
[c0,d0] = max(Z0,[],3);

figure(123)
%imshow(d0,[]);colormap('jet')
disp(['GTV only. Acc = ' num2str(mean(d0(:) == ground_truth(:)))]);
colo = colormap('jet');
overimposeimage(image,d0,seeds,0.5,colo);
%% VTV only
Z0 = segsalsa(probabilities,mu,iterations,...
    'VTV',1,'tau_vtv',5,'weight_image_vtv',weight_image,'VIS',VIS_TAG );
[c0,d0] = max(Z0,[],3);

figure(124)
%imshow(d0,[]);colormap('jet')
disp(['VTV only. Acc = ' num2str(mean(d0(:) == ground_truth(:)))]);
colo = colormap('jet');
overimposeimage(image,d0,seeds,0.5,colo);
%% GTV and VTV
Z0 = segsalsa(probabilities,mu,iterations,...
    'VTV',1,'tau_vtv',5,'weight_image_vtv',weight_image,...
    'GTV',1,'tau_gtv',[5 5 5 ],'clusters',segmentations,'VIS',1);
[c0,d0] = max(Z0,[],3);

figure(125)
%imshow(d0,[]);colormap('jet')
disp(['GTV and VTV. Acc = ' num2str(mean(d0(mask(:)) == ground_truth(mask(:))))]);
colo = colormap('jet');
overimposeimage(image,d0,seeds,0.5,colo);
%% STR only
Z0 = segsalsa(probabilities,mu,iterations,...
    'STR',1,'tau_str',5,'weight_image_str',weight_image,'gamma_str',1,'window_str',1,'VIS',VIS_TAG );
[c0,d0] = max(Z0,[],3);

figure(126)
%imshow(d0,[]);colormap('jet')
disp(['STR only. Acc = ' num2str(mean(d0(:) == ground_truth(:)))]);
colo = colormap('jet');
overimposeimage(image,d0,seeds,0.5,colo);

%% STR and VTV
Z0 = segsalsa(probabilities,mu,iterations,...
    'STR',1,'tau_str',1,'weight_image_gtv',weight_image,'gamma_str',1,'window_str',1,...
    'VTV',1,'tau_vtv',5','VIS',VIS_TAG );
[c0,d0] = max(Z0,[],3);

figure(127)
%imshow(d0,[]);colormap('jet')
disp(['STR and VTV Acc = ' num2str(mean(d0(mask(:)) == ground_truth(mask(:))))]);
colo = colormap('jet');
overimposeimage(image,d0,seeds,0.5,colo);
%% STR and GTV
Z0 = segsalsa(probabilities,mu,iterations,...
    'STR',1,'tau_str',1,'weight_image_str',weight_image,'gamma_str',1,'window_str',1,...
    'GTV',1,'tau_gtv',[5 5 5 ],'clusters',segmentations,'VIS',VIS_TAG );
[c0,d0] = max(Z0,[],3);

figure(128)
%imshow(d0,[]);colormap('jet')
disp(['STR and GTV Acc = ' num2str(mean(d0(:) == ground_truth(:)))]);
colo = colormap('jet');
overimposeimage(image,d0,seeds,0.5,colo);