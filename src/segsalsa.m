% fcondessa@gmail.com



function Zim = segsalsa(P,mu,iters,varargin)

[nr,nc,no_classes] = size(P);
n = nr*nc;



%set default flags
FLAG_GTV =0;
FLAG_VTV = 0;
FLAG_STR = 0;
FLAG_VIS = 1;
%set default arguments
no_segmentations = 0;
tau_vtv = 10;
tau_gtv = 1;
tau_str = 1;
gamma_str = 1;
Lk = 1;
L = (2*Lk+1)^2;
weight_image_vtv = ones(nr,nc);
weight_image_vtv = weight_image_vtv(:)';
weight_image_str = ones(nr,nc);
weight_image_str = weight_image_str(:)';
clusters = [];

nVarargs = length(varargin);
assert(rem(nVarargs,2) == 0,'segsalsa: arguments must be in pairs after the number of iterations')
    
for i = 1:2:nVarargs,
    arg_type = varargin{i};
    arg_val = varargin{i+1};
    switch arg_type,
        case 'GTV',
            FLAG_GTV = arg_val;
        case 'VTV',
            FLAG_VTV = arg_val;            
        case 'STR',
            FLAG_STR = arg_val;                
        case 'VIS',
            FLAG_VIS = arg_val;
        case 'tau_vtv',
            tau_vtv = arg_val;            
        case 'tau_gtv',
            tau_gtv = arg_val;
        case 'tau_str',
            tau_str = arg_val;    
        case 'gamma_str',
            gamma_str = arg_val;
        case 'window_str',
            Lk = arg_val;
        case 'weight_image_vtv',
            weight_image_vtv = arg_val;
            weight_image_vtv = weight_image_vtv(:)';
        case 'weight_image_str',
            weight_image_str = arg_val;
            weight_image_str = weight_image_str(:)';
        case 'clusters',
            clusters = arg_val;           
    end
end


no_segmentations = size(clusters,3);
dims = [no_classes,n ,no_segmentations];

if FLAG_GTV;
clusters = reshape(clusters,n,no_segmentations);
    if length(tau_gtv)<no_segmentations,
    tau_gtv = tau_gtv * ones(no_segmentations,1);
    end
end

weight_image_vtv = reshape(weight_image_vtv,nr,nc);
weight_image_vtv(1,:) = 0;
weight_image_vtv(:,1) = 0;
weight_image_str = reshape(weight_image_str,nr,nc);
weight_image_str(1,:) = 0;
weight_image_str(:,1) = 0;

aux_t1 = ones(no_classes,n,no_segmentations);
for i = 1:no_segmentations,
    aux_t1(:,:,i) = 2* aux_t1(:,:,i)*tau_gtv(i) / mu;
end

P = reshape(P, n,no_classes)';
%% defines
% define operators
ConvC = @(X,FK,nb)  reshape(real(ifft2(fft2(reshape(X', nr,nc,nb)).*repmat(FK,[1,1,nb]))), nr*nc,nb)';
conv2im  = @(X,nb)  reshape(X',nr,nc,nb);
conv2mat = @(X,nb)  reshape(X,nr*nc,nb)';
% define horiontal difference operator kernel
dh = zeros(nr,nc);
dh(1,1) = 1;
dh(1,nc) = -1;
dv = zeros(nr,nc);
dv(1,1) = 1;
dv(nr,1) = -1;
% define FFTs for filtering
FDH = fft2(dh);
FDHC = conj(FDH);
FDV = fft2(dv);
FDVC = conj(FDV);

FDH_STR = zeros(size(dh));
FDV_STR = zeros(size(dv));
if FLAG_STR,
    [X,Y] = meshgrid(-Lk:Lk);
    K = exp(-(X.^2+Y.^2)/gamma_str);
    K = K/max(K(:));
    L = length(K(:));
    for i = 1:L,
        FDH_STR = FDH_STR + tau_str * K(i) * fft2(circshift(dh,[X(i) Y(i)]));
        FDV_STR = FDV_STR + tau_str * K(i) * fft2(circshift(dv,[X(i) Y(i)]));
        FDH1_STR{i} = tau_str * K(i)*fft2(circshift(dh,[X(i) Y(i)]));
        FDV1_STR{i} = tau_str * K(i)*fft2(circshift(dv,[X(i) Y(i)]));
    end
end
FDHC_STR = conj(FDH_STR);
FDHV_STR = conj(FDV_STR);
%%  define normalization filters in the frequency domain

num_ids_ops = 3;
if FLAG_GTV,
    num_ids_ops = num_ids_ops + no_segmentations;
end

norm_filt = (FLAG_VTV*(abs(FDH).^2+ abs(FDV).^2 )+ num_ids_ops + FLAG_STR*(abs(FDV_STR).^2 + abs(FDH_STR).^2)/L);
I_DH = FDHC ./norm_filt;
I_DV = FDVC ./norm_filt;
II   = 1  ./norm_filt;
FDH = ifft2(I_DH);
FDV = ifft2(I_DV);
FI  = ifft2(II);


%% some constants
P2 =  sum(P.^2);
Iu = eye(no_classes)- 1/no_classes*ones(no_classes);
B = ones(no_classes,n)/no_classes;
% define and initialize variables
Z = P;
V_data = Z;
D_data = 0*Z;
V_non_neg_cstr = Z;
D_non_neg_cstr = 0*Z;
V_sum_to_one_cstr = Z;
D_sum_to_one_cstr = 0*Z;

if FLAG_VTV,
    V_vtv_prior_hz_comp = Z;
    D_vtv_prior_hz_comp = 0*Z;
    V_vtv_prior_vt_comp = Z;
    D_vtv_prior_vt_comp = 0*Z;
end

if FLAG_GTV,
V_gtv_prior = repmat(Z,[1,1,no_segmentations]);
D_gtv_prior = 0*V_gtv_prior;
end

if FLAG_STR,
    V_str_prior_vt_comp = 1*repmat(Z,[L,1,1]);
    V_str_prior_hz_comp = 1*repmat(Z,[L,1,1]);
    D_str_prior_vt_comp = 0*V_str_prior_vt_comp;
    D_str_prior_hz_comp = 0*V_str_prior_hz_comp;
end
for i=1:iters

% solve quadratic problem

% this is divided based on the definition of the linear operators
Z_CORE = V_data+D_data+V_non_neg_cstr+D_non_neg_cstr+V_sum_to_one_cstr+D_sum_to_one_cstr;

if FLAG_GTV,
    Z_CORE  = Z_CORE   + sum(V_gtv_prior+D_gtv_prior,3);
end
Z = ConvC(Z_CORE, II,no_classes);


if FLAG_VTV,
    Z = Z + ConvC(V_vtv_prior_hz_comp+D_vtv_prior_hz_comp, I_DH,no_classes) + ...
        ConvC(V_vtv_prior_vt_comp+D_vtv_prior_vt_comp, I_DV,no_classes);
end

if FLAG_STR,
    for j = 1:L,
%         a = circshift(ConvC(circshift(reshape(V_str_prior_vt_comp((j-1)*no_classes+(1:no_classes),:)+D_str_prior_vt_comp((j-1)*no_classes+(1:no_classes),:),[no_classes,nc*nr]),[X(j) Y(j)]), conj(FDV1_STR{j})./(norm_filt),no_classes),[X(j) -Y(j)])/L;
%         b = circshift(ConvC(circshift(reshape(V_str_prior_hz_comp((j-1)*no_classes+(1:no_classes),:)+D_str_prior_hz_comp((j-1)*no_classes+(1:no_classes),:),[no_classes,nc*nr]),[X(j) Y(j)]), conj(FDH1_STR{j})./(norm_filt),no_classes),[X(j) -Y(j)])/L;
%         subplot(1,2,1);imagesc(reshape(a(1,:),nr,nc));
%         subplot(1,2,2);imagesc(reshape(b(1,:),nr,nc));
%         pause(.1);
        Z = Z + circshift(ConvC(reshape(V_str_prior_vt_comp((j-1)*no_classes+(1:no_classes),:)+D_str_prior_vt_comp((j-1)*no_classes+(1:no_classes),:),[no_classes,nc*nr]), conj(FDV1_STR{j})./(norm_filt),no_classes),-[X(j) Y(j)])/L;
        Z = Z + circshift(ConvC(reshape(V_str_prior_hz_comp((j-1)*no_classes+(1:no_classes),:)+D_str_prior_hz_comp((j-1)*no_classes+(1:no_classes),:),[no_classes,nc*nr]), conj(FDH1_STR{j})./(norm_filt),no_classes),-[X(j) Y(j)])/L;
        %Z = Z + circshift(ConvC(circshift(reshape(V_str_prior_vt_comp((j-1)*no_classes+(1:no_classes),:)+D_str_prior_vt_comp((j-1)*no_classes+(1:no_classes),:),[no_classes,nc*nr]),[X(j) Y(j)]), conj(FDV1_STR{j})./(norm_filt),no_classes),-[X(j) Y(j)])/L;
        %Z = Z + circshift(ConvC(circshift(reshape(V_str_prior_hz_comp((j-1)*no_classes+(1:no_classes),:)+D_str_prior_hz_comp((j-1)*no_classes+(1:no_classes),:),[no_classes,nc*nr]),[X(j) Y(j)]), conj(FDH1_STR{j})./(norm_filt),no_classes),-[X(j) Y(j)])/L;
    end
end


% solve split variables
% data term
NU_data = Z-D_data;
V_data = compute_prox(NU_data,'data',mu,dims,P,P2);
D_data = -(NU_data)+V_data;

% nonnegativity constraint
V_non_neg_cstr = compute_prox(Z-D_non_neg_cstr,'nonnegativity',mu,dims);
D_non_neg_cstr = -(Z-D_non_neg_cstr) + V_non_neg_cstr;

% sum to one constraint
V_sum_to_one_cstr = compute_prox(Z-D_sum_to_one_cstr,'sumtoone',mu,dims,Iu,B);
D_sum_to_one_cstr = -(Z-D_sum_to_one_cstr) + V_non_neg_cstr;

% VTV prior 
if FLAG_VTV,
    NU_vtv_prior_hz_comp = (ConvC(Z,FDH,no_classes) - D_vtv_prior_hz_comp);
    NU_vtv_prior_vt_comp = (ConvC(Z,FDV,no_classes) - D_vtv_prior_vt_comp);
    V_VTV = compute_prox(cat(3,NU_vtv_prior_hz_comp,NU_vtv_prior_vt_comp),...
        'VTV',mu,dims,tau_vtv,weight_image_vtv(:)');
    V_vtv_prior_hz_comp = V_VTV(:,:,1);
    V_vtv_prior_vt_comp = V_VTV(:,:,2);
    D_vtv_prior_hz_comp = -NU_vtv_prior_hz_comp + V_vtv_prior_hz_comp;
    D_vtv_prior_vt_comp = -NU_vtv_prior_vt_comp + V_vtv_prior_vt_comp;
end

% GTV prior
if FLAG_GTV,
    NU_gtv_prior  = repmat(Z,[1,1,no_segmentations])-D_gtv_prior;
    V_gtv_prior = compute_prox(NU_gtv_prior,...
        'GTV',mu,dims,aux_t1,clusters);
    D_gtv_prior = -NU_gtv_prior + V_gtv_prior;
end

% STR prior
if FLAG_STR,
    [outx,outy] = patch_jacobian(Z,nr,nc,no_classes,Lk,gamma_str);
    NU_str_prior_vt_comp = reshape(outx,n,no_classes*L)' - D_str_prior_vt_comp;
    NU_str_prior_hz_comp = reshape(outy,n,no_classes*L)' - D_str_prior_hz_comp;
    V_STR = compute_prox(cat(3,NU_str_prior_hz_comp',NU_str_prior_vt_comp'),...
        'STR',mu,dims,tau_str,weight_image_str(:));
    V_str_prior_hz_comp = V_STR(:,:,1)';
    V_str_prior_vt_comp = V_STR(:,:,2)';
    
%     plot(NU_str_prior_hz_comp(:)-V_str_prior_hz_comp(:));
%     pause;
    D_str_prior_vt_comp = -NU_str_prior_vt_comp + V_str_prior_vt_comp;
    D_str_prior_hz_comp = -NU_str_prior_hz_comp + V_str_prior_hz_comp;
end

if FLAG_VIS
    Zim = conv2im(Z,no_classes);
    [c0,d0] = max(Zim,[],3);
    imagesc(d0);
    title(['iteration: ' num2str(i)]);pause(0.01);
end
end
Zim = conv2im(Z,no_classes);