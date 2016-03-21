% fcondessa@gmail.com
% 2016.03.20

function Vout = compute_prox(NUinp,prox_type,mu,dims,varargin)
no_classes = dims(1);
n = dims(2);
no_segmentations = dims(3);
switch prox_type
    case 'data',
        P = varargin{1};
        P2 = varargin{2};
        pnu = sum(P.*NUinp);
        aux = (pnu + sqrt(pnu.^2 +4*P2/mu))/2;
        Vout = NUinp + P./(repmat(aux*mu,no_classes,1));
    case 'nonnegativity',
        Vout = max(0,NUinp);
    case 'sumtoone',
        Vout = (eye(no_classes)- 1/no_classes*ones(no_classes)) *NUinp + ones(no_classes,n)/no_classes;
    case 'VTV',
        [Vout_h,Vout_v] = vector_soft_col_iso(NUinp(:,:,1),NUinp(:,:,2),varargin{1}/mu*varargin{2}); 
        Vout = cat(3,Vout_h,Vout_v);
    case 'STR',
        [Vout_h,Vout_v] = schatten_reg(NUinp(:,:,1),NUinp(:,:,2),varargin{1}/mu*varargin{2}); 
        Vout = cat(3,Vout_h,Vout_v);        
    case 'GTV'
        FNUaux = zeros(size(NUinp));
        for cluster_id = 1:no_segmentations,
            FNUaux(:,:,cluster_id) = mean_clusters(NUinp(:,:,cluster_id),varargin{2}(:,cluster_id));
        end
        Vout = (NUinp + varargin{1}.* FNUaux)./(1+varargin{1});
    otherwise
        error('compute_prox: invalid prox_type');
end

