%fcondessa
%2016.02.11
function out = compute_mean_clusters(NU,clusters)
NUM_CLUSTERS_TRADEOFF = 300;
if max(clusters(:)) >= NUM_CLUSTERS_TRADEOFF,
    vals = zeros(max(clusters),size(NU,1));
    cts = zeros(max(clusters),1);
    out = (zeros(size(NU,1),size(NU,2)));
    for i = 1:size(NU,2),
        vals(clusters(i),:) = vals(clusters(i),:) + NU(:,i)';
        cts(clusters(i)) = cts(clusters(i)) + 1;
    end
    vals = vals./repmat(cts,1,size(NU,1));
    for j = 1:size(NU,1),
        out(j,:) = vals(clusters,j)';
    end
else
    out = (zeros(size(NU,1),size(NU,2)));
    for j = 1:max(clusters),
        locs = clusters==j;
        out(:,locs) =  out(:,locs)  +repmat(mean(NU(:,locs),2),[1 sum(locs)]);
    end
end
