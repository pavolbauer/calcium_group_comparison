function fixed_trace=checkCorrelatedComponents(trace,thresh)
% P. Bauer 2020

dFF_out_raw=dFF_out;

%remove correlated components
corrMat=zeros(size(dFF_out_raw,1),size(dFF_out_raw,1));
for i=1:size(dFF_out_raw,1)
    for j=1:size(dFF_out_raw,1)
        r=corrcoef(dFF_out_raw(i,:),dFF_out_raw(j,:));
        corrMat(i,j)=r(1,2);
    end
    i
    ends

%remove diagonal
for i=1:size(dFF_out_raw,1)
    corrMat(i,i)=0;
end

%cluster data
D = pdist(dFF_out_raw,'spearman');
Z=linkage(D);
c=cluster(Z,'maxclust',20);
dendrogram(Z)

T=clusterdata(1-corrMat,1);
I=find(T==541);
imagesc(dFF_out_raw(I,:));

orderedMat=zeros(size(dFF_out_raw,1),size(dFF_out_raw,1));
uniqueMat=[];
cnt=1;
for i=1:max(T)
    [r,c]=find(T==i);
    if ~isempty(r)
        uniquer=unique(r);
        for j=1:length(uniquer)
           if j==1
               uniqueMat(end+1,:)=corrMat(uniquer(j),:);
           end
           orderedMat(cnt,:)=corrMat(uniquer(j),:);
           cnt=cnt+1;
        end
    end
end

%# and convert to a vector (as pdist)
%dissimilarity = 1 - corrMat(find(corrMat))';
dissimilarity = 1 - corrMat;

%# decide on a cutoff
%# remember that 0.4 corresponds to corr of 0.6!
cutoff = 0.5; 

%# perform complete linkage clustering
Z = linkage(dissimilarity,'complete');

%# group the data into clusters
%# (cutoff is at a correlation of 0.5)
groups = cluster(Z,'cutoff',cutoff,'criterion','distance');

dendrogram(Z,0,'colorthreshold',cutoff)

end

