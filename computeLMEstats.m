function pValue=computeLMEstats(stats,mode,groupanimals,test,plotit,logit,removeOut,quant)
% P. Bauer 2020

genotypes=[stats.genotype];
animals={stats.animal};
if(groupanimals)
    [C,ia,ic]=unique(animals(:));
    uniqueid=ic;
else
    for i=1:length(animals)
        uniqueid(i)=i;
    end
end

Y=[];
ids=[];
group=[];
lenghts=zeros(1,length(genotypes));
quantiles=zeros(1,length(genotypes));
maxsize=0;
for i=1:length(genotypes)
    if strcmp(mode,'mean_activiy_place')
        Y1=stats(i).meandFFstats(stats(i).placeCells);
    elseif strcmp(mode,'mean_dff_placecells')
        Y1=stats(i).PlaceScore(stats(i).placeCells);         
    elseif strcmp(mode,'mean_activiy_noplace')
        Y1=stats(i).meandFFstats(stats(i).noPlaceCells);
    elseif strcmp(mode,'mean_activiy_all')
        Y1=stats(i).meandFFstats;
    elseif strcmp(mode,'mean_dff_offield')
        Y1=stats(i).mean_dff_offield;
    elseif strcmp(mode,'std_dff_offield')
        Y1=stats(i).std_dff_offield;        
    elseif strcmp(mode,'mean_dff_field')
        Y1=stats(i).mean_dff_field;
    elseif strcmp(mode,'mean_fieldwidths')
        Y1=stats(i).fieldwidths;
    elseif strcmp(mode,'placescore')
        Y1=stats(i).PlaceScore;
    elseif strcmp(mode,'placescorepct')
        Y1=stats(i).PlaceScorePct;
    end
    if logit
        Y1=log(Y1);
    end
    dists{i}=Y1;
    lenghts(i)=length(Y1);
    quantiles(i)=quantile(Y1,quant);
    Y=[Y Y1];
    muh=repmat(uniqueid(i),length(Y1),1);
    ids=[ids muh'];
    muh=repmat(genotypes(i)+1,length(Y1),1);
    group=[group muh'];
    if length(Y1)>maxsize
        maxsize=length(Y1);
    end
end

if strcmp(test,'bootstrap')
    %additional loop to set group matrices
    group1=nan(length(find(genotypes==0)),maxsize);
    group2=nan(length(find(genotypes==1)),maxsize);
    group1cnt=1;
    group2cnt=1;
    for i=1:length(genotypes)
        if strcmp(mode,'mean_activiy_place')
            Y1=stats(i).meandFFstats(stats(i).placeCells);
        elseif strcmp(mode,'mean_dff_placecells')
            Y1=stats(i).PlaceScore(stats(i).placeCells);  
        elseif strcmp(mode,'mean_activiy_noplace')
            Y1=stats(i).meandFFstats(stats(i).noPlaceCells);
        elseif strcmp(mode,'mean_activiy_all')
            Y1=stats(i).meandFFstats;
        elseif strcmp(mode,'mean_dff_offield')
            Y1=stats(i).mean_dff_offield;
        elseif strcmp(mode,'std_dff_offield')
            Y1=stats(i).std_dff_offield;    
        elseif strcmp(mode,'mean_dff_field')
            Y1=stats(i).mean_dff_field;
        elseif strcmp(mode,'mean_fieldwidths')
            Y1=stats(i).fieldwidths;
        elseif strcmp(mode,'placescore')
            Y1=stats(i).PlaceScore;
        elseif strcmp(mode,'placescorepct')
            Y1=stats(i).PlaceScorePct;    
        end
        if logit
            Y1=log(Y1);
        end
        if genotypes(i)==0
            group1(group1cnt,1:length(Y1))=Y1;
            group1cnt=group1cnt+1;
        else
            group2(group2cnt,1:length(Y1))=Y1;
            group2cnt=group2cnt+1;
        end
    end
end

if removeOut
    inds=find(Y<-4);
    Y(inds)=[];
    ids(inds)=[];
    group(inds)=[];
end

tbl = table();
tbl.Y = Y';
tbl.animal = categorical(ids)';
tbl.group = categorical(group)';

if plotit
    %figure
    subplot(2,2,1)
    hist(Y(group==1),40)
    hold on
    hist(Y(group==2),40,'r')
    subplot(2,2,2)
    ksdensity(Y(group==1))
    hold on
    ksdensity(Y(group==2))    
    
    ax1=subplot(2,2,3);
    hist(lenghts(genotypes>0))
    ax2=subplot(2,2,4);
    hist(lenghts(genotypes<1))
    linkaxes([ax1,ax2],'x')
end

if strcmp(test,'nested_anova')
    [P,T,STATS,TERMS] = anovan(tbl.Y,{tbl.group,tbl.animal},'random',2,'varnames',{'group','animal'},'nested',[0,0;1 0]);
    pValue=P(1);
elseif strcmp(test,'fitlme')
    mdl = fitlme(tbl,'Y~group+(1|animal:group)');
    stats=anova(mdl,'dfmethod','satterthwaite');
    pValue=stats{2,5};
elseif strcmp(test,'ks')
    [h,pValue]=kstest2(Y(group==1),Y(group==2));
elseif strcmp(test,'wilcoxon')    
    pValue=ranksum(Y(group==1),Y(group==2));
elseif strcmp(test,'export')  
    pValue=0;
    csvwrite('stats.csv',[Y;ids;group]');
elseif strcmp(test,'ttest')  
    [h,pValue]=ttest2(Y(group==1),Y(group==2));
elseif strcmp(test,'quantiles')   
    [h,pValue]=ttest2(quantiles(genotypes>0),quantiles(genotypes<1));
elseif strcmp(test,'bootstrap')  
    num_reps = 10000;
    param = 'median';
    num_trials=125;
    [p_boot_amp, bootstats_amp,bootstats_median_amp, bootstats_sem_amp] = get_bootstrap_results_equalsamples(group1,group2,num_reps,num_trials,param);
    pValue=1-p_boot_amp;
    if plotit
        figure 
        cmp = [0.713725490196078,0.427450980392157,1;1,0.600000000000000,0];
        histogram(bootstats_amp(1,:),'Normalization','probability','FaceColor',cmp(1,:),'EdgeColor','none');
        hold on
        histogram(bootstats_amp(2,:),'Normalization','probability','FaceColor',cmp(2,:),'EdgeColor','none');
        title([mode ' p=' num2str(pValue)])
    end
end

end