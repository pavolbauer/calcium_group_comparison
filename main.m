user='FF';
%row=[23,26,29,32,35,38,39,40,41,84,86, ...
%    87,88,89,90,91,92,94,95,96,97,98,99,100];
row=[101, 102, 107];
plotit=0;

% read in data
clear criterion
criterion.minPShuffle=0.80;
%criterion.minRhoNorm=0.90;
criterion.maxRayleighP=0.05;
results_folder='/home/bauerp/Alzheimer_ray';
mode='dFF'; %dFF,S,C
for i=1:length(row)
    stats1b_dFF(i)=readStatistics(user,row(i),mode,results_folder,plotit,criterion);
end

mode='S';
for i=1:length(row)
    stats1b_S(i)=readStatistics(user,row(i),mode,results_folder,plotit,criterion);
end
% read in data
clear criterion
criterion.minPShuffle=0.80;
criterion.minRhoNorm=0.90;
%criterion.maxRayleighP=0.05;
results_folder='/home/bauerp/Alzheimer_rho';
mode='dFF'; %dFF,S,C
for i=1:length(row)
    stats2b_dFF(i)=readStatistics(user,row(i),mode,results_folder,plotit,criterion);
end

mode='S';
for i=1:length(row)
    stats2b_S(i)=readStatistics(user,row(i),mode,results_folder,plotit,criterion);
end

save('results5.mat')

%% get additional placeCell vectors that I forot to save before
clear criterion
criterion.minPShuffle=0.80;
criterion.minRhoNorm=0.90;
%criterion.maxRayleighP=0.05;
results_folder='/home/bauerp/Alzheimer_ray';
mode='S';
placecell_ratio_wt=[];
placecell_ratio_tg=[];
placecell_count_wt=[];
placecell_count_tg=[];
posspeed_ratio_wt=[];
posspeed_ratio_tg=[];
posspeed_count_wt=[];
posspeed_count_tg=[];
negspeed_ratio_wt=[];
negspeed_ratio_tg=[];
negspeed_count_wt=[];
negspeed_count_tg=[];
other_ratio_wt=[];
other_ratio_tg=[];
other_count_wt=[];
other_count_tg=[];

pcShuffle_wt=[];
pcShuffle_tg=[];
rhoNorm_wt=[];
rhoNorm_tg=[];

for i=1:length(row)
    [placeCells,noPlaceCells,posSpeedCells,negSpeedCells,notClear,centroidrhoNorm,PlaceScorePct]=readPlaceCells(user,row(i),mode,results_folder,0,criterion);
    stats1_S(i).placeCells=placeCells;
    stats1_S(i).noPlaceCells=noPlaceCells;
    if stats1_S(i).genotype
        placecell_ratio_tg(end+1)=length(placeCells)/(length(placeCells)+length(noPlaceCells));
        posspeed_ratio_tg(end+1)=length(posSpeedCells)/(length(placeCells)+length(noPlaceCells));
        negspeed_ratio_tg(end+1)=length(negSpeedCells)/(length(placeCells)+length(noPlaceCells));
        other_ratio_tg(end+1)=length(notClear)/(length(placeCells)+length(noPlaceCells));
        placecell_count_tg(end+1)=length(placeCells);
        posspeed_count_tg(end+1)=length(posSpeedCells);
        negspeed_count_tg(end+1)=length(negSpeedCells);
        other_count_tg(end+1)=length(notClear);
        pcShuffle_tg=[pcShuffle_tg PlaceScorePct];
        rhoNorm_tg=[rhoNorm_tg centroidrhoNorm];
    else
        placecell_ratio_wt(end+1)=length(placeCells)/(length(placeCells)+length(noPlaceCells));
        posspeed_ratio_wt(end+1)=length(posSpeedCells)/(length(placeCells)+length(noPlaceCells));
        negspeed_ratio_wt(end+1)=length(negSpeedCells)/(length(placeCells)+length(noPlaceCells));
        other_ratio_wt(end+1)=length(notClear)/(length(placeCells)+length(noPlaceCells));
        placecell_count_wt(end+1)=length(placeCells);
        posspeed_count_wt(end+1)=length(posSpeedCells);
        negspeed_count_wt(end+1)=length(negSpeedCells);
        other_count_wt(end+1)=length(notClear);  
        pcShuffle_wt=[pcShuffle_wt PlaceScorePct];
        rhoNorm_wt=[rhoNorm_wt centroidrhoNorm];
    end
end

%%
figure
subplot(1,2,1)
labels={['Place tuning (n=' num2str(sum(placecell_count_tg)) ', ' num2str(mean(placecell_ratio_tg)*100,'%2.1f') '%)'], ...
    ['Positive speed tuning  (n=' num2str(sum(posspeed_count_tg)) ', ' num2str(mean(posspeed_ratio_tg)*100,'%2.1f') '%)'], ...
    ['Negative speed tuning (n=' num2str(sum(negspeed_count_tg)) ', ' num2str(mean(negspeed_ratio_tg)*100,'%2.1f') '%)'], ...
    ['Other (n=' num2str(sum(other_count_tg)) ', ' num2str(mean(other_ratio_tg)*100,'%2.1f') '%)']};
pie([mean(placecell_ratio_tg),mean(posspeed_ratio_tg),mean(negspeed_ratio_tg),mean(other_ratio_tg)],labels);
title(['tg (Total=' num2str(sum(placecell_count_tg)+sum(posspeed_count_tg)+sum(negspeed_count_tg)+sum(other_count_tg)) ' components)']);
subplot(1,2,2)
labels={['Place tuning (n=' num2str(sum(placecell_count_wt)) ', ' num2str(mean(placecell_ratio_wt)*100,'%2.1f') '%)'], ...
    ['Positive speed tuning  (n=' num2str(sum(posspeed_count_wt)) ', ' num2str(mean(posspeed_ratio_wt)*100,'%2.1f') '%)'], ...
    ['Negative speed tuning (n=' num2str(sum(negspeed_count_wt)) ', ' num2str(mean(negspeed_ratio_wt)*100,'%2.1f') '%)'], ...
    ['Other (n=' num2str(sum(other_count_wt)) ', ' num2str(mean(other_ratio_wt)*100,'%2.1f') '%)']};
pie([mean(placecell_ratio_wt),mean(posspeed_ratio_wt),mean(negspeed_ratio_wt),mean(other_ratio_wt)],labels);
title(['wt (Total=' num2str(sum(placecell_count_wt)+sum(posspeed_count_wt)+sum(negspeed_count_wt)+sum(other_count_wt)) ' components)']);
%% criterion
figure
subplot(1,2,1)
hist(pcShuffle_tg,30);
vline(0.95,'r')
%hline(0.95,'r')
xlabel('Place Score Shuffle Test Percentile')
%ylabel('Normalized Population Vector Height')
%xlim([0 1])
%ylim([0 1])
title('tg')
subplot(1,2,2)
hist(pcShuffle_wt,30)
vline(0.95,'r')
%hline(0.95,'r')
xlabel('Place Score Shuffle Test Percentile')
%ylabel('Normalized Population Vector Height')
title('wt')
%xlim([0 1])
%ylim([0 1])

%% place cells 
pcwt=[];
pctg=[];
for i=1:length(stats3_S)
    if stats3_S(i).genotype
        pctg(end+1)=length(stats3_S(i).placeCells)/(length(stats3_S(i).placeCells) + length(stats3_S(i).noPlaceCells));
    else
        pcwt(end+1)=length(stats3_S(i).placeCells)/(length(stats3_S(i).placeCells) + length(stats3_S(i).noPlaceCells));
    end
end

ksdensity(pcwt)
hold on
ksdensity(pctg)

%% test
stats_test=readStatistics(user,row(2),'S',results_folder,plotit,criterion);

%% group statistics dFF
test='bootstrap';
results_folder='/home/bauerp/Alzheimer_groupstats1';
mode='dFF';
plottraces=1;
speedStats = runGroupStats(stats1_dFF,mode,results_folder,1,test,0,plottraces);

% group statistics events
results_folder='/home/bauerp/Alzheimer_groupstats1';
mode='S';
speedStats = runGroupStats(stats3_S,mode,results_folder,1,test,1,plottraces);

results_folder='/home/bauerp/Alzheimer_groupstats2';
mode='dFF';
speedStats = runGroupStats(stats2_dFF,mode,results_folder,1,test,0,plottraces);

% group statistics events
results_folder='/home/bauerp/Alzheimer_groupstats_new';
mode='S';
%mode='dFF';
speedStats = runGroupStats(stats3_S,mode,results_folder,1,test,1,plottraces);

%% compute stability
clear criterion
criterion.minPShuffle=0.80;
criterion.minRhoNorm=0.90;
%criterion.maxRayleighP=0.05;
results_folder='/home/bauerp/Alzheimer_ray';
mode='S';

for i=1:length(row)
    [mean_dff_field,std_dff_field,mean_dff_offield,std_dff_offield,stabs,binfields,fieldwidths,placeCells,noPlaceCells]=...
        readStability(user,row(i),mode,results_folder,0,criterion);    
    stats1_S(i).mean_dff_field=mean_dff_field;
    stats1_S(i).std_dff_field=std_dff_field;
    stats1_S(i).mean_dff_offield=mean_dff_offield;
    stats1_S(i).std_dff_offield=std_dff_offield;
    stats1_S(i).stabs=stabs;
    stats1_S(i).fieldwidths=fieldwidths;
    stats1_S(i).binfields=binfields;
    stats1_S(i).placeCells=placeCells;
    stats1_S(i).noPlaceCells=noPlaceCells;
    i
end

stats3_S=stats1_S;
stats3_S(22)=[];

%% compute stability using stability criterion
clear criterion
criterion.minPShuffle=0.80;
criterion.minRhoNorm=0.90;
%criterion.maxRayleighP=0.05;
results_folder='/home/bauerp/Alzheimer_ray';
mode='S';

for i=1:length(row)
    [mean_dff_field,std_dff_field,mean_dff_offield,std_dff_offield,stabs,binfields,fieldwidths,placeCells,noPlaceCells]=...
        readStability(user,row(i),mode,results_folder,0,criterion);    
    stats1_S(i).mean_dff_field=mean_dff_field;
    stats1_S(i).std_dff_field=std_dff_field;
    stats1_S(i).mean_dff_offield=mean_dff_offield;
    stats1_S(i).std_dff_offield=std_dff_offield;
    stats1_S(i).stabs=stabs;
    stats1_S(i).fieldwidths=fieldwidths;
    stats1_S(i).binfields=binfields;
    stats1_S(i).placeCells=placeCells;
    stats1_S(i).noPlaceCells=noPlaceCells;
    i
end

stats3_S=stats1_S;
stats3_S(22)=[];