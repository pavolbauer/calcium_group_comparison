function groupstats=runGroupStats(stats,mode,folder,plotit,test,events,plottraces)
% P. Bauer 2020

savepath = [folder '/' mode];
if ~exist(savepath,'dir') > 0
    try
        mkdir(savepath);
    catch
        error('could not create savedir.')
    end
end

groupstats.empty=0;

if events
    eventfix=2.5;
else
    eventfix=1;
end

%speed stats - positive correlation
WT=find([stats.genotype]==0);
TG=find([stats.genotype]==1);
genotypes=[stats.genotype];
groupnames=cell(1,length(genotypes));
for i=1:length(genotypes)
    if genotypes(i)
        groupnames{i}='TG';
    else
        groupnames{i}='WT';
    end
end

%% mean activity place
xq=0.01:0.1:5.01;
f_interps_tg=[];
f_interps_wt=[];
for i=1:length(genotypes)
    if genotypes(i)
        %h1=plot(stats(i).xi_mean_placecells,stats(i).f_mean_placecells,'r');
        f_interp=interp1(stats(i).xi_mean_placecells,stats(i).f_mean_placecells,xq);
        f_interps_tg=[f_interps_tg; f_interp];
        hold on
    else
        %h2=plot(stats(i).xi_mean_placecells,stats(i).f_mean_placecells,'b');
        f_interp=interp1(stats(i).xi_mean_placecells,stats(i).f_mean_placecells,xq);
        f_interps_wt=[f_interps_wt; f_interp];
        hold on
    end    
end

mean_act_tg = nanmean(f_interps_tg,1);
mean_act_wt = nanmean(f_interps_wt,1);

std_act_tg = nanstd(f_interps_tg,[],1);%./sqrt(size(f_interps_tg,1));
std_act_wt = nanstd(f_interps_wt,[],1);%./sqrt(size(f_interps_wt,1));

X=[xq';flipud(xq')];
U=[mean_act_tg+std_act_tg];
L=[mean_act_tg-std_act_tg];
Ytg=[U';flipud(L')];
U=[mean_act_wt+std_act_wt];
L=[mean_act_wt-std_act_wt];
Ywt=[U';flipud(L')];

fill(X.*eventfix,Ytg,'r','FaceAlpha',.2,'linestyle','none');   
hold on
fill(X.*eventfix,Ywt,'b','FaceAlpha',.2,'linestyle','none');    

if plottraces
    for i=1:length(genotypes)
        if genotypes(i)
            h1=plot(stats(i).xi_mean_placecells.*eventfix,stats(i).f_mean_placecells,'r');
        else
            h2=plot(stats(i).xi_mean_placecells.*eventfix,stats(i).f_mean_placecells,'b');
        end    
    end
end

h1=plot(xq.*eventfix,mean_act_tg,'r','LineWidth',4);
h2=plot(xq.*eventfix,mean_act_wt,'b','LineWidth',4);

if strcmp(mode,'S')
    xlim([0 2.5].*eventfix)
else
    xlim([0 1])
end
ylim([0 1.1])
pValue_mean_activiy_place=computeLMEstats(stats,'mean_activiy_place',0,test,0,0,0,0.3);
title(num2str(pValue_mean_activiy_place))
legend([h1(1),h2(1)],'TG','WT')
xlabel('mean weighted event rate (Hz)')
ylabel('% of place cells')
hline(0.3,'--k')
hline(0.8,'--k')
%title(['pValue=' num2str(pValue)]);
fn=sprintf('%s/mean_place2.png',savepath);
set(gcf,'Position',[0,0,600,600])
print(fn,'-dpng','-r400')
fn=sprintf('%s/mean_place2.pdf',savepath);
print(fn,'-dpdf')
fn=sprintf('%s/mean_place.fig',savepath);
savefig(fn);
close   

%% mean activity no place
xq=0.01:0.1:5.01;
f_interps_tg=[];
f_interps_wt=[];
for i=1:length(genotypes)
    if genotypes(i)
        %h1=plot(stats(i).xi_mean_placecells,stats(i).f_mean_placecells,'r');
        f_interp=interp1(stats(i).xi_mean_notplacecells,stats(i).f_mean_notplacecells,xq);
        f_interps_tg=[f_interps_tg; f_interp];
        hold on
    else
        %h2=plot(stats(i).xi_mean_placecells,stats(i).f_mean_placecells,'b');
        f_interp=interp1(stats(i).xi_mean_notplacecells,stats(i).f_mean_notplacecells,xq);
        f_interps_wt=[f_interps_wt; f_interp];
        hold on
    end    
end

mean_act_tg = nanmean(f_interps_tg,1);
mean_act_wt = nanmean(f_interps_wt,1);

std_act_tg = nanstd(f_interps_tg,[],1);%./sqrt(size(f_interps_tg,1));
std_act_wt = nanstd(f_interps_wt,[],1);%./sqrt(size(f_interps_wt,1));

X=[xq';flipud(xq')];
U=[mean_act_tg+std_act_tg];
L=[mean_act_tg-std_act_tg];
Ytg=[U';flipud(L')];
U=[mean_act_wt+std_act_wt];
L=[mean_act_wt-std_act_wt];
Ywt=[U';flipud(L')];

fill(X.*eventfix,Ytg,'r','FaceAlpha',.2,'linestyle','none');   
hold on
fill(X.*eventfix,Ywt,'b','FaceAlpha',.2,'linestyle','none');    

if plottraces
    for i=1:length(genotypes)
        if genotypes(i)
            h1=plot(stats(i).xi_mean_notplacecells.*eventfix,stats(i).f_mean_notplacecells,'r');
        else
            h2=plot(stats(i).xi_mean_notplacecells.*eventfix,stats(i).f_mean_notplacecells,'b');
        end    
    end
end

h1=plot(xq.*eventfix,mean_act_tg,'r','LineWidth',4);
h2=plot(xq.*eventfix,mean_act_wt,'b','LineWidth',4);

if strcmp(mode,'S')
    xlim([0 2.5].*eventfix)
else
    xlim([0 1])
end
ylim([0 1.1])
pValue_mean_activity_noplace=computeLMEstats(stats,'mean_activiy_noplace',0,test,0,0,0,0.3);

title(num2str(pValue_mean_activity_noplace))
legend([h1(1),h2(1)],'TG','WT')
xlabel('mean weighted event rate (Hz)')
ylabel('% of non-place cells')
hline(0.3,'--k')
hline(0.8,'--k')
%title(['pValue=' num2str(pValue)]);
fn=sprintf('%s/mean_noplace2.png',savepath);
set(gcf,'Position',[0,0,600,600])
print(fn,'-dpng','-r400')
fn=sprintf('%s/mean_noplace2.pdf',savepath);
print(fn,'-dpdf')
fn=sprintf('%s/mean_place.fig',savepath);
savefig(fn);
close

%% average stability
xq=0.01:0.01:1;
plotECDFs(stats,xq,1,'mean_stability',plottraces)

xlabel('Mean stability in field')
fn=sprintf('%s/mean_stability_in_field.png',savepath);
set(gcf,'Position',[0,0,600,600])
print(fn,'-dpng','-r400')
fn=sprintf('%s/mean_stability_in_field.pdf',savepath);
print(fn,'-dpdf')
fn=sprintf('%s/mean_stability_in_field.fig',savepath);
savefig(fn);
close  

%% mean place cell firing
xq=1:500;
plotECDFs(stats,xq,1,'mean_dff_placecells',plottraces)

xlabel('Mean activity place cells')
fn=sprintf('%s/mean_dff_placecells.png',savepath);
set(gcf,'Position',[0,0,600,600])
print(fn,'-dpng','-r400')
fn=sprintf('%s/mean_dff_placecells.pdf',savepath);
print(fn,'-dpdf')
fn=sprintf('%s/mean_dff_placecells.fig',savepath);
savefig(fn);
close  

%% average fieldwidth
xq=1:0.1:20;
plotECDFs(stats,xq,8,'mean_fieldwidths',plottraces)

pValue_mean_fieldwidths=computeLMEstats(stats,'mean_fieldwidths',0,test,0,0,0,0.3);
title(num2str(pValue_mean_fieldwidths))
xlabel('Mean fieldwidths (cm)')
fn=sprintf('%s/mean_fieldwidths.png',savepath);
set(gcf,'Position',[0,0,600,600])
print(fn,'-dpng','-r400')
fn=sprintf('%s/mean_fieldwidths.pdf',savepath);
print(fn,'-dpdf')
fn=sprintf('%s/mean_fieldwidths.fig',savepath);
savefig(fn);
close  

%% average positive speed correlation
xq=0.01:0.01:0.5;
plotECDFs(stats,xq,1,'mean_posspeed',plottraces)

xlabel('Positive speed correlation')
fn=sprintf('%s/mean_positive_speedcorr.png',savepath);
set(gcf,'Position',[0,0,600,600])
print(fn,'-dpng','-r400')
fn=sprintf('%s/mean_positive_speedcorr.pdf',savepath);
print(fn,'-dpdf')
fn=sprintf('%s/mean_positive_speedcorr.fig',savepath);
savefig(fn);
close  

%% average negative speed correlation
if 0
xq=0.01:0.01:0.5;
plotECDFs(stats,xq,1,'mean_negspeed',plottraces)

xlabel('Negative speed correlation')
fn=sprintf('%s/mean_negative_speedcorr.png',savepath);
set(gcf,'Position',[0,0,600,600])
print(fn,'-dpng','-r400')
fn=sprintf('%s/mean_negative_speedcorr.fig',savepath);
savefig(fn);
close  
end
%% average in-field dFF
xq=1:500;
plotECDFs(stats,xq,eventfix/100,'mean_dff_field',plottraces)

pvalue_mean_dff_field=computeLMEstats(stats,'mean_dff_field',0,test,0,0,0,0.3);
title(num2str(pvalue_mean_dff_field));
xlabel('mean weighted event rate (Hz) in field')
fn=sprintf('%s/mean_in_field.png',savepath);
set(gcf,'Position',[0,0,600,600])
print(fn,'-dpng','-r400')
fn=sprintf('%s/mean_in_field.pdf',savepath);
print(fn,'-dpdf')
fn=sprintf('%s/mean_in_field.fig',savepath);
savefig(fn);
close   

%% std in-field dFF
xq=1:1000;
plotECDFs(stats,xq,eventfix/100,'std_dff_field',plottraces)

xlabel('std weighted event rate (Hz) in field')
fn=sprintf('%s/std_in_field.png',savepath);
set(gcf,'Position',[0,0,600,600])
print(fn,'-dpng','-r400')
fn=sprintf('%s/std_in_field.pdf',savepath);
print(fn,'-dpdf')
fn=sprintf('%s/std_in_field.fig',savepath);
savefig(fn);
close   

%% CV in-field dFF
if 0
xq=1:0.1:100;
plotECDFs(stats,xq,eventfix,'cv_dff_field',plottraces)

xlabel('CV activity in field')
fn=sprintf('%s/cv_in_field.png',savepath);
set(gcf,'Position',[0,0,600,600])
print(fn,'-dpng','-r400')
fn=sprintf('%s/cv_in_field.fig',savepath);
savefig(fn);
close   
end
%% average off-field dFF
xq=1:150;
plotECDFs(stats,xq,eventfix/100,'mean_dff_offield',plottraces)

pValue_mean_dff_offield=computeLMEstats(stats,'mean_dff_offield',0,test,0,0,0,0.3);
title(num2str(pValue_mean_dff_offield))
xlabel('mean weighted event rate (Hz) off field')
fn=sprintf('%s/mean_off_field.png',savepath);
set(gcf,'Position',[0,0,600,600])
print(fn,'-dpng','-r400')
fn=sprintf('%s/mean_off_field.pdf',savepath);
print(fn,'-dpdf')
fn=sprintf('%s/mean_off_field.fig',savepath);
savefig(fn);
close   

%% std off-field dFF
xq=1:1000;
plotECDFs(stats,xq,eventfix/100,'std_dff_offield',plottraces)

pValue_std_dff_offield=computeLMEstats(stats,'std_dff_offield',0,test,0,0,0,0.3);
title(num2str(pValue_std_dff_offield))
xlabel('std weighted event rate (Hz) off field')
fn=sprintf('%s/std_off_field.png',savepath);
set(gcf,'Position',[0,0,600,600])
print(fn,'-dpng','-r400')
fn=sprintf('%s/std_off_field.pdf',savepath);
print(fn,'-dpdf')
fn=sprintf('%s/std_off_field.fig',savepath);
savefig(fn);
close   

%% place score
xq=0.1:0.1:5;
plotECDFs(stats,xq,eventfix/100,'placescore',plottraces)

pValue_std_dff_offield=computeLMEstats(stats,'placescore',0,test,0,0,0,0.3);
title(num2str(pValue_std_dff_offield))
xlabel('std weighted event rate (Hz) off field')
fn=sprintf('%s/placescore.png',savepath);
set(gcf,'Position',[0,0,600,600])
print(fn,'-dpng','-r400')
fn=sprintf('%s/placescore.pdf',savepath);
print(fn,'-dpdf')
fn=sprintf('%s/placescore.fig',savepath);
savefig(fn);
close

%% place score pct
xq=0.01:0.01:1;
plotECDFs(stats,xq,eventfix/100,'placescorepct',plottraces)

pValue_std_dff_offield=computeLMEstats(stats,'placescorepct',0,test,0,0,0,0.3);
title(num2str(pValue_std_dff_offield))
xlabel('std weighted event rate (Hz) off field')
fn=sprintf('%s/placescorepct.png',savepath);
set(gcf,'Position',[0,0,600,600])
print(fn,'-dpng','-r400')
fn=sprintf('%s/placescorepct.pdf',savepath);
print(fn,'-dpdf')
fn=sprintf('%s/placescorepct.fig',savepath);
savefig(fn);
close

%% CV off-field dFF
if 0
xq=1:0.1:100;
plotECDFs(stats,xq,eventfix,'cv_dff_offield',plottraces)

xlabel('CV activity off field')
fn=sprintf('%s/cv_off_field.png',savepath);
set(gcf,'Position',[0,0,600,600])
print(fn,'-dpng','-r400')
fn=sprintf('%s/cv_off_field.fig',savepath);
savefig(fn);
close   
end
%% speed stats
if 0
avgPosCorr=[stats.avgPosCorr];
avgPosSlope=[stats.avgPosSlope];
stdPosCorr=[stats.stdPosCorr];

h1=scatter(avgPosCorr(WT),avgPosSlope(WT),stdPosCorr(WT).*1000,'ob','filled','MarkerFaceAlpha',.5);
hold on
h2=scatter(avgPosCorr(TG),avgPosSlope(TG),stdPosCorr(TG).*1000,'or','filled','MarkerFaceAlpha',.5);
legend([h1(1),h2(1)],'WT','TG')
xlabel('positive correlation coefficient')
ylabel('positive fitting slope')
fn=sprintf('%s/pos_speed_std.png',savepath);
set(gcf,'Position',[0,0,1200,600])
print(fn,'-dpng','-r400')
fn=sprintf('%s/pos_speed_std.fig',savepath);
savefig(fn);
close    

scatterhist(avgPosCorr,avgPosSlope,'Group',groupnames','Kernel','on');
xlabel('positive correlation coefficient')
ylabel('positive fitting slope')
fn=sprintf('%s/pos_speed.png',savepath);
set(gcf,'Position',[0,0,1200,600])
print(fn,'-dpng','-r400')
fn=sprintf('%s/pos_speed.fig',savepath);
savefig(fn);
close    

%speed stats - negative correlation
avgNegCorr=[stats.avgNegCorr];
avgNegSlope=[stats.avgNegSlope];
stdNegCorr=[stats.stdNegCorr];

h1=scatter(avgNegCorr(WT),avgNegSlope(WT),stdNegCorr(WT).*1000,'ob','filled','MarkerFaceAlpha',.5);
hold on
h2=scatter(avgNegCorr(TG),avgNegSlope(TG),stdNegCorr(TG).*1000,'or','filled','MarkerFaceAlpha',.5);
legend([h1(1),h2(1)],'WT','TG')
xlabel('positive correlation coefficient')
ylabel('positive fitting slope')
fn=sprintf('%s/neg_speed_std.png',savepath);
set(gcf,'Position',[0,0,1200,600])
print(fn,'-dpng','-r400')
fn=sprintf('%s/neg_speed_std.fig',savepath);
savefig(fn);
close   

scatterhist(avgNegCorr,avgNegSlope,'Group',groupnames','Kernel','on');
xlabel('positive correlation coefficient')
ylabel('positive fitting slope')
fn=sprintf('%s/neg_speed.png',savepath);
set(gcf,'Position',[0,0,1200,600])
print(fn,'-dpng','-r400')
fn=sprintf('%s/neg_speed.fig',savepath);
savefig(fn);
close   


%lick score
lickscore=ones(1,length(genotypes));
for i=1:length(genotypes)
    lickscore(i)=mean(stats(i).dFFlickScore);
end
ksdensity(lickscore(WT))
hold on
ksdensity(lickscore(TG))
legend('WT','TG')
xlabel('lick vs no-lick activity')
fn=sprintf('%s/lickscore.png',savepath);
set(gcf,'Position',[0,0,1200,600])
print(fn,'-dpng','-r400')
fn=sprintf('%s/lickscore.fig',savepath);
savefig(fn);
close   


%rest vs run
runscore=ones(1,length(genotypes));
for i=1:length(genotypes)
    runscore(i)=nanmean(stats(i).dFFRest./(stats(i).dFFRest+stats(i).dFFRun));
end
ksdensity(runscore(WT))
hold on
ksdensity(runscore(TG))
legend('WT','TG')
xlabel('rest vs run activity')
fn=sprintf('%s/runscore.png',savepath);
set(gcf,'Position',[0,0,1200,600])
print(fn,'-dpng','-r400')
fn=sprintf('%s/runscore.fig',savepath);
savefig(fn);
close   
end