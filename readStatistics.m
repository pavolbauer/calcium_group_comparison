function stats=readStatistics(user,row,mode,folder,plotit,criterion)
% P. Bauer 2020

addpath('/groups/ag-remy-2/Imaging/AnalysisTools');
addpath(genpath('/groups/ag-remy-2/Imaging/AnalysisTools'))
path = '/groups/ag-remy-2/Imaging';
table = sprintf('%s/%s/Data/Datatable.xlsx',path,user);

%% query data from excel
cnt=1;
num=row;
    %metadata
readfields = sprintf('B%d:X%d',num,num);
[~,~,readline] = xlsread(table,readfields);
setup = readline{1,1};
folderImg = readline{1,2}
stats.genotype= readline{1,end};
splitted=strsplit(folderImg,'_');
stats.animal=splitted{1};
stats.recording=splitted{2};

%data
if strcmp(mode,'dFF')
    loadpath = sprintf('%s/%s/Results/%s/Combined/%s/analysis.mat',path,user,setup,folderImg);
    savepath = [folder '/' folderImg];
elseif strcmp(mode,'C')
    loadpath = sprintf('%s/%s/Results/%s/Combined/%s_C/analysis_C.mat',path,user,setup,folderImg);
    savepath = [folder '/' folderImg '_C'];
elseif strcmp(mode,'S')
    loadpath = sprintf('%s/%s/Results/%s/Combined/%s_S/analysis_S.mat',path,user,setup,folderImg);
    savepath = [folder '/' folderImg '_S'];
else
    error('wrong mode.')
end

%% initialize save path
if ~exist(savepath,'dir') > 0
    try
        mkdir(savepath);
    catch
        error('could not create savedir.')
    end
end

%% read in and plot
if exist(loadpath,'file')
    load(loadpath,'SpeedSlope','SpeedRMSE','ConfintSlopes','SlopeSigs','SpeedCorr','SpeedIQR', ...
        'PlaceScore','PlaceScorePct','centroidBelt','centroidrho','meandFFs','rayleighP','rayleighZ', ...
        'dFFRest','dFFRun','dFFlick','dFFnolick','dFF_out','dFFInBin','vFT','dFF_samp');
    %speed Stats
    stats.dFF_samp=dFF_samp;
    stats.vFT=vFT;
    stats.SpeedSlopes=SpeedSlope;
    stats.SpeedRMSEs=SpeedRMSE;
    stats.ConfintSlopess=ConfintSlopes;
    stats.SlopeSigs=SlopeSigs;
    stats.SpeedCorr=SpeedCorr;
    stats.SpeedIQR=SpeedIQR;
    %place stats
    stats.PlaceScore=PlaceScore;
    stats.PlaceScorePct=PlaceScorePct;
    stats.centroidBelt=centroidBelt;
    stats.centroidrho=centroidrho;
    stats.meandFFs=meandFFs; %to compute placefield widhs
    if strcmp(mode,'S') %% try to fix rayleigh P
        firstmeandFF=meandFFs{1};
        binsize=round(360/length(firstmeandFF(2:end)));        
        for k=1:length(PlaceScore)
            meandFF=meandFFs{k}; %compute new centers!
            meandFF(1)=[]; %first one is wrtong
            centerBins=(binsize/2):binsize:360-(binsize/2);
            centerBins(meandFF==0)=[];
            meandFF(meandFF==0)=[];
            if length(meandFF>0)
                [rayleighP_temp,rayleighZ_temp] = circ_rtest(deg2rad(centerBins),meandFF);
                rayleighP(k)=rayleighP_temp;
                rayleighZ(k)=rayleighZ_temp;
            else
                rayleighP(k)=1;
                rayleighZ(k)=1000;
            end
        end       
    else
        stats.rayleighP=rayleighP;
        stats.rayleighZ=rayleighZ;
    end
    centroidrhoNorm=zeros(1,length(centroidrho));
    for i=1:length(centroidrho)
        centroidrhoNorm(i)=centroidrho(i)/nanmax(dFF_out(i,:),[],2);
    end   
    stats.centroidrho_norm=centroidrhoNorm;
    %other  satats - use dff_samp    
    for k=1:length(PlaceScore)    
        stats.dFFRest(k)=nanmean(dFF_samp(k,vFT<2));
        stats.dFFRun(k)=nanmean(dFF_samp(k,vFT>=2));
    end
    stats.vFT=vFT;
    
    stats.dFFlick=dFFlick;
    stats.dFFnolick=dFFnolick;
    stats.dFF_out=dFF_out;
    
    if plotit
        %generate speed figures
        groupnames=cell(1,length(SlopeSigs));
        for i=1:length(SlopeSigs)
            if SlopeSigs(i)
                groupnames{i}='significant';
            else
                groupnames{i}='not significant';
            end
        end
        scatterhist(SpeedSlope,SpeedCorr,'Group',groupnames','Kernel','on');
        fn=sprintf('%s/speed.png',savepath);
        set(gcf,'Position',[0,0,1200,600])
        print(fn,'-dpng','-r400')
        fn=sprintf('%s/speed.fig',savepath);
        savefig(fn);
        close
        
        %generate place stats figures
        groupnames=cell(1,length(SlopeSigs));
        for i=1:length(SlopeSigs)
            if rayleighP(i)<=0.05 && PlaceScorePct(i)>=0.80
                groupnames{i}='criterion';
            else
                groupnames{i}='no criterion';
            end
        end
        scatterhist(rayleighP,PlaceScorePct,'Group',groupnames','Kernel','on');
        fn=sprintf('%s/place.png',savepath);
        set(gcf,'Position',[0,0,1200,600])
        print(fn,'-dpng','-r400')
        fn=sprintf('%s/place.fig',savepath);
        savefig(fn);
        close
        
        %test centroid as stats
        scatterhist(centroidrho,PlaceScorePct,'Kernel','on');
        fn=sprintf('%s/centroid.png',savepath);
        set(gcf,'Position',[0,0,1200,600])
        print(fn,'-dpng','-r400')
        fn=sprintf('%s/centroid.fig',savepath);
        savefig(fn);
        close
        
    end
    
    %addtional speed stats
    inds=find(SlopeSigs & SpeedCorr<0);
    stats.avgNegCorr = mean(SpeedCorr(inds));
    stats.stdNegCorr = std(SpeedCorr(inds));
    inds=find(SlopeSigs & SpeedCorr>0);
    stats.avgPosCorr = mean(SpeedCorr(inds));
    stats.stdPosCorr = std(SpeedCorr(inds));
    inds=find(SlopeSigs & SpeedSlope<0);
    stats.avgNegSlope = mean(SpeedSlope(inds));
    stats.stdNegSlope = std(SpeedSlope(inds));
    inds=find(SlopeSigs & SpeedSlope>0);
    stats.avgPosSlope = mean(SpeedSlope(inds));
    stats.stdPosSlope = std(SpeedSlope(inds));
    
    %general firing stats
    
    %place statistics
    if isfield(criterion,'minPShuffle') && isfield(criterion,'minRhoNorm')
        placeCells = find(PlaceScorePct>criterion.minPShuffle & centroidrhoNorm>criterion.minRhoNorm);
        noPlaceCells = find(not(PlaceScorePct>criterion.minPShuffle & centroidrhoNorm>criterion.minRhoNorm));
    elseif isfield(criterion,'minPShuffle') && isfield(criterion,'maxRayleighP')
        placeCells = find(PlaceScorePct>criterion.minPShuffle & rayleighP<criterion.maxRayleighP);
        noPlaceCells = find(not(PlaceScorePct>criterion.minPShuffle & rayleighP<criterion.maxRayleighP));
    else
        error('no placecell criterion supplied.')
    end
    
    %general firing statistics
    
    meandFFstats=zeros(1,length(PlaceScore));
    stddFFstats=zeros(1,length(PlaceScore));
    for k=1:length(PlaceScore)    
        meandFFstats(k)=nanmean(dFF_samp(k,:));
        stddFFstats(k)=nanstd(dFF_samp(k,:));
    end
    if strcmp(mode,'S')
       meandFFstats=meandFFstats+0.01;
       stddFFstats=stddFFstats+0.01;
    end
    stats.meandFFstats=meandFFstats;
    stats.stddFFstats=stddFFstats;
    
    [f_mean_all,xi_mean_all]=ksdensity(meandFFstats,'support','positive','function','cdf');
    [f_mean_placecells,xi_mean_placecells]=ksdensity(meandFFstats(placeCells),'support','positive','function','cdf');
    [f_mean_notplacecells,xi_mean_notplacecells]=ksdensity(meandFFstats(noPlaceCells),'support','positive','function','cdf');
    stats.f_mean_all=f_mean_all;
    stats.xi_mean_all=xi_mean_all;
    stats.f_mean_placecells=f_mean_placecells;
    stats.xi_mean_placecells=xi_mean_placecells;
    stats.f_mean_notplacecells=f_mean_notplacecells;
    stats.xi_mean_notplacecells=xi_mean_notplacecells;
    
    
    [f_std_all,xi_std_all]=ksdensity(stddFFstats,'support','positive','function','cdf');
    [f_std_placecells,xi_std_placecells]=ksdensity(stddFFstats(placeCells),'support','positive','function','cdf');
    [f_std_notplacecells,xi_std_notplacecells]=ksdensity(stddFFstats(noPlaceCells),'support','positive','function','cdf');  
    stats.f_std_all=f_std_all;
    stats.xi_std_all=xi_std_all;
    stats.f_std_placecells=f_std_placecells;
    stats.xi_std_placecells=xi_std_placecells;
    stats.f_std_notplacecells=f_std_notplacecells;
    stats.xi_std_notplacecells=xi_std_notplacecells;    
    
    
    if plotit
        subplot(2,3,1)
        plot(xi_mean_all,f_mean_all)
        title('mean all')
        subplot(2,3,2)
        plot(xi_mean_placecells,f_mean_placecells)
        title('mean placecells')        
        subplot(2,3,3)
        plot(xi_mean_notplacecells,f_mean_notplacecells)
        title('mean no placecells')   
        
        subplot(2,3,4)
        plot(xi_std_all,f_std_all)
        title('std all')
        subplot(2,3,5)
        plot(xi_std_placecells,f_std_placecells)
        title('std placecells')        
        subplot(2,3,6)
        plot(xi_std_notplacecells,f_std_notplacecells)
        title('std no placecells')   
        fn=sprintf('%s/activity.png',savepath);
        set(gcf,'Position',[0,0,1200,600])
        print(fn,'-dpng','-r400')
        fn=sprintf('%s/activity.fig',savepath);
        savefig(fn);
        close  
        
        figure
        subplot(1,2,1)
        belt=linspace(0,360,45);
        plotMat=zeros(length(placeCells),45);
        [S,I]=sort(centroidBelt(placeCells));
        for i=1:length(placeCells)
            temp=meandFFs{placeCells(I(i))};
            plotMat(i,:)=temp(2:end)./max(temp(2:end));
        end
        imagesc(belt,1:length(placeCells),plotMat)
        set(gca,'CLIM',[0,prctile(plotMat(:),99)])
        hold on
        plot(sort(centroidBelt(placeCells)),1:length(centroidBelt(placeCells)),'r','LineWidth',2)
        plot([0,360],[1,length(centroidBelt(placeCells))],'w-.','LineWidth',2)
        plot([36,36],[0,length(centroidBelt(placeCells))],'g:','LineWidth',2)
        plot([44,44],[0,length(centroidBelt(placeCells))],'g:','LineWidth',2)
        xlabel('Belt position (cm)')
        title('place cells')
        subplot(1,2,2)
        plotMat=zeros(length(noPlaceCells),45);
        [S,I]=sort(centroidBelt(noPlaceCells));
        for i=1:length(noPlaceCells)
            temp=meandFFs{noPlaceCells(I(i))};
            plotMat(i,:)=temp(2:end)./max(temp(2:end));
        end
        imagesc(belt,1:length(noPlaceCells),plotMat)
        set(gca,'CLIM',[0,prctile(plotMat(:),99)])
        hold on
        plot(sort(centroidBelt(noPlaceCells)),1:length(centroidBelt(noPlaceCells)),'r','LineWidth',2)
        plot([0,360],[1,length(centroidBelt(noPlaceCells))],'w-.','LineWidth',2)
        plot([36,36],[0,length(centroidBelt(noPlaceCells))],'g:','LineWidth',2)
        plot([44,44],[0,length(centroidBelt(noPlaceCells))],'g:','LineWidth',2)
        xlabel('Belt position (cm)')
        title('no place cells')
        fn=sprintf('%s/plotMat.png',savepath);
        set(gcf,'Position',[0,0,1200,600])
        print(fn,'-dpng','-r400')
        fn=sprintf('%s/plotMat.fig',savepath);
        savefig(fn);
        close  
    end
    
    
    [f,xi]=ksdensity(centroidrho(placeCells),'support','positive','function','cdf');
    stats.density_f=f;
    stats.density_xi=xi;
    
    fieldwidth=[];
    mean_dff_field=[];
    std_dff_field=[];
    mean_dff_offield=[];
    std_dff_offield=[];
    pad=10;
    for i=1:length(placeCells)
            temp=meandFFs{placeCells(i)};
            temp(1)=[];
            temp=smooth(temp)';
            [m,ind]=max(temp);
            temp_padded=[temp(end-pad:end) temp temp(1:pad)];
            ind=ind+pad;
            field_start=findonset_backward(temp_padded,ind,pad);
            field_end=findonset_forward(temp_padded,ind,pad);
            fieldwidth(end+1)=field_end-field_start;   
            muh=squeeze(dFFInBin.ALL(i,:,:));
            binfield=zeros(1,pad+45+pad);
            binfield(field_start:field_end)=1;
            binfield_shift=circshift(binfield',-pad+1)';
            binfield_shift(47:end)=[];
            
            muhuh=nanmean(muh(:,binfield_shift>0),2);
            mean_dff_field(end+1)=nanmean(muhuh);
            muhuh=nanstd(muh(:,binfield_shift>0),[],2);
            std_dff_field(end+1)=nanstd(muhuh);            
            
            %invert mask
            binfield_shift=~binfield_shift;
            binfield_shift(1)=0;            
            muhuh=nanmean(muh(:,binfield_shift>0),2);
            mean_dff_offield(end+1)=nanmean(muhuh);
            muhuh=nanstd(muh(:,binfield_shift>0),[],2);
            std_dff_offield(end+1)=nanstd(muhuh);  
    end    
    stats.fieldwidth=fieldwidth;
    stats.avgfieldwidth=mean(fieldwidth);
    stats.mean_dff_field=mean_dff_field;
    stats.std_dff_field=std_dff_field;
    stats.mean_dff_offield=mean_dff_offield;
    stats.std_dff_offield=std_dff_offield;    
        
    if plotit        
        %test centroid as stats
        groupnames=cell(1,length(SlopeSigs));
        for i=1:length(SlopeSigs)
            if any(placeCells==i)
                groupnames{i}='criterion';
            else
                groupnames{i}='no criterion';
            end
        end
        scatterhist(centroidrhoNorm,PlaceScorePct,'Group',groupnames','Kernel','on');        
        fn=sprintf('%s/centroid_norm.png',savepath);
        set(gcf,'Position',[0,0,1200,600])
        print(fn,'-dpng','-r400')
        fn=sprintf('%s/centroid_norm.fig',savepath);
        savefig(fn);
        close        
        
        figure
        subplot(1,2,1)
        plot(sort(centroidBelt(placeCells)))
        hold on
        plot([0,length(centroidBelt(placeCells))],[0,360],'r--')
        subplot(1,2,2)
        plot(f,xi)
        fn=sprintf('%s/placecells.png',savepath);
        set(gcf,'Position',[0,0,1200,600])
        print(fn,'-dpng','-r400')
        fn=sprintf('%s/placecells.fig',savepath);
        savefig(fn);
        close               
        
        figure
        subplot(2,2,1)
        hist(mean_dff_field)
        title('Mean activity field')
        subplot(2,2,2)
        hist(std_dff_field)
        title('Std activity field')
        subplot(2,2,3)
        hist(mean_dff_offield)
        title('Mean activity off-field')
        subplot(2,2,4)
        hist(std_dff_offield)
        title('Std activity off-field')
        
        fn=sprintf('%s/placefields.png',savepath);
        set(gcf,'Position',[0,0,1200,600])
        print(fn,'-dpng','-r400')
        fn=sprintf('%s/placefields.fig',savepath);
        savefig(fn);
        close 
    end
    
    %rest/run lick/nolick statistics
    stats.dFFRestRunScore=dFFRun./(dFFRest+dFFRun);
    % FIX!!
    stats.dFFlickScore=dFFlick./(dFFlick+dFFnolick);
    %disp(dFFnolick)
    
else
    disp('Analysis file not found.');
    fieldStats=NaN;
end
end
