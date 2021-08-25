function [placeCells,noPlaceCells,posSpeedCells,negSpeedCells,notClear,centroidrhoNorm,PlaceScorePct]=readPlaceCells(user,row,mode,folder,plotit,criterion)
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

%% read in and plot
if exist(loadpath,'file')
    load(loadpath,'SpeedSlope','SpeedRMSE','ConfintSlopes','SlopeSigs','SpeedCorr','SpeedIQR', ...
        'PlaceScore','PlaceScorePct','centroidBelt','centroidrho','meandFFs','rayleighP','rayleighZ', ...
        'dFFRest','dFFRun','dFFlick','dFFnolick','dFF_out','dFFInBin','vFT','dFF_samp');

        centroidrhoNorm=zeros(1,length(centroidrho));
    for i=1:length(centroidrho)
        centroidrhoNorm(i)=centroidrho(i)/nanmax(dFF_out(i,:),[],2);
    end 
    
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
        %print(fn,'-dpng','-r400')
        %fn=sprintf('%s/speed.fig',savepath);
        %savefig(fn);
        %close
        
        figure
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
        %print(fn,'-dpng','-r400')
        %fn=sprintf('%s/place.fig',savepath);
        %savefig(fn);
        %close
        
        %test centroid as stats
        %scatterhist(centroidrho,PlaceScorePct,'Kernel','on');
        %fn=sprintf('%s/centroid.png',savepath);
        %set(gcf,'Position',[0,0,1200,600])
        %print(fn,'-dpng','-r400')
        %fn=sprintf('%s/centroid.fig',savepath);
        %savefig(fn);
        %close
        
    end
        
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
    
    %get speed stats for non-place cells
    posSpeedCells=[];
    negSpeedCells=[];
    notClear=[];
    for i=1:length(noPlaceCells)
        if SlopeSigs(noPlaceCells(i))
            if SpeedCorr(noPlaceCells(i))>0
                posSpeedCells(end+1)=noPlaceCells(i);
            else
                negSpeedCells(end+1)=noPlaceCells(i);
            end
        else
            notClear(end+1)=noPlaceCells(i);
        end
    end
    
else
    disp('Analysis file not found.');
    fieldStats=NaN;
end
end
