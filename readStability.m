function [mean_dff_field,std_dff_field,mean_dff_offield,std_dff_offield,stabs,binfields,fieldwidths,placeCells,noPlaceCells]=readStability(user,row,mode,folder,plotit,criterion)
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
    
    pad=10;
    mean_dff_field=[];
    std_dff_field=[];
    mean_dff_offield=[];
    std_dff_offield=[];
    stabs=[];
    binfields=[];
    fieldwidths=[];
    
    for i=1:length(placeCells)
            
            muh=squeeze(dFFInBin.ALL(placeCells(i),:,:));
            temp=nanmean(muh,1);
            %temp2=meandFFs{placeCells(i)};
            temp(1)=[];
            temp=smooth(temp,7)';
            [m,ind]=max(temp);
            temp_padded=[temp(end-pad:end) temp temp(1:pad)];
            ind=ind+pad;
            field_start=findonset_backward(temp_padded,ind+1,pad);
            field_end=findonset_forward(temp_padded,ind+1,pad);
            fieldwidths(end+1)=field_end-field_start;   
                  
                        
            binfield=zeros(1,pad+45+pad);
            binfield(field_start:field_end)=1;                       
            binfield_shift=circshift(binfield',-pad+1)';
            binfield_shift(47:end)=[];
            binfields(i,:)=binfield_shift;
            
            stability=0;
            for j=1:size(muh,1)
                [mx,ix]=max(muh(j,:));
                if ix>(field_start-pad) && ix>(field_end-pad)
                    stability=stability+1;
                end
            end
            stabs(end+1)=stability/size(muh,1);
            
            fieldmean=nanmean(muh(:,binfield_shift>0),2);
            fieldmeanmean=nanmean(fieldmean);
            mean_dff_field(end+1)=fieldmeanmean;
            fieldstd=nanstd(muh(:,binfield_shift>0),[],2);
            fieldstdmean=nanmean(fieldstd);
            std_dff_field(end+1)=fieldstdmean;            
            
            %invert mask
            binfield_shift=~binfield_shift;
            binfield_shift(1)=0;            
            muhuh=nanmean(muh(:,binfield_shift>0),2);
            mean_dff_offield(end+1)=nanmean(muhuh);
            muhuh=nanstd(muh(:,binfield_shift>0),[],2);
            std_dff_offield(end+1)=nanmean(muhuh);   
    end    
    
else
    disp('Analysis file not found.');
    fieldStats=NaN;
end
end
