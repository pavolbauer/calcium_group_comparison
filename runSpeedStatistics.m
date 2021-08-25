function speedStats=runSpeedStatistics(user,row,mode)
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
folderIgor = readline{1,8};
wavenum = readline{1,9};
genotypes(cnt)= readline{1,end};
cut = readline{1,end-1};
seg_perc = readline{1,end-2};    
splitted=strsplit(folderImg,'_');
if strcmp(splitted{2},'SR')
    SROLM(cnt)=1;
else
    SROLM(cnt)=0;
end
fieldStats.animal=splitted{1};
fieldStats.recording=splitted{2};

%data
if strcmp(mode,'dFF')
    loadpath = sprintf('%s/%s/Results/%s/Combined/%s/analysis.mat',path,user,setup,folderImg);
elseif strcmp(mode,'C')
    loadpath = sprintf('%s/%s/Results/%s/Combined/%s_C/analysis_C.mat',path,user,setup,folderImg);
elseif strcmp(mode,'S')
    loadpath = sprintf('%s/%s/Results/%s/Combined/%s_S/analysis_S.mat',path,user,setup,folderImg);
else
    error('wrong mode.')
end

if exist(loadpath,'file')
    load(loadpath,'SpeedSlope','SpeedRMSE','ConfintSlopes','SlopeSigss','SpeedCorr','SpeedIQR','dFF_out')
    SpeedSlopes{cnt}=SpeedSlope;
    SpeedRMSEs{cnt}=SpeedRMSE;
    ConfintSlopess{cnt}=ConfintSlopes;
    SlopeSigss{cnt}=SlopeSigss;
    SpeedCorrs{cnt}=SpeedCorr;
    SpeedIQRs{cnt}=SpeedIQR;
else
    disp('Analysis file not found.');
    fieldStats=NaN;
end
end
