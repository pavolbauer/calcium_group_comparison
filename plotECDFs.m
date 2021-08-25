function plotECDFs(stats,xq,eventfix,mode,plottraces)
% P. Bauer 2020
genotypes=[stats.genotype];
f_interps_tg=[];
f_interps_wt=[];
for i=1:length(genotypes)
    if genotypes(i)
        if strcmp(mode,'mean_dff_field')
            [f,xi]=ksdensity(stats(i).mean_dff_field(stats(i).mean_dff_field>0),'support','positive','function','cdf');
        elseif strcmp(mode,'mean_dff_placecells')
            meanplace=stats(i).PlaceScore(stats(i).placeCells);
            [f,xi]=ksdensity(meanplace(meanplace>0),'support','positive','function','cdf');
        elseif strcmp(mode,'std_dff_field')
            [f,xi]=ksdensity(stats(i).std_dff_field(stats(i).std_dff_field>0),'support','positive','function','cdf');
        elseif strcmp(mode,'mean_dff_offield')
            [f,xi]=ksdensity(stats(i).mean_dff_offield(stats(i).mean_dff_offield>0),'support','positive','function','cdf');
        elseif strcmp(mode,'std_dff_offield')
            [f,xi]=ksdensity(stats(i).std_dff_offield(stats(i).std_dff_offield>0),'support','positive','function','cdf');
        elseif strcmp(mode,'mean_stability')
            [f,xi]=ksdensity(stats(i).stabs(stats(i).stabs>0),'support','positive','function','cdf');
        elseif strcmp(mode,'placescore')
            [f,xi]=ksdensity(stats(i).PlaceScore(stats(i).PlaceScore>0),'support','positive','function','cdf'); 
        elseif strcmp(mode,'placescorepct')
            [f,xi]=ksdensity(stats(i).PlaceScorePct(stats(i).PlaceScorePct>0),'support','positive','function','cdf');             
        elseif strcmp(mode,'mean_fieldwidths')
            [f,xi]=ksdensity(stats(i).fieldwidths(stats(i).fieldwidths>0),'support','positive','function','cdf');            
        elseif strcmp(mode,'mean_posspeed')
            signi=stats(i).SlopeSigs;
            SpeedCorr=stats(i).SpeedCorr;
            distribution=[];
            for k=1:length(SpeedCorr)
                if signi(k) && SpeedCorr(k)>0
                    distribution(end+1)=SpeedCorr(k);
                end
            end
            [f,xi]=ksdensity(distribution(distribution>0),'support','positive','function','cdf');    
        elseif strcmp(mode,'mean_negspeed')
            signi=stats(i).SlopeSigs;
            SpeedCorr=stats(i).SpeedCorr;
            distribution=[];
            for k=1:length(SpeedCorr)
                if signi(k) && SpeedCorr(k)<0
                    distribution(end+1)=SpeedCorr(k);
                end
            end
            if ~isempty(distribution)
                [f,xi]=ksdensity(distribution(distribution>0),'support','positive','function','cdf');              
            else
                xi=0;
                f=0;
            end
        elseif strcmp(mode,'cv_dff_field')
            stdvec=stats(i).std_dff_field;
            meanvec=stats(i).mean_dff_field;
            cv_dff_field=(stdvec./meanvec);
            [f,xi]=ksdensity(cv_dff_field(cv_dff_field>0),'support','positive','function','cdf');
        elseif strcmp(mode,'cv_dff_offield')
            stdvec=stats(i).std_dff_offield;
            meanvec=stats(i).mean_dff_offield;
            cv_dff_offield=(stdvec./meanvec);
            [f,xi]=ksdensity(cv_dff_offield(cv_dff_offield>0),'support','positive','function','cdf');
        end
        f_interp=interp1(xi,f,xq);
        f_interps_tg=[f_interps_tg; f_interp];
        hold on
    else
        if strcmp(mode,'mean_dff_field')
            [f,xi]=ksdensity(stats(i).mean_dff_field(stats(i).mean_dff_field>0),'support','positive','function','cdf');
        elseif strcmp(mode,'mean_dff_placecells')
            meanplace=stats(i).PlaceScore(stats(i).placeCells);
            [f,xi]=ksdensity(meanplace(meanplace>0),'support','positive','function','cdf');
        elseif strcmp(mode,'std_dff_field')
            [f,xi]=ksdensity(stats(i).std_dff_field(stats(i).std_dff_field>0),'support','positive','function','cdf');
        elseif strcmp(mode,'mean_dff_offield')
            [f,xi]=ksdensity(stats(i).mean_dff_offield(stats(i).mean_dff_offield>0),'support','positive','function','cdf');
        elseif strcmp(mode,'std_dff_offield')
            [f,xi]=ksdensity(stats(i).std_dff_offield(stats(i).std_dff_offield>0),'support','positive','function','cdf');    
        elseif strcmp(mode,'mean_stability')
            [f,xi]=ksdensity(stats(i).stabs(stats(i).stabs>0),'support','positive','function','cdf');     
        elseif strcmp(mode,'placescore')
            [f,xi]=ksdensity(stats(i).PlaceScore(stats(i).PlaceScore>0),'support','positive','function','cdf'); 
        elseif strcmp(mode,'placescorepct')
            [f,xi]=ksdensity(stats(i).PlaceScorePct(stats(i).PlaceScorePct>0),'support','positive','function','cdf');             
        elseif strcmp(mode,'mean_fieldwidths')
            [f,xi]=ksdensity(stats(i).fieldwidths(stats(i).fieldwidths>0),'support','positive','function','cdf');     
        elseif strcmp(mode,'mean_posspeed')
            signi=stats(i).SlopeSigs;
            SpeedCorr=stats(i).SpeedCorr;
            distribution=[];
            for k=1:length(SpeedCorr)
                if signi(k) && SpeedCorr(k)>0
                    distribution(end+1)=SpeedCorr(k);
                end
            end
            [f,xi]=ksdensity(distribution(distribution>0),'support','positive','function','cdf');    
        elseif strcmp(mode,'mean_negspeed')
            signi=stats(i).SlopeSigs;
            SpeedCorr=stats(i).SpeedCorr;
            distribution=[];
            for k=1:length(SpeedCorr)
                if signi(k) && SpeedCorr(k)<0
                    distribution(end+1)=SpeedCorr(k);
                end
            end
            [f,xi]=ksdensity(distribution(distribution>0),'support','positive','function','cdf');             
        elseif strcmp(mode,'cv_dff_field')
            stdvec=stats(i).std_dff_field;
            meanvec=stats(i).mean_dff_field;
            cv_dff_field=(stdvec./meanvec);
            [f,xi]=ksdensity(cv_dff_field(cv_dff_field>0),'support','positive','function','cdf');
        elseif strcmp(mode,'cv_dff_offield')
            stdvec=stats(i).std_dff_offield;
            meanvec=stats(i).mean_dff_offield;
            cv_dff_offield=(stdvec./meanvec);
            [f,xi]=ksdensity(cv_dff_offield(cv_dff_offield>0),'support','positive','function','cdf');         
        end
        f_interp=interp1(xi,f,xq);
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

tgcnt=1;
wtcnt=1;
if plottraces
    for i=1:length(genotypes)
        if genotypes(i)
            h1=plot(xq.*eventfix,f_interps_tg(tgcnt,:),'r');
            tgcnt=tgcnt+1;
        else
            h2=plot(xq.*eventfix,f_interps_wt(wtcnt,:),'b');
            wtcnt=wtcnt+1;
        end    
    end
end

h1=plot(xq.*eventfix,mean_act_tg,'r','LineWidth',4);
h2=plot(xq.*eventfix,mean_act_wt,'b','LineWidth',4);

hline(0.3,'--k')
hline(0.8,'--k')

legend([h1(1),h2(1)],'TG','WT')
ylim([0 1.1])