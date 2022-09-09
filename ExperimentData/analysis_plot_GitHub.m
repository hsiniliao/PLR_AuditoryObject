
clear

data_cond={'E1Norm.mat','E2SNorm.mat','E2GNorm.mat','E3Norm.mat','E4Norm.mat'};
name_cond = {{'Experiment 1:' 'Mathôt et al. (2013)''s paradigm'},{'Experiment 2:' 'Attend-to-Location'},{'Experiment 2:' 'Attend-to-Gender'},{'Experiment 3:' 'Mixed with diotic trials'},{'Experiment 4:' 'Loudspeakers'}};
exp_cond={'1','CRM','CRM','CRM','CRM'};

main_plot='on';
switch main_plot
    case 'on'
        for iFig=1:6
            figure(iFig),clf
        end
        h_p=NaN(length(data_cond),1);
        load('JackKnife.mat')
        latencyMean = mean(latency_matrix,'omitnan');
        latencyStd = NaN(size(latencyMean));
    case 'off'
        for iFig=3:4
            figure(iFig),clf
        end
end

iCond=1;
for iCond=1:length(data_cond)

%% Load data
filename=split(data_cond{iCond},'.');
load(data_cond{iCond});
pMatrix = readmatrix(strcat('tables/', filename{1}, '_p.csv'));
gLumMatrix = readmatrix(strcat('tables/', filename{1}, '_gLum.csv'));

variables = {'display','cue','target','subj'};

%% Setup color code and label
cmap1=[0 0 1;1 0 0;0 0 0];
ltype1={'-',':'};
cmap2=[0.8 0.8 1;1 0.8 0.8;0.8 0.8 0.8];
cmap3=[0 0 1;0 0 1;1 0 0;1 0 0];
cmap4=[0.8 0.8 1;0.8 0.8 1;1 0.8 0.8;1 0.8 0.8];
ltype2={'-',':',':','-'};
cmapLines=colormap('lines');

sig_cond = [-0.2 -0.05 -0.05 0.05 -0.05];
plot_FontSize_main = 24;
plot_FontSize_group = 14;
plotPos_cond = [1,3,2,4,6];

DisplayName = {'Black-Left','Black-Right'};
if iCond==4
    CueName = {'Left','Right','Center'};
else
    CueName = {'Left','Right'};
end
LumName = {'dark','bright'};

cue_number=length(CueName);
lum_number=length(LumName);

switch exp_cond{iCond}
    case '1'
        xlim_range = [-0.5 5];
    case 'CRM'
        xlim_range = [-0.5 7];
end

if strcmp(main_plot,'on')
%% Plot pupil result across time
tmpMatrix=NaN(length(t_trial),lum_number,length(data));

for iLum=1:lum_number
    for iSubj=1:length(data)
        I=intersect(find(pMatrix(:,length(t_trial)+length(variables)+1)==iLum),find(pMatrix(:,length(t_trial)+length(variables))==iSubj));
        tmpMatrix(:,iLum,iSubj)=mean(pMatrix(I,1:length(t_trial)),'omitnan');
    end
end

type=mean(tmpMatrix,3,'omitnan');
typeStd=std(tmpMatrix,0,3,'omitnan')/sqrt(length(data));

cond1_p = squeeze(tmpMatrix(:,1,:));
cond2_p = squeeze(tmpMatrix(:,2,:));

cfg = [];
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.numrandomization = 5000;
cfg.correctm         = 'cluster';
cfg.method           = 'montecarlo';
cfg.tail             = 0;
cfg.alpha            = 0.05;
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsize';
cfg.design           = [1:length(data) 1:length(data) % subject number
    ones(1,length(data)) 2*ones(1,length(data))];  % condition number
cfg.uvar = 1;        % "subject" is unit of observation
cfg.ivar = 2;        % "condition" is the independent variable
cfg.dimord = 'time';
cfg.dim=[1,numel(t_trial)];
cfg.connectivity =1;
cfg.spmversion = 'spm12';

stat = ft_statistics_montecarlo(cfg, [cond1_p cond2_p],cfg.design);

% Find indices of significant clusters
pos_p=[];
if isfield(stat,'posclusters')
    if ~isempty(stat.posclusters)
        pos_cluster_pvals = [stat.posclusters(:).prob];
        pos_signif_clust = find(pos_cluster_pvals < cfg.alpha);
        poss = ismember(stat.posclusterslabelmat, pos_signif_clust);
        if size(find(diff([0; poss])==-1),1) ~= size(find(diff([0; poss])==1),1)
            aa = [find(diff([0; poss])==-1); length(poss)];
            pos_p = [find(diff([0; poss])==1) aa];
        else
            pos_p = [find(diff([0; poss])==1) find(diff([0; poss])==-1)];
        end
    end
end

if iCond==4
    cond3=NaN(length(t_trial),length(data));
    for iSubj=1:length(data)
        I=intersect(find(pMatrix(:,length(t_trial)+2)==3),find(pMatrix(:,length(t_trial)+length(variables))==iSubj));
        cond3(:,iSubj)=mean(pMatrix(I,1:length(t_trial)),'omitnan');
    end

    typeCond3 = mean(cond3,2,'omitnan');
    typeStdCond3 = std(cond3,0,2,'omitnan')/sqrt(length(data));
    
    stat = ft_statistics_montecarlo(cfg, [cond1_p cond3],cfg.design);
    
    % Find indices of significant clusters
    pos13=[];
    if isfield(stat,'posclusters')
        if ~isempty(stat.posclusters)
            pos_cluster_pvals = [stat.posclusters(:).prob];
            pos_signif_clust = find(pos_cluster_pvals < cfg.alpha);
            poss = ismember(stat.posclusterslabelmat, pos_signif_clust);
            if size(find(diff([0; poss])==-1),1) ~= size(find(diff([0; poss])==1),1)
                aa = [find(diff([0; poss])==-1); length(poss)];
                pos13 = [find(diff([0; poss])==1) aa];
            else
                pos13 = [find(diff([0; poss])==1) find(diff([0; poss])==-1)];
            end
        end
    end
    
    stat = ft_statistics_montecarlo(cfg, [cond3 cond2_p],cfg.design);
    
    % Find indices of significant clusters
    pos32=[];
    if isfield(stat,'posclusters')
        if ~isempty(stat.posclusters)
            pos_cluster_pvals = [stat.posclusters(:).prob];
            pos_signif_clust = find(pos_cluster_pvals < cfg.alpha);
            poss = ismember(stat.posclusterslabelmat, pos_signif_clust);
            if size(find(diff([0; poss])==-1),1) ~= size(find(diff([0; poss])==1),1)
                aa = [find(diff([0; poss])==-1); length(poss)];
                pos32 = [find(diff([0; poss])==1) aa];
            else
                pos32 = [find(diff([0; poss])==1) find(diff([0; poss])==-1)];
            end
        end
    end
end

figure(1),clf
hold on
h=NaN(lum_number,1);
for iLum=1:lum_number
    patch([t_trial(1:end-1);t_trial(1:end-1);t_trial(2:end);t_trial(2:end)],...
        [type(1:end-1,iLum)+typeStd(1:end-1,iLum),...
        type(1:end-1,iLum)-typeStd(1:end-1,iLum),...
        type(2:end,iLum)-typeStd(2:end,iLum),...
        type(2:end,iLum)+typeStd(2:end,iLum)]',...
        ones(4,length(t_trial)-1),'FaceColor',cmap1(iLum,:),'EdgeColor','none','FaceAlpha',0.1);
    h(iLum)=plot(t_trial,type(:,iLum),'Color',cmap1(iLum,:),'LineWidth',3);
end
for i = 1:size(pos_p,1);line([t_trial(pos_p(i,1)) t_trial(pos_p(i,2))],[sig_cond(iCond) sig_cond(iCond)],'LineWidth',3,'Color',[0 0 0]); end
if iCond==4
    patch([t_trial(1:end-1);t_trial(1:end-1);t_trial(2:end);t_trial(2:end)],...
        [typeCond3(1:end-1)+typeStdCond3(1:end-1),...
        typeCond3(1:end-1)-typeStdCond3(1:end-1),...
        typeCond3(2:end)-typeStdCond3(2:end),...
        typeCond3(2:end)+typeStdCond3(2:end)]',...
        ones(4,length(t_trial)-1),'FaceColor',cmap1(3,:),'EdgeColor','none','FaceAlpha',0.1);
    h(3)=plot(t_trial,typeCond3,'Color',cmap1(3,:),'LineWidth',3);
    for i = 1:size(pos13,1);line([t_trial(pos13(i,1)) t_trial(pos13(i,2))],[sig_cond(iCond)+0.2 sig_cond(iCond)+0.2],'LineWidth',3,'Color',cmap1(1,:)); end
    for i = 1:size(pos32,1);line([t_trial(pos32(i,1)) t_trial(pos32(i,2))],[sig_cond(iCond)+0.1 sig_cond(iCond)+0.1],'LineWidth',3,'Color',cmap1(2,:)); end
end
hold off
if iCond==1
    legend(h,strcat('Attend-',LumName, '-side ear'),'Location','NorthWest')
elseif iCond==2
    legend(h,strcat('Attend-',LumName, '-side ear'),'Location','SouthWest')
elseif iCond==3
    legend(h,strcat('Target-at-',LumName, '-side ear'),'Location','SouthWest')
elseif iCond==4
    legend(h,cat(2,strcat('Target-at-',LumName, '-side ear'),'Diotic trials'),'Location','SouthWest')
elseif iCond==5
    legend(h,strcat('Target-at-',LumName, '-side speaker'),'Location','SouthWest')
end
legend boxoff
xlim(xlim_range)
if iCond==1
    ylim([-0.5 1])
else
    ylim([-0.5 1.2])
end
xlabel('Time from cue (s)')
ylabel('Norm. pupil diameter [z]');
set(gca,'FontSize',plot_FontSize_main)
set(gcf,'Color',[1 1 1])

%% Plot mean pupil result
if iCond==1
    aveRange_pupil = [0 5];
else
    aveRange_pupil = [1 7];
end
timeRange_pupil = (aveRange_pupil(1)-trange(1))*Fs+1:(aveRange_pupil(2)-trange(1))*Fs;

tmpMatrix=NaN(lum_number,length(data));

for iLum=1:lum_number
    for iSubj=1:length(data)
        I=intersect(find(pMatrix(:,length(t_trial)+length(variables)+1)==iLum),find(pMatrix(:,length(t_trial)+length(variables))==iSubj));
        tmpMatrix(iLum,iSubj)=mean2(pMatrix(I,timeRange_pupil));
    end
end

if iCond==4
    tmpMatrix(3,:)=tmpMatrix(2,:);
    for iSubj=1:length(data)
        I=intersect(find(pMatrix(:,length(t_trial)+2)==3),find(pMatrix(:,length(t_trial)+length(variables))==iSubj));
        tmpMatrix(2,iSubj)=mean2(pMatrix(I,timeRange_pupil));
    end
end

Y = mean(tmpMatrix,2);
YStd = std(tmpMatrix,0,2)/sqrt(length(data));

if iCond==4
    plot_color_order=[1 3 2];
else
    plot_color_order=[1 2];
end

figure(2),clf
hold on
for i=1:length(Y)
    h = bar(i,Y(i));
    set(h,'FaceColor',cmap1(plot_color_order(i),:),'LineWidth',1,'EdgeColor','none');
    errorbar(i,Y(i),YStd(i),'Color',cmap1(plot_color_order(i),:),'LineWidth',2);
end
for iSubj=1:length(data)
    plot(1:size(tmpMatrix,1),tmpMatrix(:,iSubj),'Color',[0.5 0.5 0.5],'LineWidth',1,'Marker','none','MarkerFaceColor','k','MarkerSize',3)
end
hold off

xlim([0 length(Y)+1]);
if iCond==1
    ylim([-0.5 1])
else
    ylim([-0.5 1.5]);
end
if iCond==4
    xtick = text(1:length(Y),-0.5*ones(1,length(Y)),{'Dark','Diotic','Bright'});
else
    xtick = text(1:length(Y),-0.5*ones(1,length(Y)),{'Dark','Bright'});
end
set(xtick,'HorizontalAlignment','center','VerticalAlignment','top','FontSize',plot_FontSize_main);
set(gca,'xticklabel','');
ylabel(['Mean pupil diameter [z] of ' num2str(aveRange_pupil(1)) '–' num2str(aveRange_pupil(end)) ' s']);
set(gca,'FontSize',plot_FontSize_main)
set(gcf,'Color',[1 1 1])

tmpStat=tmpMatrix';

end
%% Plot gaze-luminance result segregated by luminance
tmpMatrix=NaN(length(t_trial),lum_number,length(data));

for iLum=1:lum_number
    for iSubj=1:length(data)
        I=intersect(find(gLumMatrix(:,length(t_trial)+length(variables)+1)==iLum),find(gLumMatrix(:,length(t_trial)+length(variables))==iSubj));
        tmpMatrix(:,iLum,iSubj)=smooth(mean(gLumMatrix(I,1:length(t_trial)),'omitnan'),100);
    end
end

type=mean(tmpMatrix,3,'omitnan');
typeStd=std(tmpMatrix,0,3,'omitnan')/sqrt(length(data));

cond1 = squeeze(tmpMatrix(:,1,:));
cond2 = squeeze(tmpMatrix(:,2,:));

cfg = [];
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.numrandomization = 5000;
cfg.correctm         = 'cluster';
cfg.method           = 'montecarlo';
cfg.tail             = 0;
cfg.alpha            = 0.05;
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsize';
cfg.design           = [1:length(data) 1:length(data) % subject number
    ones(1,length(data)) 2*ones(1,length(data))];  % condition number
cfg.uvar = 1;        % "subject" is unit of observation
cfg.ivar = 2;        % "condition" is the independent variable
cfg.dimord = 'time';
cfg.dim=[1,numel(t_trial)];
cfg.connectivity =1;
cfg.spmversion = 'spm12';

stat = ft_statistics_montecarlo(cfg, [cond2 cond1],cfg.design);

% Find indices of significant clusters
pos=[];
if isfield(stat,'posclusters')
    if ~isempty(stat.posclusters)
        pos_cluster_pvals = [stat.posclusters(:).prob];
        pos_signif_clust = find(pos_cluster_pvals < cfg.alpha);
        poss = ismember(stat.posclusterslabelmat, pos_signif_clust);
        if size(find(diff([0; poss])==-1),1) ~= size(find(diff([0; poss])==1),1)
            aa = [find(diff([0; poss])==-1); length(poss)];
            pos = [find(diff([0; poss])==1) aa];
        else
            pos = [find(diff([0; poss])==1) find(diff([0; poss])==-1)];
        end
    end
end

figure(3)
subplot(3,2,plotPos_cond(iCond))
hold on
h=NaN(lum_number,1);
for iLum=1:lum_number
    patch([t_trial(1:end-1);t_trial(1:end-1);t_trial(2:end);t_trial(2:end)],...
        [type(1:end-1,iLum)+typeStd(1:end-1,iLum),...
        type(1:end-1,iLum)-typeStd(1:end-1,iLum),...
        type(2:end,iLum)-typeStd(2:end,iLum),...
        type(2:end,iLum)+typeStd(2:end,iLum)]',...
        ones(4,length(t_trial)-1),'FaceColor',cmap1(iLum,:),'EdgeColor','none','FaceAlpha',0.1);
    h(iLum)=plot(t_trial,type(:,iLum),'Color',cmap1(iLum,:),'LineWidth',3);
end
for i = 1:size(pos,1);line([t_trial(pos(i,1)) t_trial(pos(i,2))],[0.33 0.33],'LineWidth',3,'Color',[0 0 0]); end
hold off
legend(h,strcat('Attend-',LumName),'Location','NorthWest')
legend boxoff
xlim(xlim_range)
ylim([0.25 0.75])
xlabel('Time from cue (s)')
ylabel('Gazed luminance contrast');
title(name_cond{iCond});
set(gca,'FontSize',plot_FontSize_group)
set(gcf,'Color',[1 1 1])

%% Plot linear mixed model result
tArray_gLum = readmatrix(strcat('tables/', filename{1}, '_p_tArray_gLum_GitHub.csv'));
tArray_lum = readmatrix(strcat('tables/', filename{1}, '_p_tArray_lum_GitHub.csv'));

tArray_matrix = cat(2,tArray_gLum(:,2),tArray_lum(:,2));
tArray_matrix = tArray_matrix*-1;
tArrayName = {'Gazed luminance','Attended luminance'};

sig_matrix=cell(size(tArray_matrix,2),1);
for iVar=1:size(tArray_matrix,2)
    tmp=zeros(length(t_trial),1);
    tmp((abs(tArray_matrix(:,iVar))>2))=1;
    if size(find(diff([0; tmp])==-1),1) ~= size(find(diff([0; tmp])==1),1)
        aa = [find(diff([0; tmp])==-1); length(tmp)];
        sig = [find(diff([0; tmp])==1) aa];
    else
        sig = [find(diff([0; tmp])==1) find(diff([0; tmp])==-1)];
    end
    sig(:,3)=sig(:,2)-sig(:,1);
    sig(sig(:,3)<500,:)=[];
    
    sig_matrix{iVar}=sig;
end

figure(4)
subplot(3,2,plotPos_cond(iCond))
h=NaN(length(tArrayName),1);
hold on
for iVar=1:length(tArrayName)
    h(iVar)=plot(t_trial,smooth(tArray_matrix(:,iVar),200),'Color',cmapLines(iVar,:),'LineWidth',3);
    sig=sig_matrix{iVar};
    for i = 1:size(sig,1);line([t_trial(sig(i,1)) t_trial(sig(i,2))],[-4+(iVar-1)*1 -4+(iVar-1)*1],'LineWidth',2,'Color',cmapLines(iVar,:)); end
end
hold off
legend(h,tArrayName,'Location','SouthWest')
legend boxoff
xlim(xlim_range)
ylim([-10 10])
% xticks('')
xlabel('Time from cue (s)')
ylabel('t value');
title(name_cond{iCond});
set(gca,'FontSize',plot_FontSize_group)
set(gcf,'Color',[1 1 1])

if strcmp(main_plot,'on')
%% Plot summary result
type = mean(cond1_p,2,'omitnan')-mean(cond2_p,2,'omitnan');
if iCond==1
    type(5001:end)=NaN;
    pos_p(pos_p>5000)=5000;
end

figure(5)
hold on
h_p(iCond)=plot(t_trial,type,'Color',cmapLines(iCond,:),'LineWidth',3);
hold off

latencyStd(iCond)=std(latency_matrix(:,iCond),'omitnan')/sqrt(length(find(isfinite(latency_matrix(:,iCond)))));
figure(6)
hold on
for i = 1:size(pos_p,1);line([t_trial(pos_p(i,1)) t_trial(pos_p(i,2))],[-1*iCond -1*iCond],'LineWidth',12,'Color',cmapLines(iCond,:),'Tag','main'); end
boxplot(latency_matrix(:,iCond)/1000,'orientation','horizontal','boxstyle','outline','colors',[0 0 0],'medianstyle','line','position',-1*iCond,'widths',0.6,'symbol','+');
h=findobj(gcf,'Type','line','-not','Tag','main'); % Get handles for all lines except for the main lines
set(h,'LineWidth',1); % Change line width
hold off
end

end
%%
if strcmp(main_plot,'on')
figure(5)
legend(h_p,'Exp 1: Mathôt et al. (2013)''s paradigm','Exp 2: Attend-to-location','Exp 2: Attend-to-gender','Exp 3: With diotic trials','Exp 4: Loudspeakers','Location','NorthWest');
legend boxoff
xlim([-0.5 7])
ylim([-0.2 0.7])
% xlabel('Time from cue onset (s)')
ylabel('Diff. norm. pupil diameter [z]');
set(gca,'FontSize',20,'XTickLabel',[])
set(gcf,'Color',[1 1 1])
box off

figure(6)
xlim([-0.5 7])
ylim([-6 0])
xlabel('Time from cue onset (s)')
set(gca,'FontSize',20,'YTickLabel',[])
set(gcf,'Color',[1 1 1])
box off
end