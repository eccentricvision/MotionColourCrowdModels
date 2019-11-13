%CombinedMCCrowdPopModel_NoNeg
%model to predict Michael's motion-colour conjoint crowding data
%v5.2 for both colour and motion
%population processing version based on MotionCrowdPopModel.m and ColourCrowdPopModel.m, both v5.2
%
% NB this version essentially identical to Independent model but only max of WF is used
%
%J Greenwood September 2019

%% data values

clear all;
tic

%direction - mean data
DataPoolPropMean(1,:,1) = [0.8771    NaN NaN NaN];
DataPoolPropMean(2,:,1) = [0.9118    0.3563    0.8931    0.3847];
DataPoolPropMean(3,:,1) = [0.9111    0.4271    0.9208    0.4292];
DataPoolPropMean(4,:,1) = [0.7451    0.9062    0.7493    0.9076];
%colour - mean data
DataPoolPropMean(1,:,2) = [0.9396    NaN NaN NaN];
DataPoolPropMean(2,:,2) = [0.9562    0.9437    0.2465    0.2569];
DataPoolPropMean(3,:,2) = [0.8694    0.8604    0.8722    0.8708];
DataPoolPropMean(4,:,2) = [0.9507    0.9458    0.3528    0.3382];

%% experiment parameters

ExParam.NumFlanks        = 2;
ExParam.NumExpts         = 4; %uncrowded + strong-strong/weak-strong/strong-weak crowding experiments for motion&colour
ExParam.NumConds         = [1 4 4 4]; %number of conditions per expt - bothsame/ordiff/posdiff/bothdiff

ExParam.DirVals = [8 15 165]; %in order: target direction, small flanker dir, large flanker dir
ExParam.ColVals = [5 30 150]; %the same for colour

ExParam.TarDir        = [1 1 1 1].*ExParam.DirVals(1); %same for each experiment (mean of all observers' individual values)
ExParam.FlankDir(1,:) = [NaN NaN NaN NaN];%uncrowded
ExParam.FlankDir(2,:) = [1 -1 1 -1].*ExParam.DirVals(2); %strong-strong: both same / colour differs / motion differs / both differ
ExParam.FlankDir(3,:) = [1 -1 1 -1].*ExParam.DirVals(2); %strong-weak: both same / colour differs / motion differs / both differ
ExParam.FlankDir(4,:) = [1 -1 1 -1].*ExParam.DirVals(3); %weak-strong

ExParam.TarCol        = [1 1 1 1].*ExParam.ColVals(1); %mean of all observers values
ExParam.FlankCol(1,:) = [NaN NaN NaN NaN];%uncrowded
ExParam.FlankCol(2,:) = [1 1 -1 -1].*ExParam.ColVals(2);%strong-strong: both same / colour differs / motion differs / both differ
ExParam.FlankCol(3,:) = [1 1 -1 -1].*ExParam.ColVals(3);%strong-weak: both same / colour differs / motion differs / both differ
ExParam.FlankCol(4,:) = [1 1 -1 -1].*ExParam.ColVals(2);%weak-strong: both same / colour differs / motion differs / both differ

%% model parameters

%fixed parameters
ModParam.NumTrialsPerStim = 1024;

ModParam.DirExciteSD  = 65.3691;
ModParam.DirNegSD     = 0.01;

ModParam.DirWFPosPeak = 0.5069;
ModParam.DirWFNegPeak = 0;
ModParam.DirWFPosSD   = 70.8827;
ModParam.DirWFNegSD   = 0.01;
ModParam.DirWFNegDelT = 0;

ModParam.ColExciteSD = 65.2458;
ModParam.ColNegSD    = 0.01;

ModParam.ColWFPosPeak = 0.2358;
ModParam.ColWFNegPeak = 0;
ModParam.ColWFPosSD   = 124.6079;
ModParam.ColWFNegSD   = 0.01;
ModParam.ColWFNegDelT = 0;

%free parameters
ModParam.DirNoise  = 0.0192;
ModParam.ColNoise  = 0.0090;
ModParam.LateNoise = 0.1879;

NumFreeParam = 3;

%% run the model

[modout,errval] = CombinedMCCrowdPopModelFun_MaxProbNoNeg(ExParam,ModParam,DataPoolPropMean);
toc

[AICval,AICc] = ComputeAIC(errval,26,NumFreeParam);
fprintf('Sum of squares error value: %4.3f \n',errval);
fprintf('AIC: %4.3f \n',AICval);

%% plot the simulation + data

markertype = {'o','^','s','*','v','+','d','x','o','^','s','*','v','+','d','x'};
marface = [1 0 0; 0 0.8 0; 0 0 1; 0.8 0.1 0.8; 0 0 0; 0.5 0.5 0; 0 0.5 0.5; 0 0.25 0.5; 0.25 0.25 0.75; 0.25 0.25 0.75; 0 0.25 0.5; 0 0.5 0.5; 0.5 0.5 0; 0.5 0.5 0.5; 0 0 0;1 0 0; 0 0.8 0; 0 0 1; 0.8 0 0.8];

CondLabels = {'Both Match','Motion Differs','Colour Differs','Both Differ'};
ExptLabels = {'Uncrowded','Small Motion-Small Colour','Small Motion-Large Colour','Large Motion-Small Colour'};

figure

for expt=1:4
    subplot(2,2,expt)
    clear h;
    if expt==1 %uncrowded
        xlim([0 1]);
        xtick([0 0.5 1]);
        xlabel('Prop Corr Motion');
        ylim([0 1]);
        ytick([0 0.5 1]);
        ylabel('Prop Corr Colour');
        
        h(1) = plot(modout.ModelDirPropCorr(expt,:),modout.ModelColPropCorr(expt,:),'Marker',markertype{3},'LineStyle','none','color',marface(1,:),'MarkerFaceColor',[1 1 1],'MarkerSize',10);
        hold on;
        plot(DataPoolPropMean(expt,:,1),DataPoolPropMean(expt,:,2),'Marker',markertype{1},'LineStyle','none','color',marface(1,:),'MarkerFaceColor',marface(1,:),'MarkerSize',10);
        plot([0 1],[0.5 0.5],'k--'); plot([0.5 0.5],[0 1],'k--');
        text(0.25,0.25,'Both Error');
        text(0.25,0.75,'Motion Error');
        text(0.75,0.25,'Colour Error');
        text(0.75,0.75,'Both Correct');
        
        xlim([0 1]);
        xtick([0 0.5 1]);
        xlabel('Prop Corr Motion');
        ylim([0 1]);
        ytick([0 0.5 1]);
        ylabel('Prop Corr Colour');
        
        legend(h,'Uncrowded');
        title(ExptLabels{expt})
        box on;
        axis square;
        
    else %crowded conditions
        clear h;
        
        xlim([0 1]);
        xtick([0 0.5 1]);
        xlabel('Prop Corr Motion');
        ylim([0 1]);
        ytick([0 0.5 1]);
        ylabel('Prop Corr Colour');
        
        for cond=1:4
            h(cond) = plot(modout.ModelDirPropCorr(expt,cond),modout.ModelColPropCorr(expt,cond),'Marker',markertype{3},'LineStyle','none','color',marface(cond,:),'MarkerFaceColor',[1 1 1],'MarkerSize',10);
            hold on;
            plot(DataPoolPropMean(expt,cond,1),DataPoolPropMean(expt,cond,2),'Marker',markertype{1},'LineStyle','none','color',marface(cond,:),'MarkerFaceColor',marface(cond,:),'MarkerSize',10);
        end
        
        plot([0 1],[0.5 0.5],'k--'); plot([0.5 0.5],[0 1],'k--');
        text(0.25,0.25,'Both Error');
        text(0.25,0.75,'Motion Error');
        text(0.75,0.25,'Colour Error');
        text(0.75,0.75,'Both Correct');
        
        xlim([0 1]);
        xtick([0 0.5 1]);
        xlabel('Prop Corr Motion');
        ylim([0 1]);
        ytick([0 0.5 1]);
        ylabel('Prop Corr Colour');
        
        legend(h,CondLabels);
        title(ExptLabels{expt})
        box on;
        axis square;
    end
end

%% plot weighting field values

figure
clear h;
subplot(2,2,1);
h(1)=plot(modout.WFdim,modout.DirWFposDist,'r-');
hold on;
h(2)=plot(modout.WFdim,modout.DirWFnegDist,'b-');
title('Direction weighting fields');
xlabel('Flanker direction difference (deg.)');
ylabel('Weight of flanker responses');
ylim([0 1]);
axis square;
box off;
%axis square;
legend(h,{'Positive','Negative'});

subplot(2,2,2);
clear h;
h(1)=plot(modout.WFdim,modout.ColWFposDist,'r-');
hold on;
h(2)=plot(modout.WFdim,modout.ColWFnegDist,'b-');
title('Colour weighting fields');
xlabel('Flanker hue difference (deg.)');
ylabel('Weight of flanker responses');
ylim([0 1]);
axis square;
box off;
%axis square;
legend(h,{'Positive','Negative'});

subplot(2,2,3) %get a sense of how the flanker response population profile changes with direction differences
clear h;
% colmap(:,1) = [ones(128,1)./2;linspace(0.5,1,127)'];
% colmap(:,2) = [linspace(1,0.5,127)';ones(128,1)./2];
% colmap(:,3) = ones(255,1).*0.5;
ExParam.DirAxis = -180:1:180;
PlotFlankVals   = -180:1:180;%ExParam.FlankDirDiffs;%[-90 -30 0 30 90];%[0 30 90 150 180]; %just a taste

clear FlankResp;
for ff=1:numel(PlotFlankVals)
    %calculate weighting field values for all flanker differences
    [~,WFind]=min(abs(modout.WFdim-PlotFlankVals(ff)));
    
    WFposVal(ff)  = modout.DirWFposDist(WFind(1));
    WFnegVal(ff)  = modout.DirWFnegDist(WFind(1));
    
    %calculate base flanker response without noise
    temppos = DrawGaussian(ExParam.DirAxis,PlotFlankVals(ff),ModParam.DirExciteSD,1,0).*WFposVal(ff); %FlankPeakResp+FlankNegMax,0); %(x,uEst,varEst,scaleEst,offsetEst)
    tempneg = DrawGaussian(ExParam.DirAxis,PlotFlankVals(ff),ModParam.DirNegSD,1,0).*WFnegVal(ff);
    FlankResp(ff,:) = temppos-tempneg+0; %prob1 - prob2 + offsetEst;
end
colormap(parula);%colmap);
contourf(ExParam.DirAxis,PlotFlankVals,FlankResp,10);
xlabel('Detector direction (deg)');
ylabel('Flanker direction difference (deg)');%('Flanker response profile');
xtick(-180:90:180);
ytick(-180:90:180);
caxis([-1 1]);
title('Flanker response after weighting');
colorbar('Ticks',[-1:0.25:1]);
axis square;

subplot(2,2,4)
clear FlankResp;
for ff=1:numel(PlotFlankVals)
    %calculate weighting field values for all flanker differences
    [~,WFind]=min(abs(modout.WFdim-PlotFlankVals(ff)));
    
    WFposVal(ff)  = modout.ColWFposDist(WFind(1));
    WFnegVal(ff)  = modout.ColWFnegDist(WFind(1));
    
    %calculate base flanker response without noise
    temppos = DrawGaussian(ExParam.DirAxis,PlotFlankVals(ff),ModParam.ColExciteSD,1,0).*WFposVal(ff); %FlankPeakResp+FlankNegMax,0); %(x,uEst,varEst,scaleEst,offsetEst)
    tempneg = DrawGaussian(ExParam.DirAxis,PlotFlankVals(ff),ModParam.ColNegSD,1,0).*WFnegVal(ff);
    FlankResp(ff,:) = temppos-tempneg+0; %prob1 - prob2 + offsetEst;
end
colormap(parula);%colmap);
contourf(ExParam.DirAxis,PlotFlankVals,FlankResp,10);
xlabel('Detector hue (deg)');
ylabel('Flanker hue difference (deg)');%('Flanker response profile');
xtick(-180:90:180);
ytick(-180:90:180);
caxis([-1 1]);
title('Flanker response after weighting');
colorbar('Ticks',[-1:0.25:1]);
axis square;

%% plot mean distributions
% 
%CondLabels = {'Both Match','Motion Differs','Colour Differs','Both Differ'};
%ExptLabels = {'Uncrowded','Small Motion-Small Colour','Small Motion-Large Colour','Large Motion-Small Colour'};

%plot for direction
colvals = cool(3);
cnt=1;
figure;%first plot uncrowded response
clear h;
for expt=1:4
    for cond=1:4
        subplot(4,4,cnt);
        h(1) = plot(ExParam.DirAxis,modout.MeanTarRespDir(expt,:,1),'-','Color',colvals(1,:)); %target response
        xlim([-180 180]);
        xtick(-180:90:180);
        ylim([-1.5 1.5]);
        ytick([-1 0 1]);
        hold on;
        plot([0 0],ylim,'k--'); %draw decision boundary
        if expt>1
            h(2) = plot(ExParam.DirAxis,modout.MeanFlankRespDir(expt,:,cond),'-','Color',colvals(2,:)); %target response
            h(3) = plot(ExParam.DirAxis,modout.MeanComboRespDir(expt,:,cond),'-','Color',colvals(3,:)); %target response
            maxloc = find(modout.MeanComboRespDir(expt,:,cond)==max(All(modout.MeanComboRespDir(expt,:,cond))));
            plot(ExParam.DirAxis(maxloc),max(All(modout.MeanComboRespDir(expt,:,cond))),'v','Color',colvals(3,:)); %plot a triangle for the peak of the combo distribution
            if ~mod(cnt,4)
            legend(h,{'Target','Flankers','Combined'});
            end
        end
        title(strcat(ExptLabels{expt},'-',CondLabels{cond},'- Direction'));
        xlabel('Direction (deg.)');
        ylabel('Response');
        cnt=cnt+1;
    end
end

%plot for colour
figure;%first plot uncrowded response
cnt=1;
clear h;
for expt=1:4
    for cond=1:4
        subplot(4,4,cnt);
        h(1) = plot(ExParam.DirAxis,modout.MeanTarRespCol(expt,:,1),'-','Color',colvals(1,:)); %target response
        xlim([-180 180]);
        xtick(-180:90:180);
        ylim([-1.5 1.5]);
        ytick([-1 0 1]);
        hold on;
        plot([0 0],ylim,'k--'); %draw decision boundary
        if expt>1
            h(2) = plot(ExParam.DirAxis,modout.MeanFlankRespCol(expt,:,cond),'-','Color',colvals(2,:)); %target response
            h(3) = plot(ExParam.DirAxis,modout.MeanComboRespCol(expt,:,cond),'-','Color',colvals(3,:)); %target response
                         maxloc = find(modout.MeanComboRespCol(expt,:,cond)==max(All(modout.MeanComboRespCol(expt,:,cond))));
            plot(ExParam.DirAxis(maxloc),max(All(modout.MeanComboRespCol(expt,:,cond))),'v','Color',colvals(3,:)); %plot a triangle for the peak of the combo distribution
            if ~mod(cnt,4)
            legend(h,{'Target','Flankers','Combined'});
            end
        end
        title(strcat(ExptLabels{expt},'-',CondLabels{cond},'- Colour'));
        xlabel('Hue (deg.)');
        ylabel('Response');
        cnt=cnt+1;
    end
end

