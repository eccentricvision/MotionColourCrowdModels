%ColourCrowdPopModel
%model to predict Michael's colour crowding data with responses at the population level
%
%Uses PopModel v5.0 - same as MotionCrowdPopModel (see there for details)
%
%J Greenwood August 2018


clear all;
tic
%% experimental parameters and data

ExParam.DirAxis       = -180:1:180; %all possible directions (used for population response plotting)

ExParam.TarDir(1,:)  = [-12    -9    -6    -3     0     3     6     9    12]; %uncrowded
ExParam.TarDir(2,:)  = [-18.0 -13.5 -9.0 -4.5 0.0 4.5 9.0 13.5 18.0]; %crowded
ExParam.BaseDir       = 0;
ExParam.TarDir       = ExParam.TarDir-ExParam.BaseDir; %easier to do the directions if it's ±0 for the sign
ExParam.FlankDirDiffs = [-180  -150  -135   -45   -30   -15     0    15    30    45   135   150   180];

ExParam.NumTarVals   = numel(ExParam.TarDir(1,:));
ExParam.NumFlankVals = numel(ExParam.FlankDirDiffs);
ExParam.NumFlanks    = 2;

%mean data to simulate
MeanBias = [1.1364    3.8719    5.0929   10.2376    6.8444    7.3398    2.4732   -8.1357   -7.6515  -11.3655   -3.4262   -2.1587    1.1364]; %mean of n=6 bias values

MeanUCbias   = 1.0342;
MeanUCthresh = 2.4770;
MeanCthresh  = [2.5248    3.0371    3.7093    5.4851    4.5019    4.8636    5.1070    6.0277    5.8941    8.0284    4.0793    4.7310    2.5248];

%% model parameters

%free parameters
ModParam.ExciteSD  = 65.2458;
ModParam.NegSD     = 0.01; %NB need SD to be >0 to avoid errors in Gaussian function, but peak is 0 anyway

ModParam.WFPosPeak = 0.2358;
ModParam.WFNegPeak = 0;

ModParam.WFPosSD   = 124.6079;
ModParam.WFNegSD   = 0.01;

ModParam.DirNoise  = 0.0044;
ModParam.LateNoise = 0.0802;

%fixed parameters
ModParam.TarPeak = 1;
ModParam.TarNeg  = 0;

ModParam.NumTrialsPerStim = 1024;

%% run the model

[modout,errval] = ColourCrowdPopModelFun(ExParam,ModParam,1);
toc

fprintf('Sum of squares error value: %4.3f \n',errval);


%% plot weighting field values

figure
subplot(1,2,1);
h(1)=plot(modout.WFdim,modout.WFposDist,'r-');
hold on;
h(2)=plot(modout.WFdim,modout.WFnegDist,'b-');
xlim([-180 180]);
xtick([-180:45:180]);
title('Crowding weighting fields');
xlabel('Flanker hue difference (deg.)');
ylabel('Weight of flanker responses');
ylim([0 1]);
ytick(0:0.2:1);
axis square;
box off;
%axis square;
legend(h,{'Positive','Negative'});

subplot(1,2,2) %get a sense of how the flanker response population profile changes with direction differences
clear h; clear colmapvals;
midcolval = 0.65;
colmapvals(:,1) = [linspace(0,midcolval,240)';linspace(midcolval,1,160)'];
colmapvals(:,2) = [linspace(0,midcolval,240)';linspace(midcolval,0,160)'];
colmapvals(:,3) = [linspace(1,midcolval,240)';linspace(midcolval,0,160)'];


PlotFlankVals = -180:1:180;
clear FlankResp;
for ff=1:numel(PlotFlankVals)
    %calculate weighting field values for all flanker differences
    [~,WFind]=min(abs(modout.WFdim-PlotFlankVals(ff)));
    
    WFposVal(ff)  = modout.WFposDist(WFind(1));
    WFnegVal(ff)  = modout.WFnegDist(WFind(1));
    
    %calculate base flanker response without noise
    temppos = DrawGaussian(ExParam.DirAxis,PlotFlankVals(ff),ModParam.ExciteSD,1,0).*WFposVal(ff); 
    tempneg = DrawGaussian(ExParam.DirAxis,PlotFlankVals(ff),ModParam.NegSD,1,0).*WFnegVal(ff);
    FlankResp(ff,:) = temppos-tempneg+0; %prob1 - prob2 + offsetEst;
end
colormap(colmapvals);%(parula);%colmap);
contourf(ExParam.DirAxis,PlotFlankVals,FlankResp,9,'LineStyle','--');
xlabel('Detector direction (deg)');
ylabel('Flanker direction difference (deg)');%('Flanker response profile');
xtick(-180:45:180);
ytick(-180:45:180);
caxis([-0.75 0.5]);
title('Flanker response after weighting');
colorbar('Ticks',[-0.75:0.25:0.75]);
axis square;
box off;
hold on;
plot(ExParam.FlankDirDiffs,ExParam.FlankDirDiffs,'ko','MarkerFaceColor',[0 0 0]);
%print(gcf,'-depsc','-painters','ColourFlankerRespWithWeights.eps');%code to get around the vector export issue in matlab - prints with distinct elements

%% plot mean distributions

colvals = cool(2);

figure;%first plot uncrowded response
h=plot(ExParam.DirAxis,modout.UCmeanDist,'-');
xlim([-180 180]);
xtick(-180:90:180);
ylim([-1 1]);

hold on;
plot(xlim,[0 0],'k--');
plot([0 0],ylim,'k--');
axis square;
box off;
title('Uncrowded TarResp');
legend(h,num2str(ExParam.TarDir(1,:)'));

figure; %plot flanker responses
GMFlankResp = squeeze(mean(modout.MeanFlankResp,1)); %average over target orientations
h=plot(ExParam.DirAxis,GMFlankResp,'-');
hold on;
xlim([-180 180]);
xtick(-180:90:180);
ylim([-1 1]);
plot(xlim,[0 0],'k--');
plot([0 0],ylim,'k--');
axis square;
box off;
title('FlankerResp');
legend(h,num2str(ExParam.FlankDirDiffs'));

figure; %now plot crowded response distributions (combined)
cnt=1;
for ff=1:ExParam.NumFlankVals
    subplot(3,ceil(ExParam.NumFlankVals/3),cnt)
    h=plot(ExParam.DirAxis,modout.CmeanDist(:,:,ff),'-');
    hold on
    xlim([-180 180]);
    xtick(-180:90:180);
    ylim([-1 1]);
    plot(xlim,[0 0],'k--');
    plot([0 0],ylim,'k--');
    axis square;
    box off;
    if cnt==ExParam.NumFlankVals
        legend(h,num2str(ExParam.TarDir(2,:)'));
    end
    title(strcat(num2str(ExParam.FlankDirDiffs(ff)),'deg ComboResp'));
    cnt=cnt+1;
end

%% plot example combinations of target and flanker response

TarDirCond   = 4; %-3deg CW
FlankDirCond = 9; %30deg CCW

clear h;
%small-diff assimilation
figure;%first plot target response
h(1)=plot(ExParam.DirAxis,modout.MeanTarResp(TarDirCond,:,FlankDirCond),'r-'); %format is (target orient, detector, flanker cond)
hold on
h(2)=plot(ExParam.DirAxis,modout.MeanFlankResp(TarDirCond,:,FlankDirCond),'g-'); %format is (target orient, detector, flanker cond)
h(3)=plot(ExParam.DirAxis,modout.CmeanDist(TarDirCond,:,FlankDirCond),'b-'); %format is (target orient, detector, flanker cond)
%plot means
plot(ExParam.TarDir(2,TarDirCond),0.05,'rv');
plot(ExParam.FlankDirDiffs(FlankDirCond),0.05,'gv');
plot(ExParam.DirAxis(modout.CmeanDist(TarDirCond,:,FlankDirCond)==max(modout.CmeanDist(TarDirCond,:,FlankDirCond))),0.05,'bv');

xlim([-180 180]);
xtick([-180:45:180]);
ylim([-0.6 1]);
ytick(-1:0.2:1);
plot(xlim,[0 0],'k--');
plot([0 0],ylim,'k--');
axis square;
box off;
title('Example distribution comparison (NB weights already applied)');
legend(h,{'Target','Flankers','Combined'});

%large-diff (still assimilation) example
TarDirCond   = 4; %-3deg CW
FlankDirCond = 11; %120deg CCW

clear h;
%small-diff assimilation
figure;%first plot target response
h(1)=plot(ExParam.DirAxis,modout.MeanTarResp(TarDirCond,:,FlankDirCond),'r-'); %format is (target orient, detector, flanker cond)
hold on
h(2)=plot(ExParam.DirAxis,modout.MeanFlankResp(TarDirCond,:,FlankDirCond),'g-'); %format is (target orient, detector, flanker cond)
h(3)=plot(ExParam.DirAxis,modout.CmeanDist(TarDirCond,:,FlankDirCond),'b-'); %format is (target orient, detector, flanker cond)
%plot means
plot(ExParam.TarDir(2,TarDirCond),0.05,'rv');
plot(ExParam.FlankDirDiffs(FlankDirCond),0.05,'gv');
plot(ExParam.DirAxis(modout.CmeanDist(TarDirCond,:,FlankDirCond)==max(modout.CmeanDist(TarDirCond,:,FlankDirCond))),0.05,'bv');

xlim([-180 180]);
xtick([-180:45:180]);
ylim([-0.6 1]);
ytick(-1:0.2:1);
plot(xlim,[0 0],'k--');
plot([0 0],ylim,'k--');
axis square;
box off;
title('Example distribution comparison (NB weights already applied)');
legend(h,{'Target','Flankers','Combined'});


%% plot psychometric functions

markers = {'o','^','s','*','v','+','d','x','o','^','s','*','v','+','d','x'};
marface = lines(20);

figure
%first plot uncrowded response
subplot(3,round((ExParam.NumFlankVals+1)/3),1)

Xaxis = ExParam.TarDir(1,:)+ExParam.BaseDir; %directions tested in target
Xfine = min(Xaxis):0.001:max(Xaxis);
h=plot(Xaxis,modout.UncrowdProbCCW,'o','Color',marface(1,:),'MarkerFaceColor',marface(1,:)); %plot data
hold on
[UCu,UCv,UCkp,UCcuts,UCfb] = FitCumuGaussian(Xaxis,modout.UncrowdProbCCW,ModParam.NumTrialsPerStim.*ones(1,ExParam.NumTarVals),0,0,[1 1 1],[],0,1);
curvefit=DrawCumuGaussian(Xfine,UCu,UCv,UCkp,0,UCfb);
plot(Xfine,curvefit,'-','Color',marface(1,:)); %plot curvefit
xlim([min(Xaxis)-3 max(Xaxis)+3])
xtick(Xaxis);
ylim([0 1]);
ytick([0 0.25 0.5 0.75 1]);
xlabel('Target hue (deg)');
ylabel('Prop. CCW (pink)');
title('Uncrowded responses');
axis square;

%now plot the crowded ones
Xaxis = ExParam.TarDir(2,:)+ExParam.BaseDir; %directions tested in target;
Xfine = min(Xaxis):0.001:max(Xaxis);

cnt=2;
for ff=1:ExParam.NumFlankVals
    subplot(3,round((ExParam.NumFlankVals+1)/3),cnt)
    
    h(ff)=plot(Xaxis,modout.CrowdProbCCW(ff,:),'o','Color',marface(ff+1,:),'MarkerFaceColor',marface(ff+1,:)); %plot data
    hold on
    [Cu(ff) Cv(ff) Ckp(ff) Ccuts(ff,:) Cfb(ff)] = FitCumuGaussian(Xaxis,modout.CrowdProbCCW(ff,:),ModParam.NumTrialsPerStim.*ones(1,ExParam.NumTarVals),0,0,[1 1 1],[],0,1);
    curvefit=DrawCumuGaussian(Xfine,Cu(ff),Cv(ff),Ckp(ff),0,Cfb(ff));
    plot(Xfine,curvefit,'-','Color',marface(ff+1,:)); %plot curvefit
    
    xlim([min(Xaxis)-3 max(Xaxis)+3])
    xtick(Xaxis);
    ylim([0 1]);
    ytick([0 0.25 0.5 0.75 1]);
    xlabel('Target hue (deg)');
    ylabel('Prop. CCW (pink)');
    title(strcat('Crowded with flankers at',num2str(ExParam.FlankDirDiffs(ff)),'deg'));
    axis square;
    cnt=cnt+1;
end

%% plot the midpoints, thresholds, and threshold elevation values

figure
clear h;

%plot midpoints
subplot(1,2,1);
Xaxis = ExParam.FlankDirDiffs;
h(1)  = plot(Xaxis,modout.Cmidpoint,'-','Color',marface(1,:),'MarkerFaceColor',marface(1,:)); %plot simulated midpoint data
hold on;
h(2)  = plot(Xaxis,MeanBias,'o','Color',marface(3,:),'MarkerFaceColor',marface(3,:));      %plot actual midpoint data
plot(xlim,[modout.UCmidpoint modout.UCmidpoint],'-','Color',marface(1,:));%,'MarkerFaceColor',[1 1 1]); %plot simulated UC data
plot(xlim,[MeanUCbias MeanUCbias],':','Color',marface(3,:));
plot(xlim,[90 90],'k--');
xtick(-180:30:180);
xlabel('Flanker hue difference (deg.)')
ylabel('Midpoint (deg.)');
axis square
title('Midpoints');
MaxDiff = max(abs(modout.Cmidpoint-90));
ylim([-20 20]);
ytick(-20:5:20);
xlim([min(Xaxis)-10 max(Xaxis)+10]);
text(-160,ExParam.BaseDir+MaxDiff,'Assimilation');
text(-160,ExParam.BaseDir-(MaxDiff),'Repulsion');
text(100,ExParam.BaseDir+MaxDiff,'Repulsion');
text(100,ExParam.BaseDir-(MaxDiff),'Assimilation');

%plot thresholds

subplot(1,2,2);
h(1)  = plot(Xaxis,modout.Cthreshold,'-','Color',marface(1,:),'MarkerFaceColor',marface(1,:)); %plot threshold data
hold on;
h(2)  = plot(Xaxis,MeanCthresh,'o','Color',marface(3,:),'MarkerFaceColor',marface(3,:));      %plot actual TE data
plot(xlim,[modout.UCthreshold modout.UCthreshold],'-','Color',marface(1,:));%,plot simulated UC threshold
plot(xlim,[MeanUCthresh MeanUCthresh],':','Color',marface(3,:)); %measured UC threshold
plot(xlim,[0 0],'k--');
xtick(-180:30:180);
xlabel('Flanker hue difference (deg.)')
ylabel('Threshold (deg)');
axis square
title('Threshold values');
legend(h,{'Model','Data'});
ylim([0 max([modout.Cthreshold(:)' MeanCthresh(:)']+1)]);
xlim([min(Xaxis)-10 max(Xaxis)+10]);


