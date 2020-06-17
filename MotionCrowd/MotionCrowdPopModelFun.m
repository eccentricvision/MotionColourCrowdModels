function [modout,errorval] = MotionCrowdPopModelFun(EP,MP,CalcErrVal)
%MotionCrowdPopModelFun
%
%function [modout] = MotionCrowdPopModelFun(ExParam,ModParam,CalcErrVal)
%model to predict the motion crowding data with responses at the population level
%
% v5.2
% new in 5.2 - turning inhibition in the target back on and posSD+negSD of WF are collapsed to be the same
%
%J Greenwood September 2019

modout.version = 5.2; %denote version of model here to be output with files for future reference

%% model parameters

if nargin<3
    CalcErrVal = 1; %if not specified then calculate error relative to data
end

%free parameters all come from ModParam (MP) structure
%except fixed parameters:
TarBaseResp    = 0;
FlankBaseResp  = 0;

DirAxis       = -180:1:180; %all possible directions (used for population response plotting)

%% generate weighting fields and specific weighting values for this configuration

WFdim  = min(EP.FlankDirDiffs(:)):0.1:max(EP.FlankDirDiffs(:)); %values between -180 to 180 to generate the Gaussian weighting fields

%weighting field distributions (across direction) - NB same peak height for both +ve and -ve components, different bandwidths
WFposDist = DrawGaussian(WFdim,0,MP.WFPosSD,MP.WFPosPeak,0); %values from 0 - WFPeak %MP.WFPosSD,1,0);
WFnegDist = DrawBimodalGaussian(WFdim,0-(0.5*MP.WFNegDelT),0+(0.5*MP.WFNegDelT),MP.WFNegSD,MP.WFNegPeak,MP.WFNegPeak,0);

for ff=1:EP.NumFlankVals
    [~,WFind]=min(abs(WFdim-EP.FlankDirDiffs(ff)));
    
    WFposVal(ff)  = WFposDist(WFind);
    WFnegVal(ff)  = WFnegDist(WFind);
end

%% simulate the trials

% Unflanked Experiment
for tar=1:EP.NumTarVals %NB. here we skip the 361 detectors and jump straight to a Gaussian population response
    %repeat for number of trials and add noise to each point
    temppos = repmat(DrawGaussianWrap(DirAxis,EP.TarDir(1,tar),MP.ExciteSD,MP.TarPeak,0),[MP.NumTrialsPerStim 1]) + (randn(MP.NumTrialsPerStim,numel(DirAxis)).*MP.DirNoise);
    tempneg = repmat(DrawGaussianWrap(DirAxis,EP.TarDir(1,tar),MP.NegSD,MP.TarNeg,0),[MP.NumTrialsPerStim 1]) + (randn(MP.NumTrialsPerStim,numel(DirAxis)).*MP.DirNoise);
    temp    = temppos-tempneg + 0; %prob1 - prob2 + offsetEst;
    
    UCmeanDist(tar,:) = mean(temp,1); %store the mean distribution to see the model workings
    
    [maxVal,maxInd] = max(temp,[],2); %get the maximum response value in each noise-corrupted population response
    UCdirmax(tar,:) = DirAxis(maxInd)'; %find the x-axis locations of each maximum
    
    %work out proportion counterclockwise based on trials
    UncrowdProbCCW(tar) = sum((UCdirmax(tar,:))>0)./numel(UCdirmax(tar,:));%sum(temp(DirAxis>0))./sum(temp);%(sum(UncrowdModelResp,2)./NumTrialsPerStim)'; %makes final data in terms of proportion upwards
end

%Flanked experiment
for ff=1:EP.NumFlankVals
    
    DirDiffs(ff,:) = (EP.TarDir(2,:)-EP.FlankDirDiffs(ff)); %difference between target and the flanker Dir in this condition (assume identical flankers here)
    
    for tar=1:EP.NumTarVals %NB. here we skip the 361 detectors and jump straight to a Gaussian population response
        %target and flanker responses - now with target weights (1-WF) as well as flanker weights (WF)
        temppos = (repmat(DrawGaussianWrap(DirAxis,EP.TarDir(2,tar),MP.ExciteSD,MP.TarPeak,0),[MP.NumTrialsPerStim 1]) + (randn(MP.NumTrialsPerStim,numel(DirAxis)).*MP.DirNoise)).*((1-WFposVal(ff)).*ones(MP.NumTrialsPerStim,numel(DirAxis)));
        tempneg = (repmat(DrawGaussianWrap(DirAxis,EP.TarDir(2,tar),MP.NegSD,MP.TarNeg,0),[MP.NumTrialsPerStim 1]) + (randn(MP.NumTrialsPerStim,numel(DirAxis)).*MP.DirNoise)).*((1-WFnegVal(ff)).*ones(MP.NumTrialsPerStim,numel(DirAxis)));
        TarResp    = temppos-tempneg + 0; %prob1 - prob2 + offsetEst;
        
        %generate flanker response (same general shape for each target Dir - NB. make flanker distributions with peak of 1 and then multiply by flanker weights to determine magnitude
        temppos = (repmat(DrawGaussianWrap(DirAxis,EP.FlankDirDiffs(ff),MP.ExciteSD,1,0),[MP.NumTrialsPerStim 1 1]) + (randn(MP.NumTrialsPerStim,numel(DirAxis)).*MP.LateNoise)).*((WFposVal(ff)).*ones(MP.NumTrialsPerStim,numel(DirAxis))); 
        tempneg = (repmat(DrawGaussianWrap(DirAxis,EP.FlankDirDiffs(ff),MP.NegSD,1,0),[MP.NumTrialsPerStim 1 1]) + (randn(MP.NumTrialsPerStim,numel(DirAxis)).*MP.LateNoise)).*((WFnegVal(ff)).*ones(MP.NumTrialsPerStim,numel(DirAxis)));
        FlankResp = temppos-tempneg+0; %prob1 - prob2 + offsetEst;
        
        %combine them
        ComboResp = TarResp+FlankResp;
        %ComboResp = MaxMin(ComboResp,-1,1); %round to 0-1 now to ensure positive responses;
        
        %store means
        CmeanDist(tar,:,ff) = mean(ComboResp,1); %store the mean distribution to see the model workings
        MeanTarResp(tar,:,ff)   = mean(TarResp,1); %just the target response
        MeanFlankResp(tar,:,ff) = mean(FlankResp,1); %just the flanker response
        
        [maxVal,maxInd] = max(ComboResp,[],2); %get the maximum response value in each noise-corrupted population response
        Cdirmax(tar,:,ff) = DirAxis(maxInd)'; %find the x-axis locations of each maximum
        
        %work out proportion counterclockwise based on trials
        CrowdProbCCW(ff,tar) = sum((Cdirmax(tar,:,ff))>0)./numel(Cdirmax(tar,:,ff)); %makes final data in terms of proportion upwards
    end
end


%% organise outputs

%weighting field characteristics
modout.WFposDist = WFposDist;
modout.WFnegDist = WFnegDist;
modout.WFdim     = WFdim;
modout.WFposVal  = WFposVal;
modout.WFnegVal  = WFnegVal;

%model population distributions and responses
modout.UCmeanDist     = UCmeanDist;
modout.CmeanDist      = CmeanDist;
modout.UncrowdProbCCW = UncrowdProbCCW;
modout.MeanTarResp    = MeanTarResp;
modout.MeanFlankResp  = MeanFlankResp;
modout.CrowdProbCCW   = CrowdProbCCW;

%% fit psychometric functions and get threshold/midpoint values

%plot crowded responses
Cutlevels = [0.25 0.5 0.75]; %take three cut points to get midpoint and threshold

Xaxis = EP.TarDir(1,:)+EP.BaseDir; %directions tested in target
Xfine = min(Xaxis):0.001:max(Xaxis);

[modout.UCu modout.UCv modout.UCkp modout.UCcuts modout.UCfb] = FitCumuGaussian(Xaxis,modout.UncrowdProbCCW,MP.NumTrialsPerStim.*ones(1,EP.NumTarVals),0,0,[1 1 1],[],Cutlevels,1);

modout.UCmidpoint = modout.UCcuts(2);
modout.UCthreshold = abs(modout.UCcuts(3)-modout.UCcuts(2));

Xaxis = EP.TarDir(2,:)+EP.BaseDir; %directions tested in target;

for ff=1:EP.NumFlankVals
    [modout.Cu(ff) modout.Cv(ff) modout.Ckp(ff) modout.Ccuts(ff,:) modout.Cfb(ff)] = FitCumuGaussian(Xaxis,modout.CrowdProbCCW(ff,:),MP.NumTrialsPerStim.*ones(1,EP.NumTarVals),0,0,[1 1 1],[],Cutlevels,1);
    
    modout.Cmidpoint(ff) = modout.Ccuts(ff,2);
    modout.Cthreshold(ff) = abs(modout.Ccuts(ff,3)-modout.Ccuts(ff,2));
end

modout.ThreshElev = modout.Cthreshold./modout.UCthreshold;

%% work out squared error

if CalcErrVal
    %using TE with UC threshold/bias
    data(1,:) = ([90.1699 89.5355   83.3052   79.4524   73.5981   76.7235   90.4823  104.9497  105.5996   90.7504   79.2415   81.1134   95.6976  108.5280  108.5469 101.6390   96.3529   89.5355]-EP.BaseDir)./18.5469;
    data(2,:)   = ([4.4075 1.6464    1.7543    2.0779    2.8139    2.7016    2.6777    3.0878    2.9969    3.7984    3.0996    2.7666    2.6332    2.5842    2.5812  2.1407    1.9172    1.6464])./3.7984;
    
    %using TE with UC threshold/bias
    datasim(1,:) = ([modout.UCmidpoint modout.Cmidpoint]-EP.BaseDir)./18.5469; %divide by max difference in actual data to put data in the same range ±1
    datasim(2,:) = [modout.UCthreshold modout.ThreshElev]./3.7984; %divide by max threshold elevation from real data to put in range 0-1
    
    errorval=sum((datasim(:)-data(:)).^2);
else
    errorval = NaN;
end

