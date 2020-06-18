function [modout,errorval] = CombinedMCCrowdPopModelFun_MinProbNoNeg(EP,MP,data)
%CombinedMCCrowdPopModelFun
%
%function [modout,errorval] = CombinedMCCrowdPopModelFun(ExParam,ModParam)
%population model to predict the motion crowding data with responses at the population level
%v5.2
%
% NB this version essentially identical to Independent model but only max of WF is used
%
%J Greenwood October 2019

modout.version = 5.2; %denote version of model here to be output with files for future reference

%% model parameters

%free parameters all come from ModParam (MP) structure
%except fixed parameters:
TarBaseResp       = 0;
FlankBaseResp     = 0; 
DirTarPeakTotal   = 1;
DirTarNegTotal    = 0;
DirFlankPeakTotal = 1;
DirFlankNegTotal  = 0;
ColTarPeakTotal   = 1;
ColTarNegTotal    = 0;
ColFlankPeakTotal = 1;
ColFlankNegTotal  = 0;

DirAxis       = -180:1:180; %all possible directions (used for population response plotting)

%% generate weighting fields and specific weighting values for this configuration

WFdim  = -180:0.1:180; %values between -180 to 180 to generate the Gaussian weighting fields

%weighting field distributions (across direction)
DirWFposDist = DrawGaussian(WFdim,0,MP.DirWFPosSD,MP.DirWFPosPeak,0);
DirWFnegDist = DrawBimodalGaussian(WFdim,0-(0.5*MP.DirWFNegDelT),0+(0.5*MP.DirWFNegDelT),MP.DirWFNegSD,MP.DirWFNegPeak,MP.DirWFNegPeak,0);

%weighting field distributions (across hue)
ColWFposDist = DrawGaussian(WFdim,0,MP.ColWFPosSD,MP.ColWFPosPeak,0);
ColWFnegDist = DrawBimodalGaussian(WFdim,0-(0.5*MP.ColWFNegDelT),0+(0.5*MP.ColWFNegDelT),MP.ColWFNegSD,MP.ColWFNegPeak,MP.ColWFNegPeak,0);

%pre-generate NaN values for weights (used in non-crowded conditions)
DirWFposVal = NaN(EP.NumExpts,max(EP.NumConds));
DirWFnegVal = NaN(EP.NumExpts,max(EP.NumConds));
ColWFposVal = NaN(EP.NumExpts,max(EP.NumConds));
ColWFnegVal = NaN(EP.NumExpts,max(EP.NumConds));

%get the weights
for expt=2:EP.NumExpts %should be 4 experiments (but first is uncrowded so skip - Just the 3 strong/weak crowded combos)
    for cond = 1:EP.NumConds(expt) %should be 1 4 4 4 for each experiment in terms of flanker condition
        %direction
        [~,DirWFind]=min(abs(WFdim-EP.FlankDir(expt,cond)));
        %hue
        [~,ColWFind]=min(abs(WFdim-EP.FlankCol(expt,cond)));
        
        %find the smallest and use the same for both
        DirWFposVal(expt,cond)  = min([DirWFposDist(DirWFind) ColWFposDist(ColWFind)]);
        DirWFnegVal(expt,cond)  = min([DirWFnegDist(DirWFind) ColWFnegDist(ColWFind)]);
        
        ColWFposVal(expt,cond)  = min([DirWFposDist(DirWFind) ColWFposDist(ColWFind)]);
        ColWFnegVal(expt,cond)  = min([DirWFnegDist(DirWFind) ColWFnegDist(ColWFind)]);
    end
end

%% simulate the trials

ModelDirPropCorr = NaN(EP.NumExpts,max(EP.NumConds));
ModelColPropCorr = NaN(EP.NumExpts,max(EP.NumConds));

for expt=1:EP.NumExpts %should be 4 experiments (uncrowded + 3 strong/weak crowded combos)
    for cond = 1:EP.NumConds(expt) %should be 1 4 4 4 for each experiment in terms of flanker condition
        if expt==1% Unflanked Experiment  %NB. here we skip the 361 detectors and jump straight to a Gaussian population response
            %target population for direction
            temppos    = repmat(DrawGaussianWrap(DirAxis,EP.TarDir(expt),MP.DirExciteSD,DirTarPeakTotal,0),[MP.NumTrialsPerStim 1]) + (randn(MP.NumTrialsPerStim,numel(DirAxis)).*MP.DirNoise);
            tempneg    = repmat(DrawGaussianWrap(DirAxis,EP.TarDir(expt),MP.DirNegSD,DirTarNegTotal,0),[MP.NumTrialsPerStim 1]) + (randn(MP.NumTrialsPerStim,numel(DirAxis)).*MP.DirNoise);
            TarRespDir = temppos-tempneg + 0; %prob1 - prob2 + offsetEst;
            %target population for colour
            temppos    = repmat(DrawGaussianWrap(DirAxis,EP.TarCol(expt),MP.ColExciteSD,ColTarPeakTotal,0),[MP.NumTrialsPerStim 1]) + (randn(MP.NumTrialsPerStim,numel(DirAxis)).*MP.ColNoise);
            tempneg    = repmat(DrawGaussianWrap(DirAxis,EP.TarCol(expt),MP.ColNegSD,ColTarNegTotal,0),[MP.NumTrialsPerStim 1]) + (randn(MP.NumTrialsPerStim,numel(DirAxis)).*MP.ColNoise);
            TarRespCol = temppos-tempneg + 0; %prob1 - prob2 + offsetEst;
            
            %direction
            %repeat for number of trials and add noise to each point
            MeanTarRespDir(expt,:,cond) = mean(TarRespDir,1); %store the mean distribution to see the model workings
            %temp = MaxMin(temp,0,1); %round to 0-1 to only use positive responses
            [maxVal,maxInd]     = max(TarRespDir,[],2); %get the maximum response value in each noise-corrupted population response
            DirMax(expt,:,cond) = DirAxis(maxInd)'; %find the x-axis locations of each maximum
            %colour
            %repeat for number of trials and add noise to each point
            MeanTarRespCol(expt,:,cond) = mean(TarRespCol,1); %store the mean distribution to see the model workings
            %temp = MaxMin(temp,0,1); %round to 0-1 to only use positive responses
            [maxVal,maxInd]     = max(TarRespCol,[],2); %get the maximum response value in each noise-corrupted population response
            ColMax(expt,:,cond) = DirAxis(maxInd)'; %find the x-axis locations of each maximum
            
        else %Flanked experiment  %NB. here we skip the 361 detectors and jump straight to a Gaussian population response
            %target population for direction
            temppos    = (repmat(DrawGaussianWrap(DirAxis,EP.TarDir(expt),MP.DirExciteSD,DirTarPeakTotal,0),[MP.NumTrialsPerStim 1]) + (randn(MP.NumTrialsPerStim,numel(DirAxis)).*MP.DirNoise)).*((1-DirWFposVal(expt,cond)).*ones(MP.NumTrialsPerStim,numel(DirAxis)));
            tempneg    = (repmat(DrawGaussianWrap(DirAxis,EP.TarDir(expt),MP.DirNegSD,DirTarNegTotal,0),[MP.NumTrialsPerStim 1]) + (randn(MP.NumTrialsPerStim,numel(DirAxis)).*MP.DirNoise)).*((1-DirWFnegVal(expt,cond)).*ones(MP.NumTrialsPerStim,numel(DirAxis)));
            TarRespDir = temppos-tempneg + 0; %prob1 - prob2 + offsetEst;
            
            %direction flanker response (same general shape for each target Dir
            temppos = (repmat(DrawGaussianWrap(DirAxis,EP.FlankDir(expt,cond),MP.DirExciteSD,DirFlankPeakTotal,0),[MP.NumTrialsPerStim 1 1]) + (randn(MP.NumTrialsPerStim,numel(DirAxis)).*MP.LateNoise)).*((DirWFposVal(expt,cond)).*ones(MP.NumTrialsPerStim,numel(DirAxis))); %FlankPeakResp+FlankNegMax,0); %(x,uEst,varEst,scaleEst,offsetEst)
            tempneg = (repmat(DrawGaussianWrap(DirAxis,EP.FlankDir(expt,cond),MP.DirNegSD,DirFlankNegTotal,0),[MP.NumTrialsPerStim 1 1]) + (randn(MP.NumTrialsPerStim,numel(DirAxis)).*MP.LateNoise)).*((DirWFnegVal(expt,cond)).*ones(MP.NumTrialsPerStim,numel(DirAxis)));
            FlankRespDir = temppos-tempneg+0; %prob1 - prob2 + offsetEst;
            %combine them
            ComboRespDir = TarRespDir+FlankRespDir;
            
            %target population for colour
            temppos    = (repmat(DrawGaussianWrap(DirAxis,EP.TarCol(expt),MP.ColExciteSD,ColTarPeakTotal,0),[MP.NumTrialsPerStim 1]) + (randn(MP.NumTrialsPerStim,numel(DirAxis)).*MP.ColNoise)).*((1-ColWFposVal(expt,cond)).*ones(MP.NumTrialsPerStim,numel(DirAxis)));
            tempneg    = (repmat(DrawGaussianWrap(DirAxis,EP.TarCol(expt),MP.ColNegSD,ColTarNegTotal,0),[MP.NumTrialsPerStim 1]) + (randn(MP.NumTrialsPerStim,numel(DirAxis)).*MP.ColNoise)).*((1-ColWFnegVal(expt,cond)).*ones(MP.NumTrialsPerStim,numel(DirAxis)));
            TarRespCol = temppos-tempneg + 0; %prob1 - prob2 + offsetEst;
            
            %colour flanker response (same general shape for each target Dir
            temppos = (repmat(DrawGaussianWrap(DirAxis,EP.FlankCol(expt,cond),MP.ColExciteSD,ColFlankPeakTotal,0),[MP.NumTrialsPerStim 1 1]) + (randn(MP.NumTrialsPerStim,numel(DirAxis)).*MP.LateNoise)).*((ColWFposVal(expt,cond)).*ones(MP.NumTrialsPerStim,numel(DirAxis))); %FlankPeakResp+FlankNegMax,0); %(x,uEst,varEst,scaleEst,offsetEst)
            tempneg = (repmat(DrawGaussianWrap(DirAxis,EP.FlankCol(expt,cond),MP.ColNegSD,ColFlankNegTotal,0),[MP.NumTrialsPerStim 1 1]) + (randn(MP.NumTrialsPerStim,numel(DirAxis)).*MP.LateNoise)).*((ColWFnegVal(expt,cond)).*ones(MP.NumTrialsPerStim,numel(DirAxis)));
            FlankRespCol = temppos-tempneg+0; %prob1 - prob2 + offsetEst;
            %combine them
            ComboRespCol = TarRespCol+FlankRespCol;
            
            %store means
            MeanTarRespDir(expt,:,cond)      = mean(TarRespDir,1); %just the target response
            MeanFlankRespDir(expt,:,cond)    = mean(FlankRespDir,1); %just the flanker response
            MeanComboRespDir(expt,:,cond)    = mean(ComboRespDir,1); %store the mean distribution to see the model workings
            
            MeanTarRespCol(expt,:,cond)      = mean(TarRespCol,1); %just the target response
            MeanFlankRespCol(expt,:,cond)    = mean(FlankRespCol,1); %just the flanker response
            MeanComboRespCol(expt,:,cond)    = mean(ComboRespCol,1); %store the mean distribution to see the model workings
            
            %work out model responses/prop correct
            [maxVal,maxInd] = max(ComboRespDir,[],2); %get the maximum response value in each noise-corrupted population response
            DirMax(expt,:,cond) = DirAxis(maxInd)'; %find the x-axis locations of each maximum
            [maxVal,maxInd] = max(ComboRespCol,[],2); %get the maximum response value in each noise-corrupted population response
            ColMax(expt,:,cond) = DirAxis(maxInd)'; %find the x-axis locations of each maximum
        end
        %work out proportion correct based on trials - same process for all conditions/experiments
        ModelDirPropCorr(expt,cond) = sum((DirMax(expt,:,cond))>0)./numel(DirMax(expt,:,cond));%>0 is correct
        ModelColPropCorr(expt,cond) = sum((ColMax(expt,:,cond))>0)./numel(ColMax(expt,:,cond));%>0 is correct
    end
end

%% organise outputs

%weighting field characteristics
modout.DirWFposDist = DirWFposDist;
modout.DirWFnegDist = DirWFnegDist;
modout.ColWFposDist = ColWFposDist;
modout.ColWFnegDist = ColWFnegDist;
modout.WFdim        = WFdim;
modout.DirWFposVal  = DirWFposVal;
modout.DirWFnegVal  = DirWFnegVal;
modout.ColWFposVal  = ColWFposVal;
modout.ColWFnegVal  = ColWFnegVal;

%model population distributions and responses
modout.MeanTarRespDir = MeanTarRespDir;
modout.MeanFlankRespDir = MeanFlankRespDir;
modout.MeanComboRespDir = MeanComboRespDir;

modout.MeanTarRespCol = MeanTarRespCol;
modout.MeanFlankRespCol = MeanFlankRespCol;
modout.MeanComboRespCol = MeanComboRespCol;

%model Prop Corr
modout.ModelDirPropCorr = ModelDirPropCorr;
modout.ModelColPropCorr = ModelColPropCorr;

%% work out squared error
%note data for both TE and bias are equalised to run between -1 to 1 now,
%sim data is equalised to the same range (limited by actual data so those
%values can go larger than ±1)

datasim(:,:,1) = ModelDirPropCorr;
datasim(:,:,2) = ModelColPropCorr;

diffvals = ((data(:)-datasim(:)).^2);
errorval = sum(diffvals(~isnan(diffvals))); %remove the NaNs from the uncrowded expt

