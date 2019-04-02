%% Function to analyze instances
%
%    v = newView;
%    v = viewSet(v, 'curGroup', 'Concatenation'); v = viewSet(v, 'curScan', 1);
%    v = loadAnalysis(v, 'erAnal/erAnal.mat');
%
%    analyzeTexInstances(v, rois);
%
function rois = analyzeTexInstances(v, rois, varargin)

getArgs(varargin, {'doZscore=0', 'useMean=0', 'texIdx', [1 2 3 4]}, 'verbose=1');

% Parse args
if ieNotDefined('v')
  v = newView;
  v = viewSet(v, 'curGroup', 'Concatenation');
  v = viewSet(v, 'curScan', 1);
  v = loadAnalysis(v, 'erAnal/erAnal.mat');
end

if ieNotDefined('rois')
  roiNames = {'V1', 'V2', 'V3', 'V4', 'LO1', 'LO2'}; rois = {};
  for i = 1:length(roiNames)
    rois{i} = loadROITSeries(v, roiNames{i});
    rois{i}.concatInfo = viewGet(v, 'concatInfo');
    rois{i}.framePeriod = viewGet(v, 'framePeriod');
    rois{i}.nFrames = viewGet(v, 'nFrames');
  end
end

%% Get condition names and r2
d = viewGet(v, 'd');
overlays = viewGet(v, 'overlays');
r2 = overlays(1).data{1}; % Get r2 map

stimNames = d.stimNames;
[stimValues, condNames] = parseConditionName(stimNames);
texFams = unique(stimValues(1,:));
%% Get the instances
rois = getSortIndex(v, rois, r2);
if useMean
  rois = getInstances(v, rois, d.stimvol, 'startLag=2','blockLen=15', 'n=200', 'r2cutoff=.1');
else
  rois = getInstances(v,rois, d.stimvol, 'type=glm','canonicalType=allfit2', 'n=200', 'r2cutoff=.1');
end
%keyboard
%%
corr_mtx = ER_getCorrMtx(rois, stimNames, stimValues, condNames, 'doZscore', doZscore);
%%
corrs = ER_getCorrStruct(rois, stimValues, condNames, corr_mtx, 1);

stimfile = viewGet(v, 'stimfile');
if iscell(stimfile)
  stimfile = stimfile{1};
end
imNames = stimfile.stimulus.imNames;

%
roiNames = cellfun(@(x) x.name, rois, 'un', 0);
plotCorrelations(corrs, roiNames, texFams, imNames);
%%
%figure; for i = 1:6; subplot(3,2,i); cm = squeeze(corr_mtx(i,:,:)); imagesc(cm); caxis([-.5 1]); colorbar; end
%keyboard

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to get full correlation matrix
% - can call either with a r2cutoff or a number of voxels.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function corr_mtx = ER_getCorrMtx(rois, stimNames, stimValues, condNames, varargin)
% Options you can set:
%    - r2cutoff
%    - nVoxels
%    - subtracted
%    - useAvgTrace
%    - texIdx   % whichsamples are textures (textures=1-4)
getArgs(varargin, {'texIdx=[1 2 3 4]', 'r2cutoff=.2', 'doZscore=0', 'subtracted=1'}, 'verbose=1');

% Extract conditions and values
isTexture = arrayfun(@(x) any(x==texIdx), stimValues(2,:));

texFams = unique(stimValues(1,:));
texSmps = find(isTexture(stimValues(2,:)));
cond_pairs = combvecs(texSmps, texSmps);
corr_mtx = zeros(length(rois), length(stimNames), length(stimNames));
corr_cond = zeros(length(rois), 4);

nSmps = floor(size(rois{1}.instance.instances{1}, 1)/2);

for ri = 1:length(rois)
  
  trialAmps = cellfun(@(x) x.', rois{ri}.instance.instances, 'un', 0);
  
  if (doZscore)
    allTrials = cat(2, trialAmps{isTexture});
    meanAmps = mean(allTrials,2); % take mean of each voxel over all presentations
    stdAmps = std(allTrials, [], 2); % std of each voxel over all presentations

    % Select whether to z-score, or just subtract the mean, or just divide by standard deviation.
    if doZscore==1
      zscoreFunc = @(x) (x-repmat(meanAmps, 1, size(x,2)))./repmat(stdAmps, 1, size(x,2));
    elseif doZscore==2 % if doZscore=2, just substract the mean
      zscoreFunc = @(x) (x-repmat(meanAmps, 1, size(x,2)));
    elseif doZscore==3
      zscoreFunc = @(x) (x./repmat(stdAmps, 1, size(x,2)));
    end
    %trialAmps = cellfun(@(x) (x - repmat(meanAmps, 1, size(x,2)))./repmat(stdAmps, 1, size(x,2)), trialAmps, 'un', 0);
    %trialAmps = cellfun(@(x) (x - repmat(meanAmps, 1, size(x,2))), trialAmps, 'un', 0);

    trialAmps = cellfun(zscoreFunc, trialAmps, 'un', 0);

  end
  
  disppercent(-inf, sprintf('Correlation Matrix of ROI: %s', rois{ri}.name));
  for ci = 1:size(cond_pairs,1)
    cond = cond_pairs(ci,:);

    t1 = trialAmps{cond(1)}; % nVoxels x nRepeats
    t2 = trialAmps{cond(2)}; % nVoxels x nRepeats
    
    if subtracted
      % Get the corresponding phase-scrambled conditions
      ps1 = cat(2, trialAmps{stimValues(1,:) == stimValues(1,cond(1)) & ~isTexture(stimValues(2,:))}); %nVoxels x 3*nRepeats
      ps2 = cat(2, trialAmps{stimValues(1,:) == stimValues(1,cond(2)) & ~isTexture(stimValues(2,:))}); %nVoxels x 3*nRepeats

      % Take the average over all presentations of the corresponding phase-scrambled condition.
      n1 = mean(ps1,2);
      n2 = mean(ps1,2);
    else
      n1 = 0; n2 = 0;
    end
    
    % Draw 100 samples of prespecified size and take correlations.
    corrs = [];
    for i = 1:100
      idx = randperm(min(size(t2,2), size(t1,2))); % randomly permute the number of observations

      % Take the average over nSmps presentations, and sort voxels by r2.
      % Subtract out the response to the corresponding phase-scrambled images.
      a1 = mean(t1(:, idx(1:nSmps)),2) - n1; % take nSmps observations of condition 1
      a2 = mean(t2(:, idx(nSmps+1 : min(length(idx),2*nSmps))), 2) - n2; % and a different nSmps observations of condition 2

      % Take the nVox voxels with highest r2.
      c = corrcoef(a1, a2); % and correlate them.
      corrs(i) = c(1,2);
    end
    corr_mtx(ri, cond(1), cond(2)) = mean(corrs);

    disppercent(ci / size(cond_pairs,1));
  end
  disppercent(inf);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to get within/between family correlations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function corrs = ER_getCorrStruct(rois, stimValues, condNames, corr_mtx, includeSameSmp)
if ieNotDefined('includeSameSmp'); includeSameSmp=0; end

texIdx = [1 2 3 4]; % which samples are textures (textures = 1-4, noise = 5-7)
isTexture = arrayfun(@(x) any(x==texIdx), stimValues(2,:));

texFams = unique(stimValues(1,:));

%
corrs = {};
for ri = 1:length(rois)
  cm = squeeze(corr_mtx(ri,:,:));
  
  roi = {};
  % Loop through each row of the correlation matrix
  for ti = 1:size(cm,1)
    cmI = cm(ti,:);
    thisFam = stimValues(1,ti);
    thisSmp = stimValues(2,ti);
    thisSmpType = isTexture(ti);
    
    sameFam = (stimValues(1,:) == thisFam);
    sameSmpType = (isTexture == thisSmpType);
    sameSmp = (stimValues(2,:) == thisSmp);
    
    if includeSameSmp
      sameSmp = 0; 
    end

    if thisSmpType == 1
      roi.(sprintf('fam%i_within_textex', thisFam)) = appendToStruc(roi,sprintf('fam%i_within_textex', thisFam), cmI(sameFam & sameSmpType & ~sameSmp));
      roi.(sprintf('fam%i_between_textex', thisFam)) = appendToStruc(roi,sprintf('fam%i_between_textex', thisFam), cmI(~sameFam & sameSmpType));
    
      roi.(sprintf('fam%i_within_texnoise', thisFam)) = appendToStruc(roi,sprintf('fam%i_within_texnoise', thisFam), cmI(sameFam & ~sameSmpType));
      roi.(sprintf('fam%i_between_texnoise', thisFam)) = appendToStruc(roi,sprintf('fam%i_between_texnoise', thisFam), cmI(~sameFam & ~sameSmpType));
    end
    
  end
  corrs.(rois{ri}.name) = roi;
end

corrTypes = {'within_textex', 'between_textex', 'within_texnoise', 'between_texnoise'};

for ci = 1:length(corrTypes)
  corrType = corrTypes{ci};
  corrs.(corrType) = zeros(length(rois), length(texFams));
  corrs.(sprintf('%s_se', corrType)) = zeros(length(rois), length(texFams));
  for ri = 1:length(rois)
    corrStruc = corrs.(rois{ri}.name);
    for i = 1:length(texFams)
      cor = corrStruc.(sprintf('fam%i_%s', i, corrType));
      corrs.(corrType)(ri, i) = mean(cor);
      corrs.(sprintf('%s_se', corrType))(ri,i) = 1.96*std(cor)./length(cor);    
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% function to parse name
%%%%%%%%%%%%%%%%%%%%%%%%%
function [stimValues, condNames] = parseConditionName(stimNames)

if strcmp(class(stimNames), 'char')
  stimNames = {stimNames};
end

nConds = length(strsplit(stimNames{1}, ' and '));
stimValues = zeros(nConds, length(stimNames));
condNames = {};

for ci = 1:length(stimNames)
  conds = strsplit(stimNames{ci}, ' and ');
  for i = 1:nConds
    splits = strsplit(conds{i}, '=');
    stimValues(i,ci) = str2num(splits{2});
    condNames{i} = splits{1};
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to get all combinations of two vectors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cv = combvecs(A,B)
[a,b] = meshgrid(A,B);
c = cat(2,a',b');
cv = reshape(c, [], 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to append vectors to struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function retval = appendToStruc(structure, field, vector)

if isfield(structure, field)
  retval = [structure.(field) vector];
else
  retval = vector;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to plot within/between family Correlations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotCorrelations(corrs, rois, texFams, imNames)
colors = brewermap(length(rois)+1, 'YlGnBu');
colors = colors(2:end, :);

f=figure;
% First plot bar plots averaged across texture families
subplot(2,1,1);
col = colors(3,:);
bar((1:length(rois))-.125, mean(corrs.within_textex,2), .25, 'FaceColor', col);  hold on;
bar((1:length(rois))+.125, mean(corrs.between_textex,2), .25, 'FaceColor', col+(1-col)*(.5)); 
myerrorbar((1:length(rois))-.125, mean(corrs.within_textex,2), 'yError', mean(corrs.within_textex_se,2), 'Symbol', '.');
myerrorbar((1:length(rois))+.125, mean(corrs.between_textex,2), 'yError', mean(corrs.between_textex_se,2), 'Symbol', '.');
ylim([-.1, .7]);
set(gca, 'XTick', 1:length(rois));
set(gca, 'XTickLabel', rois);
title('Within vs Between family texture correlation');
legend({'Within', 'Between'}, 'Location', 'Northwest');
legend('boxoff');
drawPublishAxis('labelFontSize', 14);

% Next plot point plots split out by texture family.
subplot(2,1,2);
handles = [];
for ti = 1:length(texFams)
  col = colors(ti,:);
  h = myerrorbar((1:length(rois))-0.1, corrs.('within_textex')(:,ti), 'yError', corrs.('within_textex_se')(:,ti), 'Symbol', 'o', 'Color', col, 'MarkerSize', 10); hold on;
  handles = [handles h];
  h = myerrorbar((1:length(rois))+0.1, corrs.('between_textex')(:,ti), 'yError', corrs.('between_textex_se')(:,ti), 'Symbol', 'o', 'Color', col+(1-col)*(.5), 'MarkerSize', 10); hold on;
  handles = [handles h];
end
ylim([-.1, .7]);
set(gca, 'XTick', 1:length(rois));
set(gca, 'XTickLabel', rois);
title('Within vs Between family texture correlation, split by texture family');
xlim([0 length(rois)+1]);
labels=[cellfun(@(x) [x ' Within'], imNames, 'un', 0), cellfun(@(x) [x ' Between'], imNames, 'un', 0)];
legend(handles([1:2:end, 2:2:end]), labels, 'Location', 'Northwest');
legend('boxoff');
drawPublishAxis('labelFontSize', 14);
%set(f, 'Position', [21.4, 16.2, 23.1, 21.2]);
set(f, 'Position', [37.65 3.698 29.44 35.18]);
