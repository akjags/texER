function computeStimTraces(rois, saveName, justGetTraces)

if ieNotDefined('rois')
  rois = {'V1', 'V2', 'V3', 'V4'};
end

if ieNotDefined('saveName')
  saveName = '~/proj/texER/stimTraceAmps.mat';
end

if ieNotDefined('justGetTraces')
  justGetTraces = 0;
end

%% Set up view
v = newView;
v = viewSet(v, 'curGroup', 'Concatenation');
v = viewSet(v, 'curScan', 1);

%% Load analysis and get d variable
v = loadAnalysis(v, 'erAnal/erAnal2.mat');

d = viewGet(v, 'd');
overlays = viewGet(v, 'overlays');
r2 = overlays(1).data{1}; % Get r2 map
scanDims = viewGet(v, 'scanDims');
concatInfo = viewGet(v, 'concatInfo');
framePeriod = viewGet(v, 'framePeriod');

% Get traces and condition names
ehdr = d.ehdr;
stimNames = d.stimNames;
[stimValues, condNames] = parseConditionName(stimNames);

% Set r2 cutoff and hemodynamic lag to start trace 
cutoff=.2;
hLag = 2;
hdrlen = 15;

%% Loop through ROIs and extract response to each observation.
stimTraces = [];
stimTraces.stimNames = stimNames; 
stimTraces.condNames = condNames;
stimTraces.stimValues = stimValues;
stimTraces.ehdr = ehdr;
for ri = 1:length(rois)
  roi = loadROITSeries(v, rois{ri});
  tSeries = roi.tSeries;

  % Get high r2 voxels
  ind = sub2ind(scanDims, roi.scanCoords(1,:), roi.scanCoords(2,:), roi.scanCoords(3,:));
  roi_r2 = r2(ind);
  
  if ~justGetTraces
    % Fit gamma to each voxel and get amplitude.
    allAmplitudes = {};
    amplitudes = zeros(length(stimNames), size(tSeries,1));
    ampSTE = zeros(length(stimNames), size(tSeries,1));
    disppercent(-inf, sprintf('Fitting timecourse to %i voxels in ROI %s', size(tSeries,1), roi.name));
    parfor vi = 1:size(tSeries,1)
      timecourse = tSeries(vi,:);

      ftc = fitTimecourse(timecourse, d.stimvol, framePeriod, 'concatInfo', concatInfo, 'fitType=deconv', 'amplitudeType=fit1', 'displayFit=0', 'verbose=0');
      glm = fitTimecourse(timecourse, d.stimvol, framePeriod, 'concatInfo', concatInfo, 'fitType=glm', 'displayFit=0', 'verbose=0', 'returnAllFields=1', 'option=std');
      
      allAmplitudes{vi} = glm.allAmplitudes;
      amplitudes(:,vi) = ftc.amplitude';
      ampSTE(:,vi) = ftc.amplitudeSTE';

      disppercent(vi/size(tSeries,1));
    end
    disppercent(inf);
    stimTraces.(rois{ri}).amplitudes = amplitudes;
    stimTraces.(rois{ri}).amplitudeSTE = ampSTE;
    stimTraces.(rois{ri}).trialAmps = allAmplitudes;
  else
    disp('Skipping fitting gamma functions. Just getting traces');
  end
  stimTraces.(rois{ri}).scanCoords = roi.scanCoords;
  stimTraces.(rois{ri}).r2 = roi_r2;

  
  % Loop through each condition and extract the stimvol traces.
  traces = {};
  for ci = 1:length(stimNames)
    condName = stimNames{ci};
    stimvol = d.stimvol{ci};
    
    ts = [];
    for svi = 1:length(stimvol)
      if hLag+stimvol(svi)+hdrlen > size(tSeries,2)
        continue;
      else
        ts = cat(3, ts, tSeries(:, (hLag+stimvol(svi)):(hLag+stimvol(svi)+hdrlen)));
      end
    end
    % Keep only the voxels whose r2 exceeds the specified cutoff.
    %traces{ci} = ts(roi_r2>cutoff,:,:);
    traces{ci} = ts;
  end

  stimTraces.(rois{ri}).traces = traces;
end

save(saveName, '-struct', 'stimTraces');
keyboard
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
