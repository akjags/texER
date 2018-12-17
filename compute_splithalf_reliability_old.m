function stimTraces = compute_splithalf_reliability(stimTraces)

if ieNotDefined('stimTraces')
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

  % Get traces and condition names
  ehdr = d.ehdr;
  stimNames = d.stimNames;
  [stimValues, condNames] = parseConditionName(stimNames);

  cutoff=.2;
  hdrlen = 15;
  
  %% Loop through ROIs and extract response to each observation.
  rois = {'V1', 'V2', 'V3', 'V4'};

  stimTraces = [];
  stimTraces.stimNames = stimNames; stimTraces.condNames = condNames;
  stimTraces.stimValues = stimValues;

  for ri = 1:length(rois)
    roi = loadROITSeries(v, rois{ri});
    tSeries = roi.tSeries;

    traces = {};
    
    % Get high r2 voxels
    ind = sub2ind(scanDims, roi.scanCoords(1,:), roi.scanCoords(2,:), roi.scanCoords(3,:));
    roi_r2 = r2(ind);

    % Loop through each condition and extract the stimvol traces.
    for ci = 1:length(stimNames)
      condName = stimNames{ci};
      stimvol = d.stimvol{ci};

      ts = [];
      for svi = 1:length(stimvol)
        if stimvol(svi)+hdrlen > size(tSeries,2)
          continue;
        else
          ts = cat(3, ts, tSeries(:, stimvol(svi):stimvol(svi)+hdrlen));
        end
      end
      traces{ci} = ts(roi_r2>cutoff,:,:);
    end

    stimTraces.(sprintf('%s_traces', rois{ri})) = traces;
  end
end
%save('~/proj/texER/stimTraces.mat', '-struct', 'stimTraces');
%% Compute and plot correlation matrix.
conds = combvec(1:length(stimTraces.stimNames), 1:length(stimTraces.stimNames));

for nSmps = [2,4,6,8,11]
  disp(sprintf('Computing correlation matrix with nObservations = %i', nSmps*2)); 
  rois = {'V1', 'V2', 'V3', 'V4'};
  for ri = 1:length(rois)
    traces = stimTraces.(sprintf('%s_traces', rois{ri}));

    corr_mtx = zeros(length(stimTraces.stimNames), length(stimTraces.stimNames));
    within_fam_tex_corr = []; between_fam_tex_corr = []; same_sample_tex_corr = [];
    which_tex_fams = []; which_cond_type = []; tex_corr = []; 
    within_fam_tex_noise_corr = []; between_fam_tex_noise_corr = [];
    disppercent(-inf, sprintf('ROI: %s', rois{ri}));
    for ci = 1:size(conds,1)
      cond = conds(ci,:);

      t1 = squeeze(mean(traces{cond(1)},2));
      t2 = squeeze(mean(traces{cond(2)},2));

      % Draw 1000 samples of prespecified size and take correlations.
      corrs = [];
      for i = 1:100
        idx = randperm(min(size(t2,2), size(t1,2)));

        a1 = mean(t1(:, idx(1:nSmps)),2);
        a2 = mean(t2(:, idx(nSmps+1 : 2*nSmps)), 2);

        c = corrcoef(a1, a2);
        corrs(i) = c(1,2);
        
      end
      corr_mtx(cond(1), cond(2)) = mean(corrs);
      
      which_tex_fams = [which_tex_fams; stimTraces.stimValues(1, cond(1)) stimTraces.stimValues(1, cond(2))];
      tex_corr = [tex_corr; corrs];
      % Get the within and between family correlations, for the conditions
      % in which sample shown is a texture (not noise) sample.
      if stimTraces.stimValues(1,cond(1)) == stimTraces.stimValues(1, cond(2)) && stimTraces.stimValues(2,cond(1)) < 5 && stimTraces.stimValues(2,cond(2)) < 5
        if stimTraces.stimValues(2, cond(1)) == stimTraces.stimValues(2, cond(2)) % same sample
          which_cond_type = [which_cond_type 0];
          same_sample_tex_corr = [same_sample_tex_corr; corrs];
        else % Different samples, but same family.
          which_cond_type = [which_cond_type 1];
          within_fam_tex_corr = [within_fam_tex_corr; corrs];
        end
      elseif stimTraces.stimValues(1,cond(1)) ~= stimTraces.stimValues(1, cond(2)) && stimTraces.stimValues(2,cond(1)) < 5 && stimTraces.stimValues(2,cond(2)) < 5
        which_cond_type = [which_cond_type 2];
        between_fam_tex_corr = [between_fam_tex_corr; corrs];
      elseif xor(stimTraces.stimValues(2,cond(1)) < 5, stimTraces.stimValues(2, cond(2))<5) % if one of the samples is noise.
        if stimTraces.stimValues(1,cond(1)) == stimTraces.stimValues(1,cond(2)) % same family
          which_cond_type = [which_cond_type 3];
          within_fam_tex_noise_corr = [within_fam_tex_noise_corr; corrs];
        else
          which_cond_type = [which_cond_type 4];
          between_fam_tex_noise_corr = [between_fam_tex_noise_corr; corrs];
        end
      end
      disppercent(ci / size(conds,1));
    end
    disppercent(inf);

    stimTraces.(rois{ri}).cond_types = {'same sample', 'within family', 'between family', 'within family noise', 'between family noise'};
    stimTraces.(rois{ri}).which_tex_fams = which_tex_fams;
    stimTraces.(rois{ri}).which_cond_type = which_cond_type;
    stimTraces.(rois{ri}).(sprintf('tex_corr_%i', 2*nSmps)) = tex_corr;
    stimTraces.(rois{ri}).(sprintf('corr_mtx_%i', 2*nSmps)) = corr_mtx;
    stimTraces.(rois{ri}).(sprintf('same_sample_tex_corr_%i', 2*nSmps)) = same_sample_tex_corr;
    stimTraces.(rois{ri}).(sprintf('within_fam_tex_corr_%i', 2*nSmps)) = within_fam_tex_corr;
    stimTraces.(rois{ri}).(sprintf('between_fam_tex_corr_%i', 2*nSmps)) = between_fam_tex_corr;
  end

  % Plot correlation matrix.
  figure; 
  for i = 1:4
    subplot(2,2,i);
    imagesc(stimTraces.(rois{ri}).(sprintf('corr_mtx_%i', 2*nSmps)));
    colormap('hot');
    colorbar;
    title(sprintf('%s Correlation Matrix, nObs = %g', rois{i}, nSmps*2), 'FontSize', 18, 'FontName', 'Gill Sans');
  end
  
end
%%  Plot within and between family split half correlations as a function of number of repeats.
nObsSizes = [2,4,6,8,11];
wf = []; bf = [];
for oi = 1:length(nObsSizes)
  nObs = nObsSizes(oi);
  for ri = 1:length(rois)
    roi = rois{ri};
    
    same_samp = stimTraces.(roi).(sprintf('same_sample_tex_corr_%i', nObs*2));
    within_fam = stimTraces.(roi).(sprintf('%s_within_fam_tex_corr_%i', nObs*2));
    between_fam = stimTraces.(roi).(sprintf('%s_between_fam_tex_corr_%i', nObs*2));
    
    % Get mean correlation for same-sample, within-family, and
    % between-family.
    ss(oi, ri) = mean(mean(same_samp,1));
    wf(oi, ri) = mean(mean(within_fam,1));
    bf(oi, ri) = mean(mean(between_fam,1));
    
    sort_same = sort(mean(same_samp,1),2, 'descend');
    ss_err(oi,ri) = 0.5*(sort_same(25) - sort_same(975));
    ss_errs(oi,ri,:) = [sort_same(25) - median(sort_same), median(sort_same) - sort_same(75)];
    %ss_err(oi,ri) = 1.96*std(mean(same_samp,1)) / sqrt(length(mean(same_samp,1)));
    
    sort_within = sort(mean(within_fam,1),2, 'descend');
    wf_err(oi,ri) = 0.5*(sort_within(25) - sort_within(975));
    wf_errs(oi,ri,:) = [sort_within(25) - median(sort_within), median(sort_within) - sort_within(75)];
    %wf_err(oi,ri) = 1.96*std(mean(within_fam,1)) / sqrt(length(mean(within_fam,1)));
    
    sort_between = sort(mean(between_fam,1),2, 'descend');
    bf_err(oi,ri) = 0.5*(sort_between(25) - sort_between(975));
    bf_errs(oi,ri,:) = [sort_between(25) - median(sort_between), median(sort_between) - sort_between(75)];
    %bf_err(oi,ri) = 1.96*std(mean(between_fam,1)) / sqrt(length(mean(between_fam,1)));
  end
end

%% Plot 
colors = brewermap(5, 'YlGnBu');
colors = colors(2:end, :);

f = figure; set(gcf, 'Color', [1 1 1]);
%set(gcf, 'Position', [78.019, 4.015, 27.8967, 32.7928]);

subplot(3,1,1); colormap('cool');
lineProps = []; lineProps.col = {colors(1,:), colors(2,:), colors(3,:), colors(4,:)};
lineProps.style = {':', ':', ':', ':'}; 
%mseb(2*nObsSizes, ss', ss_err', lineProps); hold on;

xJit = rand(1,4)-.5;
for i = 1:4
  plot(2*nObsSizes + xJit(i), ss(:,i), '.', 'Color', colors(i,:), 'MarkerSize', 20); hold on;
  myerrorbar(2*nObsSizes +xJit(i), ss(:,i), 'yLow', ss(:,i)-ss_errs(:,i,2), 'yHigh', ss(:,i)+ss_errs(:,i,1), 'Symbol=:', 'Color', colors(i,:));
end

xlabel('Number of Repeated Presentations', 'FontSize', 14);
ylabel('Correlation', 'FontSize', 14); ylim([-.01 .25]);
title('Same Sample Split-Half Correlation as a function of number of repeats', 'FontSize', 18, 'FontName', 'Gill Sans');
%legend(rois, 'Location', 'northwest'); 
drawPublishAxis

subplot(3,1,2); colormap('cool');

%mseb(2*nObsSizes, wf', wf_err', lineProps); hold on;
for i = 1:4
  plot(2*nObsSizes + xJit(i), wf(:,i), '.', 'Color', colors(i,:), 'MarkerSize', 20); hold on;
  myerrorbar(2*nObsSizes+xJit(i), wf(:,i), 'yLow', wf(:,i)-wf_errs(:,i,2), 'yHigh', wf(:,i)+wf_errs(:,i,1), 'Symbol=:', 'Color', colors(i,:));
end
xlabel('Number of Repeated Presentations');
ylabel('Correlation', 'FontSize', 14); ylim([-.01 .25]);
title('Within Family Split-Half Correlation as a function of number of repeats');
drawPublishAxis

subplot(3,1,3); 
%mseb(2*nObsSizes, bf', bf_err', lineProps); hold on;
for i = 1:4
  plot(2*nObsSizes+xJit(i), bf(:,i), '.', 'Color', colors(i,:), 'MarkerSize', 20); hold on;
  myerrorbar(2*nObsSizes+xJit(i), bf(:,i), 'yLow', bf(:,i)-bf_errs(:,i,2), 'yHigh', bf(:,i)+bf_errs(:,i,1),'Symbol=:', 'Color', colors(i,:));
end
xlabel('Number of Repeated Presentations');
ylabel('Correlation');  ylim([-.01 .25]);
title('Between Family Split-Half Correlation as a function of number of repeats');
colormap('Winter'); drawPublishAxis


%savepdf(f, '~/proj/texER/splithalf_corr.pdf')
%%
keyboard

%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%
function cv = combvec(A,B)
[a,b] = meshgrid(A,B);
c = cat(2,a',b');
cv = reshape(c, [], 2);

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
