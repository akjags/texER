function compute_splithalf_errorbars(stimTraces)
roiName = 'V1';
rois = {'V1', 'V2', 'V3', 'V4'};
condI = 1;
nRepeats = 100;
bootstrap_size = 100;

nObservations = [2,4,6,8,11];
%kObs = kObs(2);

%% Step 1: Simulate N traces by bootstrapping the original traces.
disp(sprintf('Creating %i simulated traces per condition by bootstrapping', bootstrap_size));
simTraces = [];
for ri = 1:length(rois)
  roiName = rois{ri};
  roi_traces = cellfun(@(x) squeeze(mean(x, 2)), stimTraces.(sprintf('%s_traces', roiName)), 'un', 0);
  simTraces.(sprintf('%s_traces', roiName)) = {};
  % Loop through each condition and sample observations with replacement.
  for ci = 1:length(roi_traces)
    trace = roi_traces{ci};
    % sample with replacement to create a simulated dataset.
    simTraces.(sprintf('%s_traces', roiName)){ci} = datasample(trace, 10, bootstrap_size, 'Replace', true);
  end
end
    


%% Step 2: Go through the simulated traces and compute the correlation
% matrix.

nConditions = length(stimTraces.stimNames);
cond_pairs = combvec(1:nConditions, 1:nConditions);

corr_mtx = zeros(length(nObservations), length(rois), nConditions, nConditions, nRepeats);
disp(sprintf('Computing %i correlation matrix', nRepeats));
delete(gcp('nocreate')); parpool(10);
tic;
parfor repI = 1:nRepeats
  c_mtx = zeros(length(nObservations), length(rois), nConditions,nConditions);
  
  for oi = 1:length(nObservations)
    kObs = nObservations(oi);

    for ri = 1:length(rois)
      % Get all pairs of conditions and loop through.


      roi_traces = simTraces.(sprintf('%s_traces', roiName));
      
      % For each pair of conditions, sample k observations.
      for ci = 1:size(cond_pairs,1)
        conds = cond_pairs(ci,:);

        condI = roi_traces{conds(1)}; condJ = roi_traces{conds(2)};

        corrs = zeros(1,nRepeats);
        for i = 1:nRepeats
          idx = randperm(min(size(condI,2), size(condJ,2))); % randomly permute the number of observations
          I_idx = idx(1:kObs); J_idx = idx(kObs+1: 2*kObs);
          % Randomly sample k observations from each conditions and correlate
          smpI = mean(condI(:,I_idx),2); 
          smpJ = mean(condJ(:,J_idx),2);

          c = corrcoef(smpI, smpJ);
          corrs(i) = c(1,2);
        end
        c_mtx(oi,ri,conds(1), conds(2)) = mean(corrs(i));
        %
      end
    end
  end
  disp(sprintf('%i complete', repI));
  corr_mtx(:,:, :,:, repI) = c_mtx;
end
toc;

%%%

%% Step 3: Extract from the corr_mtx the within and between family
% correlations and use these to compute confidence intervals.
conf_intrvl = zeros(length(nObservations), length(rois), 4,2);

for oi = 1:length(nObservations)
  kObs = nObservations(oi);

  for ri = 1:length(rois)
    % Get corr_mtx
    cm = squeeze(corr_mtx(oi,ri,:,:,:));
    cm_sort = sort(cm, 3, 'descend');
    
    cond1 = []; cond2 = []; cond3 = []; cond4 = [];
    for ci = 1:size(cond_pairs,1)
      conds = cond_pairs(ci,:);

      fam1 = stimTraces.stimValues(1, conds(1)); fam2 = stimTraces.stimValues(1, conds(2));
      smp1 = stimTraces.stimValues(2, conds(1)); smp2 = stimTraces.stimValues(2, conds(2));

      if fam1 == fam2 % Within-Family
        if smp1 < 5 && smp2 < 5 % tex-tex
          cond1 = [cond1 squeeze(cm(conds(1), conds(2), :))];

        elseif xor(smp1<5, smp2<5) % tex-noise.
          cond2 = [cond2 squeeze(cm(conds(1), conds(2), :))];
        end

      else % Between Family
        if smp1 < 5 && smp2 < 5 % tex-tex
          cond3 = [cond3 squeeze(cm(conds(1), conds(2), :))];
        elseif xor(smp1<5, smp2<5) % tex-noise.
          cond4 = [cond4 squeeze(cm(conds(1), conds(2), :))];
        end
      end
    end

    for i = 1:4
      confInt = mean(eval(sprintf('cond%i', i)),2);
      ci_sort = sort(confInt, 'descend');
      conf_intrvl(oi,ri,i,:) = [median(confInt)- ci_sort(ceil(.975*nRepeats)), ci_sort(ceil(.025*nRepeats)) - median(confInt)];
    end
  end
end

confidence_interval = [];
confidence_interval.conditions = {'Within family texture-texture', 'Within family texture-noise', 'Between family texture-texture', 'Between family texture-noise'};
confidence_interval.conf_intrvl = conf_intrvl;
save('~/proj/texER/confidence_interval.mat', '-struct', 'confidence_interval');

%%
function cv = combvec(A,B)
[a,b] = meshgrid(A,B);
c = cat(2,a',b');
cv = reshape(c, [], 2);