import numpy as np
import matplotlib.pyplot as plt 
import h5py, os
from sklearn.cross_decomposition import PLSRegression
from sklearn.decomposition import PCA


def load_rois_hdf5(hdf_file='s06420190206/rois_pool4_glm.mat'):
    if not os.path.isfile(hdf_file):
        print 'file %s not found' % hdf_file

    if not os.access(hdf_file, os.R_OK):
        print 'file %s not readable' % hdf_file

    f = h5py.File(hdf_file, 'r')
    f.keys()

    # Load and extract the Instances for each ROI
    rois_mat = f

    rois = {}
    for column in rois_mat['rois']:
        roi = rois_mat[column[0]]
        instances = roi['instance']['instances']
        amps = np.array(roi['instance']['fit']['amplitude'])
        r2 = np.array(roi['instance']['r2'])[0]


        inst = {};
        for i,col in enumerate(instances):
            inst[i] = np.array(rois_mat[col[0]])
        roi_name = rois_mat[column[0]]['name'][:].tostring().decode('utf-16')
        rois[roi_name] = {};
        rois[roi_name]['instances'] = inst;
        rois[roi_name]['amplitudes'] = amps
        rois[roi_name]['r2'] = r2
    stimValues = np.array(rois_mat['stimValues'])

    return rois, stimValues

## FUNCTIONS to get voxel responses and get gram matrix features.
gm_dir = '/Users/akshay/proj/TextureSynthesis/behav_data/gram_mtx_bw'
imNames = ['bricks', 'bark', 'rocks', 'glass']
texSmps = [1,2,3,4]
def get_voxel_responses(rois):
    responses = {};
    for roi in rois.keys():
        inst = rois[roi]['instances']

        vox_resp = [];
        for I in inst.keys():
            smp = int(stimValues[I,1])            
            if smp in texSmps:
                instance = inst[I]
                vox_resp.append(instance)

        responses[roi] = np.array(vox_resp)

    return responses

def get_gram_matrices(rois, getAll=0, gm_dir = gm_dir, obs_lay='pool4', obs_rf = '1x1'):
    for roi in rois.keys():
        inst=rois[roi]['instances']

        gms = []
        for I in inst.keys():
            smp = int(stimValues[I,1])
            fam = int(stimValues[I,0])-1
            
            if smp in texSmps:
                gm = np.load('{}/gram{}_1x1_pool4_{}_smp{}.npy'.format(gm_dir, obs_rf, imNames[fam], smp)).item()
                gms.append(gm[obs_lay].ravel())
        break
    return np.array(gms)

## Extract voxel responses of specified ROI and gram matrices of specified layer.
def get_features_and_responses(rois, obs_lay='pool2', roi='V3', avg_reps=1, runPCA=1):
    responses = get_voxel_responses(rois)
    gms = get_gram_matrices(rois, obs_lay=obs_lay)

    # Can specify whether to average repeated stimuli or to keep them separate (to have more stimuli)
    if avg_reps:
        roi_resp = np.mean(responses[roi],axis=2)
        gram_matrices = gms
    else:
        # Separate out responses to each sample
        roi_resp0 = responses[roi]
        a = [roi_resp0[:,:,i] for i in range(roi_resp0.shape[2])]
        roi_resp = np.concatenate(a, axis=0)
        gram_matrices = np.concatenate((gms,)*12,axis=0)

    ## Perform PCA on the features
    if runPCA:
        pca = np.load('/Users/akshay/proj/TextureSynthesis/texpca/pool2_dims.npy').item()['pca']
        gram_matrices = pca.transform(gram_matrices)

    return roi_resp, gram_matrices


## Run PLS regression on each voxel (leave-one-image-out)
def predict_voxel_responses(roi_resp, gm):
    nVox = roi_resp.shape[1]; nIms = roi_resp.shape[0]
    cv_r2 = np.zeros(nVox);
    cv_preds = np.zeros((nVox, nIms))
    for vi in range(nVox):
        vox_resp = roi_resp[:,vi]

        for i in range(len(vox_resp)):
            pls = PLSRegression(n_components=10)

            # Fit on all images except the i'th
            all_but1 = np.setdiff1d(np.arange(len(vox_resp)), i)
            pls.fit(gm[all_but1, :], vox_resp[all_but1])

            # Test on the i'th image
            cv_preds[vi,i] = pls.predict(gm[i,:].reshape(1,-1))

        # compute r2
        cv_r2[vi] = compute_r2(vox_resp, cv_preds[vi,:])
        
    return cv_preds, cv_r2

def compute_r2(y_true, y_pred):
    sse = ((y_true - y_pred)**2).sum()
    ssto = ((y_true - np.mean(y_true))**2).sum()
    return 1 - (sse / ssto)



