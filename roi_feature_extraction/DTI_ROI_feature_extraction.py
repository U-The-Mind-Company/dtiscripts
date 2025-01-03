# run this command when these libraries/packages are not installed
# pip install nibabel matplotlib nilearn dipy pandas numpy matplotlib

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import nibabel as nib
from nibabel.processing import resample_from_to
from nilearn import plotting

from dipy.reconst.dti import TensorModel, fractional_anisotropy, mean_diffusivity, radial_diffusivity, axial_diffusivity
from dipy.io.image import load_nifti_data
from dipy.core.gradients import gradient_table

import warnings
warnings.filterwarnings('ignore')


class FeatureExtractor:
    def __init__(self, data_dir, results_dir, atlas_dir, dti_filename, atlas_filename, atlas_labels_filename):
        self.data_dir = data_dir
        self.results_dir = results_dir
        self.atlas_dir = atlas_dir
        self.dti_filename = dti_filename
        self.atlas_filename = atlas_filename

        # set the DTI and atlas file paths -- assuming that the DTI file is preprocessed
        self.bvals_file = os.path.join(self.data_dir, f'{dti_filename}.bval') 
        self.bvecs_file = os.path.join(self.data_dir, f'{dti_filename}.bvec')
        self.img_file = os.path.join(self.data_dir, f'{dti_filename}.nii.gz')  

        self.atlas_file = os.path.join(self.atlas_dir, f'{atlas_filename}.nii.gz')  
        self.atlas_labels_file = os.path.join(self.atlas_dir, f'{atlas_labels_filename}.csv')  


    def load_data(self):
        # Load data
        self.img = nib.load(self.img_file)
        self.atlas = nib.load(self.atlas_file)
        self.bvals = np.loadtxt(self.bvals_file)
        self.bvecs = np.loadtxt(self.bvecs_file)
        self.aal_labels = pd.read_csv(self.atlas_labels_file)

        print(f"Data loaded: DTI img = {self.img.shape}, Atlas = {self.atlas.shape}, B-val = {self.bvals.shape}, B-vec = {self.bvecs.shape}")

    def def_tensor_model(self):
        # gradient_table needs to be constructed when the gtab_file is not directly usable
        gtab = gradient_table(self.bvals, self.bvecs)

        # resample the atlas to the dti image when they have different sizes
        self.atlas_resampled = resample_from_to(self.atlas, self.img.slicer[:,:,:,0], order=0)  # Using nearest-neighbor interpolation

        # Check new shapes
        # print("Resampled atlas shape:", self.atlas_resampled.shape)
        # print("DTI image shape:", self.img.shape)

        # fit the tensor model and compute the diffusion features
        tensor_model = TensorModel(gtab)
        tensor_fit = tensor_model.fit(self.img.get_fdata()) # because img already loaded as nibabel nii obj, need to convert to array

        self.FA = fractional_anisotropy(tensor_fit.evals)
        self.MD = mean_diffusivity(tensor_fit.evals)
        self.RD = radial_diffusivity(tensor_fit.evals)
        self.AD = axial_diffusivity(tensor_fit.evals)

        print("Tensor model built..")

    def extract_save_features(self):
        labels = self.atlas_resampled.get_fdata()  # Assuming atlas data is already loaded correctly
        unique_labels = np.unique(labels)[1:]  # Exclude 0, assuming it's background

        # Prepare DataFrame to hold features
        df = pd.DataFrame(columns=['ROI', 'FA_mean', 'MD_mean', 'RD_mean', 'AD_mean'])

        # Loop through each label and calculate mean DTI metrics
        for roi in unique_labels:
            mask = labels == roi
            print(roi)
            roi_name = self.aal_labels[self.aal_labels['data__label__index'] == roi]['data__label__name'].values[0]

            fa_mean = np.mean(self.FA[mask])
            md_mean = np.mean(self.MD[mask])
            rd_mean = np.mean(self.RD[mask])
            ad_mean = np.mean(self.AD[mask])
            
            # Create a DataFrame for the current ROI's metrics
            roi_df = pd.DataFrame({'ROI': [roi_name], 'FA_mean': [fa_mean], 'MD_mean': [md_mean], 'RD_mean': [rd_mean], 'AD_mean': [ad_mean]})
            
            # Concatenate the current DataFrame with the new row
            df = pd.concat([df, roi_df], ignore_index=True)

        print(f"Size of the extracted features for this DTI image: {df.shape}")
        print(df.head())
        df.to_csv(f'{self.results_dir}/SUB-{self.dti_filename}-dti_metrics_by_roi.csv', index=False)


if __name__ == "__main__":
    """
    File structure: -- main directory
                        -- DTI & Atlas data files directory
                        -- results directory 
    """
    main_dir = "/Users/mansoor/Documents/Projects/DTI-Analysis"
    data_dir = f"{main_dir}/data/preproc/sub_01_post"
    results_dir = f"{main_dir}/dti-scripts/roi_feature_extraction/results"
    atlas_dir = f"{main_dir}/data/templates/"
    
    dti_filename = "sub_01_post"
    atlas_filename = "AAL"
    atlas_labels_filename = "AAL_labels"  # This path may need adjustment

    DTI_feat_extractor = FeatureExtractor(data_dir, results_dir, atlas_dir,
                                           dti_filename, atlas_filename, atlas_labels_filename)
    DTI_feat_extractor.load_data()
    DTI_feat_extractor.def_tensor_model()
    DTI_feat_extractor.extract_save_features()