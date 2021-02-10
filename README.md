[![Abcdspec-compliant](https://img.shields.io/badge/ABCD_Spec-v1.1-green.svg)](https://github.com/brain-life/abcd-spec)
[![Run on Brainlife.io](https://img.shields.io/badge/Brainlife-brainlife.app.185-blue.svg)](https://doi.org/10.25663/brainlife.app.185)

# Tract Analysis Profiles 

This app will compute profiles of tensor (i.e. AD, FA, MD, RD) and/or NODDI (i.e. ICVF, ISOVF, OD) measures, and their inverse measures and corresponding standard deviations, over a user-specified number of nodes along segmented white matter fascicles. This app requires the white matter classification (wmc) datatype, the whole-brain tractogram tck datatype used to generate the wmc dattype, and a tensor dataytpe. Optionally, the user can input a NODDI dataype. This version of the app uses Dipy and pyAFQ to generate tract profiles. EPS and PNG images of these profiles are created for the measures, and a .csv file is outputted containing all the measures and standard deviations for each segmented fascicle. 

### Authors 

- Brad Caron (bacaron@iu.edu) 

### Contributors 

- Soichi Hayashi (hayashis@iu.edu) 

### Funding 

[![NSF-BCS-1734853](https://img.shields.io/badge/NSF_BCS-1734853-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=1734853)
[![NSF-BCS-1636893](https://img.shields.io/badge/NSF_BCS-1636893-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=1636893)
[![NSF-ACI-1916518](https://img.shields.io/badge/NSF_ACI-1916518-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=1916518)
[![NSF-IIS-1912270](https://img.shields.io/badge/NSF_IIS-1912270-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=1912270)
[![NIH-NIBIB-R01EB029272](https://img.shields.io/badge/NIH_NIBIB-R01EB029272-green.svg)](https://grantome.com/grant/NIH/R01-EB029272-01)

### Citations 

Please cite the following articles when publishing papers that used data, code or other resources created by the brainlife.io community. 

1. Yeatman JD, Dougherty RF, Myall NJ, Wandell BA, Feldman HM (2012) Tract Profiles of White Matter Properties: Automating Fiber-Tract Quantification. PLoS ONE 7(11): e49790. https://doi.org/10.1371/journal.pone.0049790
1. Garyfallidis E, Brett M, Amirbekian B, Rokem A, van der Walt S, Descoteaux M, Nimmo-Smith I and Dipy Contributors (2014). DIPY, a library for the analysis of diffusion MRI data. Frontiers in Neuroinformatics, vol.8, no.8.
1. Avesani, P., McPherson, B., Hayashi, S. et al. The open diffusion data derivatives, brain data upcycling via integrated publishing of derivatives and reproducible open cloud services. Sci Data 6, 69 (2019). https://doi.org/10.1038/s41597-019-0073-y

## Running the App 

### On Brainlife.io 

You can submit this App online at [https://doi.org/10.25663/brainlife.app.185](https://doi.org/10.25663/brainlife.app.185) via the 'Execute' tab. 

### Running Locally (on your machine) 

1. git clone this repo 

2. Inside the cloned directory, create `config.json` with something like the following content with paths to your input files. 

```json 
{
   "classification":    "testdata/wmc/classification.mat",
   "fa":    "testdata/tensor/ifa.nii.gz",
   "md":    "testdata/tensor/md.nii.gz",
   "rd":    "tesdata/tensor/rd.nii.gz",
   "ad":    "tesdata/tensor/ad.nii.gz",
   "icvf":    "testdata/noddi/icvf.nii.gz",
   "isovf":    "testdata/noddi/isovf.nii.gz",
   "od":    "tesdata/noddi/od.nii.gz",
   "track":    "testdata/track/track.tck",
   "num_nodes":    200
} 
``` 

### Sample Datasets 

You can download sample datasets from Brainlife using [Brainlife CLI](https://github.com/brain-life/cli). 

```
npm install -g brainlife 
bl login 
mkdir input 
bl dataset download 
``` 

3. Launch the App by executing 'main' 

```bash 
./main 
``` 

## Output 

The main output of this App is contains one CSV file for all measure data for each track and PNG and EPS images of the profiles. The CSV files are formatted in the following format: ad_mean, ad_sd, fa_mean, fa_sd, md_mean, md_sd, rd_mean, rd_sd, ndi_mean, ndi_sd, isovf_mean, isovf_sd, odi_mean, odi_sd.

#### Product.json 

The secondary output of this app is `product.json`. This file allows web interfaces, DB and API calls on the results of the processing. 

### Dependencies 

This App requires the following libraries when run locally. 

  - singularity: https://singularity.lbl.gov/
