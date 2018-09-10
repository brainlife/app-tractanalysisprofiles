
[![Run on Brainlife.io](https://img.shields.io/badge/Brainlife-bl.app.43-blue.svg)](https://doi.org/10.25663/bl.app.43)

# app-tractanalysisprofiles
This service computes profiles of tensor (i.e. AD, FA, MD, RD) and/or NODDI (i.e. ICVF, ISOVF, OD) measures, and their inverse measures and corresponding standard deviations, over a user-specified number of nodes along segmented white matter fascicles. Based on user-input, it computes these profiles using either dtiComputeDiffusionPropertiesAlongFG_sd from Vistasoft, which weights streamlines based on their distance away from a central core of the fascicle (i.e. 'volume based'), or Compute_FA_AlongFG, which does not weight the streamlines and is appropriate for testing models involving streamline/fiber corssings (i.e. 'fiber based'). EPS and PNG images of these profiles are created for the measures, and a .csv file is outputted containing all the measures and standard deviations for each segmented fascicle.

### Authors
- Brad Caron (bacaron@iu.edu)
- Lindsey Kitchell (kitchell@iu.edu)
- Soichi Hayashi (hayashi@iu.edu)
- Franco Pestilli (franpest@indiana.edu)

### Contributors
- Soichi Hayashi (hayashi@iu.edu)
- Franco Pestilli (franpest@indiana.edu)

### Funding
[![NSF-BCS-1734853](https://img.shields.io/badge/NSF_BCS-1734853-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=1734853)
[![NSF-BCS-1636893](https://img.shields.io/badge/NSF_BCS-1636893-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=1636893)

## Running the App 

### On Brainlife.io

You can submit this App online at [https://doi.org/10.25663/bl.app.43](https://doi.org/10.25663/bl.app.43) via the "Execute" tab.

### Running Locally (on your machine)

1. git clone this repo.
2. Inside the cloned directory, create `config.json` with something like the following content with paths to your input files.

```json
{
        "tensor": "./input/tensor/",
      	"noddi": "./input/noddi/",
        "afq": "./input/wmc/output.mat",
        "numnodes": 100,
        "fiberbased": False
}
```

### Sample Datasets

You can download sample datasets from Brainlife using [Brainlife CLI](https://github.com/brain-life/cli).

```
npm install -g brainlife
bl login
mkdir input
bl dataset download 5b09308ce81b47020407c4e9 && mv 5b09308ce81b47020407c4e9 input/tensor
bl dataset download 5b0a713ce81b47020407c739 && mv 5b0a713ce81b47020407c739 input/noddi
bl dataset download 5b0dbdcfe81b47020407ca17 && mv 5b0dbdcfe81b47020407ca17 input/wmc

```


3. Launch the App by executing `main`

```bash
./main
```

## Output

The two main outputs of this App are folders called "images" and "tractprofile". The "images" folder contains .eps and .png images of the tract profiles for every measure and fascicle. The "tractprofile" folder contains a .csv file for every fascicle. The .csv file follows this layout:

```
ad_1 (measure data) ad_2 (measure std) fa_1 (measure data) fa_2 (measure std) md_1 (measure data) md_2 (measure std) rd_1 (measure data) rd_2 (measure std) ad_inverse_1 (measure data) ad_inverse_2 (measure std) fa_inverse_1 (measure data) fa_inverse_2 (measure std) md_inverse_1 (measure data) md_inverse_2 (measure std) rd_inverse_1 (measure data) rd_inverse_2 (measure std) icvf_1 (measure data) icvf_2 (measure std) isovf_1 (measure data) isovf_2 (measure std) od_1 (measure data) od_2 (measure std) icvf_inverse_1 (measure data) icvf_inverse_2 (measure std) isovf_inverse_1 (measure data) isovf_inverse_2 (measure std) od_inverse_1 (measure data) od_inverse_2 (measure std)
```

#### Product.json
The secondary output of this app is `product.json`. This file allows web interfaces, DB and API calls on the results of the processing. 

### Dependencies

This App requires the following libraries when run locally.

  - singularity: https://singularity.lbl.gov/
  - VISTASOFT: https://github.com/vistalab/vistasoft/
  - SPM 8: https://www.fil.ion.ucl.ac.uk/spm/software/spm8/
  - jsonlab: https://github.com/fangq/jsonlab.git
