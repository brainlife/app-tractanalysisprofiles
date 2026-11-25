<h1 align="center">Tract Analysis Profiles</h1>

<p align="center">
  Tract profiles of diffusion, NODDI, myelin or quantitative MRI measures along every segmented white matter tract, with plots and a tractmeasures CSV.
</p>

<p align="center">
  <a href="https://github.com/brainlife/abcd-spec"><img src="https://img.shields.io/badge/ABCD_Spec-v1.1-green.svg" alt="ABCD Spec v1.1"></a>
  <a href="https://doi.org/10.25663/brainlife.app.345"><img src="https://img.shields.io/badge/Brainlife-brainlife.app.345-blue.svg" alt="Run on brainlife.io"></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/License-MIT-lightgrey.svg" alt="MIT license"></a>
</p>

> **This branch, `1.6`, runs 2 brainlife.io Apps** that share the same code:
>
> - Tract Analysis Profiles (NODDI) - Deprecated ([doi:10.25663/brainlife.app.345](https://doi.org/10.25663/brainlife.app.345)) (deprecated)
> - Tract Analysis Profiles - Deprecated ([doi:10.25663/brainlife.app.185](https://doi.org/10.25663/brainlife.app.185)) (deprecated)
>
> The other Apps in this repository are listed under **Versions on brainlife.io** below.

## At a glance

<!-- brainlife-app:begin (generated from the App manifest; edit the manifest, not this block) -->

### Tract Analysis Profiles (NODDI) - Deprecated

> **Create plots of diffusion metrics (i.e.**

Create plots of diffusion metrics (i.e. FA, MD, RD, AD) for each of the segmented tracts from AFQ, known as Tract Profiles. Obtains streamline positions from segmented tracts and plots the metrics of interest along "nodes" of the tract, allowing for comparison of individual subject tracts. Requires the dt6 output from dtiinit and a white matter classification output from AFQ or WMA

| | |
|---|---|
| **Inputs** | `neuro/wmc` · `neuro/noddi-deprecated` · `neuro/track/tck` |
| **Outputs** | `generic/images` · `neuro/tractprofile-deprecated` |
| **Container** | `brainlife/mcr:neurodebian1604-r2017a` |
| **Resources** | 4 CPUs · 28 GB memory · 6 h walltime |
| **Run it** | [doi.org/10.25663/brainlife.app.345](https://doi.org/10.25663/brainlife.app.345) |

#### Inputs

| Input | Datatype | Required | Description |
|---|---|---|---|
| `afq` | `neuro/wmc` | yes |  |
| `noddi` | `neuro/noddi-deprecated` | yes |  |
| `tck` | `neuro/track/tck` | yes |  |

#### Outputs

| Output | Datatype | Description |
|---|---|---|
| `images` | `generic/images` |  |
| `profiles` | `neuro/tractprofile-deprecated` |  |

<details>
<summary><b>Parameters</b> (2)</summary>

| Parameter | Type | Default | Description |
|---|---|---|---|
| `fiberbased` | boolean | `false` | to create profiles using a non-weighted, fiber-based method: "true" |
| `numnodes` | number | `200` |  |

</details>

### Tract Analysis Profiles - Deprecated

> **Create plots of diffusion metrics (i.e.**

Create plots of diffusion metrics (i.e. FA, MD, RD, AD) for each of the segmented tracts from AFQ, known as Tract Profiles. Obtains streamline positions from segmented tracts and plots the metrics of interest along "nodes" of the tract, allowing for comparison of individual subject tracts. Requires the dt6 output from dtiinit and a white matter classification output from AFQ or WMA

| | |
|---|---|
| **Inputs** | `neuro/tensor` · `neuro/wmc` · `neuro/noddi-deprecated` (optional) · `neuro/track/tck` |
| **Outputs** | `generic/images` · `neuro/tractprofile-deprecated` |
| **Container** | `brainlife/mcr:neurodebian1604-r2017a` |
| **Resources** | 4 CPUs · 28 GB memory · 6 h walltime |
| **Run it** | [doi.org/10.25663/brainlife.app.185](https://doi.org/10.25663/brainlife.app.185) |

#### Inputs

| Input | Datatype | Required | Description |
|---|---|---|---|
| `tensor` | `neuro/tensor` | yes |  |
| `afq` | `neuro/wmc` | yes |  |
| `noddi` | `neuro/noddi-deprecated` | no |  |
| `tck` | `neuro/track/tck` | yes |  |

#### Outputs

| Output | Datatype | Description |
|---|---|---|
| `images` | `generic/images` |  |
| `profiles` | `neuro/tractprofile-deprecated` |  |

<details>
<summary><b>Parameters</b> (2)</summary>

| Parameter | Type | Default | Description |
|---|---|---|---|
| `fiberbased` | boolean | `false` | to create profiles using a non-weighted, fiber-based method: "true" |
| `numnodes` | number | `200` |  |

</details>

<!-- brainlife-app:end -->

## Versions on brainlife.io

Each version is a separate App on brainlife.io, built from a branch of this repository. The Apps on branch `1.13` share the same code and differ only in the measures they ask for.

| Version | Method | Branch | Run it | Status |
|---|---|---|---|---|
| Tract Analysis Profiles - Myelin mapping | Vistasoft | [`1.13`](https://github.com/brainlife/app-tractanalysisprofiles/tree/1.13) | [doi:10.25663/brainlife.app.523](https://doi.org/10.25663/brainlife.app.523) | current |
| Tract Analysis Profiles - Dipy | DIPY | [`dipy-1.0`](https://github.com/brainlife/app-tractanalysisprofiles/tree/dipy-1.0) | [doi:10.25663/brainlife.app.474](https://doi.org/10.25663/brainlife.app.474) | available |
| Tract Analysis Profiles (NODDI) | Vistasoft | [`1.13`](https://github.com/brainlife/app-tractanalysisprofiles/tree/1.13) | [doi:10.25663/brainlife.app.362](https://doi.org/10.25663/brainlife.app.362) | available |
| Tract Analysis Profiles | Vistasoft | [`1.13`](https://github.com/brainlife/app-tractanalysisprofiles/tree/1.13) | [doi:10.25663/brainlife.app.361](https://doi.org/10.25663/brainlife.app.361) | available |
| **Tract Analysis Profiles (NODDI) - Deprecated** (this branch) | Vistasoft | [`1.6`](https://github.com/brainlife/app-tractanalysisprofiles/tree/1.6) | [doi:10.25663/brainlife.app.345](https://doi.org/10.25663/brainlife.app.345) | deprecated |
| **Tract Analysis Profiles - Deprecated** (this branch) | Vistasoft | [`1.6`](https://github.com/brainlife/app-tractanalysisprofiles/tree/1.6) | [doi:10.25663/brainlife.app.185](https://doi.org/10.25663/brainlife.app.185) | deprecated |
| Tract Analysis Profiles | Vistasoft | [`1.4`](https://github.com/brainlife/app-tractanalysisprofiles/tree/1.4) | [doi:10.25663/bl.app.43](https://doi.org/10.25663/bl.app.43) | deprecated |

- **Tract Analysis Profiles** (`brainlife.app.361`): tensor maps (FA, MD, RD, AD, and DKI maps if present), with optional NODDI, myelin-map and quantitative MRI inputs.
- **Tract Analysis Profiles (NODDI)** (`brainlife.app.362`): NODDI maps (NDI, ISOVF, ODI).
- **Tract Analysis Profiles - Myelin mapping** (`brainlife.app.523`): a myelin map.
- **Tract Analysis Profiles - Dipy** (`brainlife.app.474`): a Python reimplementation that orients the streamlines of each tract to a QuickBundles centroid and averages the measures with DIPY's Gaussian distance weights; it also needs the reference DWI.

## Running the App

### On brainlife.io

Open the App at [doi.org/10.25663/brainlife.app.523](https://doi.org/10.25663/brainlife.app.523) (or one of the variants above), choose a white matter classification (`neuro/wmc`), the whole-brain tractogram (`neuro/track/tck`) it was derived from and the measure maps, set the number of nodes and press **Execute**. The App archives a `generic/images` dataset of profile plots and a `neuro/tractmeasures` dataset to your project.

Each profile samples the measures at `numnodes` evenly spaced nodes along the tract (200 by default). By default streamlines are weighted by their distance from the tract core (Vistasoft `dtiComputeDiffusionPropertiesAlongFG`, "volume based"). Setting `fiberbased` to `true` uses `Compute_FA_AlongFG` instead, which does not weight the streamlines and is appropriate for testing models involving streamline crossings ("fiber based").

### Locally

The Vistasoft versions are compiled MATLAB code. You need [Singularity/Apptainer](https://apptainer.org/) and [`jq`](https://jqlang.github.io/jq/); `main` runs the MATLAB Compiler Runtime R2017a (`brainlife/mcr:neurodebian1604-r2017a`) and builds the summary CSV with `brainlife/freesurfer-stats:1.2`.

1. Clone this repository and check out the version you want (for example `git checkout 1.13`).
2. Create `config.json` with the paths to your inputs. For tensor measures:

   ```json
   {
       "fa": "/input/tensor/fa.nii.gz",
       "md": "/input/tensor/md.nii.gz",
       "rd": "/input/tensor/rd.nii.gz",
       "ad": "/input/tensor/ad.nii.gz",
       "afq": "/input/wmc/classification.mat",
       "tck": "/input/tck/track.tck",
       "numnodes": 200,
       "fiberbased": false
   }
   ```

   For NODDI measures, replace the tensor keys with `"ndi"`, `"isovf"` and `"odi"`; for a myelin map, use `"myelin"`.

3. Run the App:

   ```bash
   ./main
   ```

The DIPY version (`git checkout dipy-1.0`) runs in `brainlife/pyafq:dev` and reads `dwi`, `classification`, `track`, the tensor keys and `num_nodes` from `config.json`.

To run the uncompiled MATLAB code you need [Vistasoft](https://github.com/vistalab/vistasoft/), [SPM 8](https://www.fil.ion.ucl.ac.uk/spm/software/spm8/) and [jsonlab](https://github.com/fangq/jsonlab).

### Sample data

Download sample datasets from brainlife.io with the [brainlife CLI](https://github.com/brainlife/cli):

```bash
npm install -g brainlife
bl login
mkdir input
bl dataset download 5b96bdce059cf900271924fb && mv 5b96bdce059cf900271924fb input/tensor   # tensor
bl dataset download 5b96bdba059cf900271924fa && mv 5b96bdba059cf900271924fa input/noddi    # NODDI
bl dataset download 5b96bdd2059cf900271924fc && mv 5b96bdd2059cf900271924fc input/wmc
```

## Outputs

```text
images/          PNG and EPS plots of the profile of every measure for every tract
profiles/        <tract>_profiles.csv: per-node mean and standard deviation of every measure
tractmeasures/   tractmeasures.csv: all tracts, nodes and measures in one table
```

In the per-tract CSVs each measure has two columns, `<measure>_1` (mean) and `<measure>_2` (standard deviation). The DIPY version writes the same three directories.

## Troubleshooting

- **`not all output generated`** (Vistasoft versions): the number of profile images differs from the number expected. Check the log for tracts with too few streamlines or for missing measure files.
- **`something failed. checked logs`** (DIPY version): `tractmeasures/tractmeasures.csv` was not created. Check the log of `tractAnalysisProfilesDipy.py` for the failing tract or measure.
- The tractogram must be the whole-brain tractogram from which the white matter classification was made; the classification indexes its streamlines.

## Citing

> [!TIP]
> If you use this App in your work, please cite the App, the software it builds on and brainlife.io (below). GitHub's **Cite this repository** button gives the App citation in APA and BibTeX.

- Yeatman JD, Dougherty RF, Myall NJ, Wandell BA, Feldman HM. Tract profiles of white matter properties: automating fiber-tract quantification. *PLoS ONE* 7(11):e49790 (2012). [doi:10.1371/journal.pone.0049790](https://doi.org/10.1371/journal.pone.0049790)
- Yeatman JD, Wandell BA, Mezer A. Lifespan maturation and degeneration of human brain white matter. *Nature Communications* (2014).
- Yeatman JD, Weiner KS, Pestilli F, Rokem A, Mezer A, Wandell BA. The vertical occipital fasciculus: a century of controversy resolved by in vivo measurements. *PNAS* (2014).
- Tract Analysis Profiles (NODDI) - Deprecated. brainlife.io App. [doi:10.25663/brainlife.app.345](https://doi.org/10.25663/brainlife.app.345)
- Tract Analysis Profiles - Deprecated. brainlife.io App. [doi:10.25663/brainlife.app.185](https://doi.org/10.25663/brainlife.app.185)
- Garyfallidis et al. Dipy, a library for the analysis of diffusion MRI data. *Frontiers in Neuroinformatics* 8 (2014). [doi:10.3389/fninf.2014.00008](https://doi.org/10.3389/fninf.2014.00008)
- Zhang et al. NODDI: Practical in vivo neurite orientation dispersion and density imaging of the human brain. *NeuroImage* 61:1000-1016 (2012). [doi:10.1016/j.neuroimage.2012.03.072](https://doi.org/10.1016/j.neuroimage.2012.03.072)
- Kruper et al. Evaluating the Reliability of Human Brain White Matter Tractometry. *Aperture Neuro* 2021 (2021). [doi:10.52294/e6198273-b8e3-4b63-babb-6e6b0da10669](https://doi.org/10.52294/e6198273-b8e3-4b63-babb-6e6b0da10669)
- Hayashi S, Caron BA, Heinsfeld AS, et al. brainlife.io: a decentralized and open-source cloud platform to support neuroscience research. *Nature Methods* 21:809–813 (2024). [doi:10.1038/s41592-024-02237-2](https://doi.org/10.1038/s41592-024-02237-2)
- Avesani P, McPherson B, Hayashi S, et al. The open diffusion data derivatives, brain data upcycling via integrated publishing of derivatives and reproducible open cloud services. *Scientific Data* 6:69 (2019). [doi:10.1038/s41597-019-0073-y](https://doi.org/10.1038/s41597-019-0073-y)

## Authors

- Brad Caron
- Lindsey Kitchell

**Contributors:** Soichi Hayashi, Franco Pestilli

## Funding

[![Wellcome-226486/Z/22/Z](https://img.shields.io/badge/Wellcome-226486%2FZ%2F22%2FZ-purple.svg)](https://wellcome.org/)
[![NIH-NINDS-U24NS140384](https://img.shields.io/badge/NIH_NINDS-U24NS140384-green.svg)](https://reporter.nih.gov/project-details/11261037)
[![NIH-NINDS-UM1NS132207](https://img.shields.io/badge/NIH_NINDS-UM1NS132207-green.svg)](https://reporter.nih.gov/project-details/11411328)
[![NIH-NIBIB-R01EB029272](https://img.shields.io/badge/NIH_NIBIB-R01EB029272-green.svg)](https://reporter.nih.gov/project-details/10262925)
[![NIH-NIBIB-R01EB030896](https://img.shields.io/badge/NIH_NIBIB-R01EB030896-green.svg)](https://reporter.nih.gov/project-details/10428625)
[![NIH-NIMH-R01MH126699](https://img.shields.io/badge/NIH_NIMH-R01MH126699-green.svg)](https://reporter.nih.gov/project-details/10460628)
[![NIH-NIDCD-R01DC020088](https://img.shields.io/badge/NIH_NIDCD-R01DC020088-green.svg)](https://reporter.nih.gov/project-details/11140284)
[![NSF-BCS-1734853](https://img.shields.io/badge/NSF_BCS-1734853-blue.svg)](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1734853)
[![NSF-BCS-2203524](https://img.shields.io/badge/NSF_BCS-2203524-blue.svg)](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2203524)
[![NSF-IIS-1636893](https://img.shields.io/badge/NSF_IIS-1636893-blue.svg)](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1636893)
[![NSF-OAC-2148729](https://img.shields.io/badge/NSF_OAC-2148729-blue.svg)](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2148729)
[![NSF-IIS-1912270](https://img.shields.io/badge/NSF_IIS-1912270-blue.svg)](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1912270)
[![NSF-IIS-2148700](https://img.shields.io/badge/NSF_IIS-2148700-blue.svg)](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2148700)
[![NSF-OAC-1916518](https://img.shields.io/badge/NSF_OAC-1916518-blue.svg)](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1916518)

<sub>brainlife.io is publicly funded, and acknowledging that support helps keep it running. This App was developed with support from Wellcome, the U.S. National Institutes of Health (including the BRAIN Initiative BRAIN CONNECTS program: The Axonal Projectome EXchange, U24NS140384, and the Center for Mesoscale Connectomics, UM1NS132207) and the U.S. National Science Foundation.</sub>

## License

[MIT](LICENSE).
