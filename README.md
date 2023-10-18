# SysteMHC-pipeline
This is a computational pipeline for MS-based immunopeptidomics and used in [SysteMHC Atlas 2.0](https://systemhc.sjtu.edu.cn/)

# Overview
This pipeline contains several modules:  
1. Database Search: Three search engines([Comet](https://comet-ms.sourceforge.net/), [MSFragger](https://github.com/Nesvilab/MSFragger) and [MSGF+](https://github.com/MSGFPlus/msgfplus)) are used and the result of the three search engines are combined by iProphet intergrated in [Trans-Proteomics Pipeline(TPP)](http://tools.proteomecenter.org/wiki/index.php?title=Software:TPP).
2. FDR estimation: Global peptide level FDR estimation.
3. Binding affinity prediction: [NetMHCpan](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/) and [NetMHCIIpan](https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.0/) are used for MHC class I and class II, respectively.
4. Motif deconvolution: [GibbsCluster](https://services.healthtech.dtu.dk/services/GibbsCluster-2.0/) is used to perform motif deconvolution.
5. PTM visulization: Heatmap of PTM distribution.
6. Library generation: Sample-specific and allele-specific libraries are generated by SpectraST intergrated in [TPP](http://tools.proteomecenter.org/wiki/index.php?title=Software:TPP). 

# Prerequisite
To run this pipeline, we suppose that the random access memory (RAM) of your computational platform is no less than **210 GB**, due to the memory utilization in [MSFragger](https://github.com/Nesvilab/MSFragger) when performing database search in offset mode. 

# Installation
1. Install [Docker](https://docs.docker.com/get-docker/)
2. Install [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation)
3. Download SysteMHC-pipeline:  
   `git clone https://github.com/WShaoLab/SysteMHC-pipeline`  
   or  
   download the `ZIP` file.
4. Move the SysteMHC-pipeline folder to somewhere you want, such as /www/SysteMHC-pipeline
5. The `bin` file contains the scripts that will be used in the pipeline, the `Params` folder provides several default parameters and the `Fasta` provides Human and Mouse target-Decoy fasta formated sequence. For better and easy use of the pipeline, users are recommended to edit the nextflow.config file.  
   * change the memoey in param `msgf_memory` and `fragger_memory` according to the resource of your computer. This is **required**.
   * change the thread in param `comet_threads`, `msgf_threads`, `fragger_threads`, `xinteract_threads` and `iprophet_threads` to make them compatiable with your computation resource. This is **required**.
   * Also, users can change `mods`,`protein_db`,`comet_params`,`msgf_params`,`fragger_params` accordingly. This is **optional** as these params need to be changed to be compatiable with the mzML files and can be override by command line options.
   * The nextflow.config also defines the containers used in the pipeline. These containers are all created by Docker.
6. In the first time when the user runs the pipeline, the contatiners will be downloaded automatically, so there is no need to download them manually (If something wrong in the downloading process, please download the docker images manually). 

# Usage
$ nextflow run SysteMHC-pipeline/ --help  
N E X T F L O W  ~  version 21.10.6  
Launching `SysteMHC-pipeline/main.nf` [crazy_heisenberg] - revision: f88b0d422b  

-----------------------------------------
Options:
*  --help:          show this message and exit
*  --MHClass:       of which MHC class the analysis is to perform, MHC1 or MHC2. (default: MHC1)
*  --dda_folder:    folder with DDA files to be searched. eg. --dda_folder '/path/to/DDA_file_folder' 
*  --comet_params:  comet parameter file. eg. --comet_params '/path/to/Params/comet.params' 
*  --fragger_params:  msfragger parameter file. eg. --fragger_params '/path/to/Params/fragger.params'
*  --msgf_params:  msgf parameter file. eg. --msgf_params '/path/to/Params/msgf.params'
*  --mods:  modification file (default: Params/mymods)
*  --comet_threads: number of cores to be used in comet search (default: 32)
*  --fragger_threads: number of cores to be used in msfragger search (default: 32)
*  --msgf_threads: number of cores to be used in msgf search (default: 24)
*  --xinteract_threads: number of cores to be used in xinteract (default: 32)
*  --iprophet_threads: number of cores to be used in iprophet (default: 32)
*  --tpp:           options to pass to TPP xinteract (default: -OAPNd -PPM -eN)
*  --decoy:         decoy prefix used in protein_db (default: DECOY_)
*  --fdr:           peptide level fdr (default: 0.01)
*  --alleles:         alleles used for netMHCpan to predict the binders (default: noalleles)
*  --ionstype:         ions fragmentation type, such as HCD, CID, CID-QTOF and ETD.  (default: HCD)
*  --instrument:       The resolution of MS2 spectra, High or Low.   (default: High)
*  --protein_db:    fasta formatted sequence database file (default: Params/proteome.fasta)
*  --neo:       is the fasta formated sequence file a customized one? Yes or no. (default: no)
------------------------------------------

The output of this pipeline is a folder named **Results**  
###
Comet Results will be stored in **Results/Comet**  
MSGF+ Results will be stored in **Results/MSGF**  
Msfragger Results will be stored in **Results/Fragger**  
FDR calculation Results will be stored in **Results/Peplevel_FDR**  
Sample Library will be stored in **Results/SampleLib**  
Sample PTMs Results will be stored in **Results/SamplePTMs**  
Sample non-canonical peptides will be stored in **Results/SampleVariant**  
MHC specific library will be stored in **Results/MHCSpecificLib**  
MHC specific non-canonical peptides will be stored in **Results/Variant**  
MHC specific motifs will be stored in **Results/MHCmotifs**  
MHC libraries will be stored in **Results/MHCLib**  
MHC PTMs will be stored in **Results/MHCPTMs**  
A brief summary will be stored in **Results/Summary**  
###
------------------------------------------

# Example
nextflow run /www/SysteMHC-pipeline/ --dda_folder /path/to/dda_folder/ \  
--comet_params /path/to/comet-hh-nofixmod.params \  
--fragger_params /path/to/fragger-classI-hh-offset.params \  
--msgf_params /path/to/msgf-QE-classI.params \  
--mods Params/mymods --alleles HLA-A02_01,HLA-B35_01 \  
--ionstype HCD --instrument High --MHClass MHC1 \  
--protein_db Params/uniprot-human-reviewed_202108_iRT_plus_cRAP_targDecoy.fasta \   
--fdr 0.01 --decoy DECOY_ --neo no

### if you configured the params in nextflow.config accordingly, then the following command will work the same as above  
nextflow run /www/SysteMHC-pipeline/ --dda_folder /path/to/dda_folder/ \  
--alleles HLA-A02_01,HLA-B35_01  

### Note 
* `--dda_folder` is the folder path of the DDA mzML files. Raw files must be converted into mzML in centroid mode by tools, such as MSconvert. 
* `--decoy` is the flag of reversed sequence, it must be 'DECOY_', 'rev_' or something else according to the fasta formated sequence used in database search. 
* For the params of `--comet_params`, `--fragger_params` and `--msgf_params`, 'hh' represents high resolution of MS1 and MS2, while 'hl' represents high resolution of MS1 and low resolution of MS2. 'classI' is for MHC class I, whereas 'classII' is for MHC class II. Users can refer to [Comet](https://comet-ms.sourceforge.net/), [MSFragger](https://github.com/Nesvilab/MSFragger) and [MSGF+](https://github.com/MSGFPlus/msgfplus) to get the detailed information.
* For `--msgf_params`, 'QE' means 'Q-Exactive', 'Q-Exactive HF' or 'Q-Exactive Plus'; 'Lumos' means 'Orbitrap Fusion Lumos'; 'TOF' means 'Triple-TOF'; 'LTQ' means 'LTQ-FT', and 'LTQ-Orbitrap'. For more detailed information, please refer to [MSGF+](https://github.com/MSGFPlus/msgfplus).
* `--alleles` indicates the alleles used in [NetMHCpan](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/) and [NetMHCIIpan](https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.0/). Here `:` was changed to `_` to be compatiable with the commandline.
* `--ionstype` means the fragmentation method in mzML file.
* If a customized fasta file is used, the `--neo=Yes` means to find peptides not contained in the uniprot-reviewed proteins. 
* In total, The params must be compatiable with each other.
###

# Contact Us
For issues in using **SysteMHC-pipeline**, please report to this GitHub repository.
