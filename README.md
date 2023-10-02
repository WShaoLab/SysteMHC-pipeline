# SysteMHC-pipeline
This is a  computational pipeline for MS-based immunopeptidomics. 

# Installation
1. Install [Docker](https://docs.docker.com/get-docker/)
2. Install [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation)

# Usage
$ nextflow run /SysteMHC-pipeline/ --help  
N E X T F L O W  ~  version 21.10.6  
Launching `/SysteMHC-pipeline/main.nf` [crazy_heisenberg] - revision: f88b0d422b  

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
*  --alleles:         alleles used for netMHCpan to predict the binders (default: noalleles)
*  --ionstype:         ions fragmentation type, such as HCD, CID, CID-QTOF and ETD.  (default: HCD)
*  --instrument:       The resolution of MS2 spectra, High or Low.   (default: High)
*  --protein_db:    fasta formatted sequence database file (default: Params/proteome.fasta)
*  --neo:       is the fasta formated sequence file a customized one? Yes or no. (default: no)
------------------------------------------

The output of this pipeline is a folder named 'Results'  
##########################################  
Comet Results will be stored in Results/Comet  
MSGF+ Results will be stored in Results/MSGF  
Msfragger Results will be stored in Results/Fragger  
FDR calculation Results will be stored in Results/Peplevel_FDR  
Sample Library will be stored in Results/SampleLib  
Sample PTMs Results will be stored in Results/SamplePTMs  
Sample non-canonical peptides will be stored in Results/SampleVariant  
MHC specific library will be stored in Results/MHCSpecificLib  
MHC specific non-canonical peptides will be stored in Results/Variant  
MHC specific motifs will be stored in Results/MHCmotifs  
MHC libraries will be stored in Results/MHCLib  
MHC PTMs will be stored in Results/MHCPTMs  
A brief summary will be stored in Results/Summary  
##########################################  

# Example
nextflow run SysteMHC-pipeline --dda_folder /path/to/dda_folder --comet_params Params/comet-hh-nofixmod.params \  
--fragger_params Params/fragger-classI-hh-offset.params --msgf_params Params/msgf-QE-classI.params \  
--mods Params/mymods --alleles HLA-A02_01,HLA-B35_01 --ionstype HCD --instrument High --MHClass MHC1 \  
--protein_db Params/uniprot-human-reviewed_202108_iRT_plus_cRAP_targDecoy.fasta \   
--neo no -resume  

