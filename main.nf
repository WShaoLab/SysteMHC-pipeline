#!/usr/bin/env nextflow

if(params.help) {
    log.info ""
    log.info "SysteMHC computational pipeline"
    log.info "-----------------------------------------"
    log.info "Options:"
    log.info "*  --help:          show this message and exit"
    log.info "*  --MHClass:       of which MHC class the analysis is to performe, MHC1 or MHC2. (default: $params.MHClass)"
    log.info "*  --dda_folder:    folder with DDA files to be searched. eg. --dda_folder '/path/to/DDA_file_folder' "
    log.info "*  --comet_params:  comet parameter file. eg. --comet_params '/path/to/Params/comet.params' "
    log.info "*  --fragger_params:  msfragger parameter file. eg. --fragger_params '/path/to/Params/fragger.params'"
    log.info "*  --msgf_params:  msgf parameter file. eg. --msgf_params '/path/to/Params/msgf.params'"
    log.info "*  --mods:  modification file (default: $params.mods)"
    log.info "*  --comet_threads: number of cores to be used in comet search (default: $params.comet_threads)"
    log.info "*  --fragger_threads: number of cores to be used in msfragger search (default: $params.fragger_threads)"
    log.info "*  --msgf_threads: number of cores to be used in msgf search (default: $params.msgf_threads)"
    log.info "*  --xinteract_threads: number of cores to be used in xinteract (default: $params.xinteract_threads)"
    log.info "*  --iprophet_threads: number of cores to be used in iprophet (default: $params.iprophet_threads)"
    log.info "*  --tpp:           options to pass to TPP xinteract (default: $params.tpp)"
    log.info "*  --decoy:         decoy prefix used in protein_db (default: $params.decoy)"
    log.info "*  --alleles:         alleles used for netMHCpan to predict the binders (default: $params.alleles)"
    log.info "*  --ionstype:         ions fragmentation type, such as HCD, CID, CID-QTOF and ETD.  (default: $params.ionstype)"
    log.info "*  --instrument:       The resolution of MS2 spectra, High or Low.   (default: $params.ionstype)"
    log.info "*  --protein_db:    fasta formatted sequence database file (default: $params.protein_db)"
    log.info "*  --neo:       is the fasta formated sequence file a customized one? Yes or no. (default: $params.ionstype)"
    log.info ""
    log.info "Comet Results will be stored in Results/Comet"
    log.info "Msfragger Results will be stored in Results/Fragger"
    log.info ""
    exit 1
}


process CometSearch {

    cache "lenient"
    cpus params.comet_threads
    //memory 30.GB

    publishDir 'Results/Comet', mode: 'link'

    input:
    file mzML from Channel.fromPath("${params.dda_folder}/*.mzML").collect()
    file comet_params from file(params.comet_params)
    file protein_db from file(params.protein_db)

    output:
    set val ('comet'), file ('*.pep.xml') into cometsearchOut
    file mzML into cometMzMLOut

    script:
    """
    # Set proteins DB
    sed -i 's,database_name.*,database_name = $protein_db,' $comet_params
    sed -i 's,num_threads.*,num_threads = $params.comet_threads,' $comet_params

    comet -P$comet_params $mzML
    sed -ri 's|/tmp/nxf.{11}|${workflow.launchDir}/Results/Comet/|'g *.pep.xml
    sed -ri 's|<search_database local_path="|<search_database local_path="${workflow.launchDir}/Results/Comet/|'g *.pep.xml
    """
}
cometMzMLOut.into{ cometMzMLOut1; cometMzMLOut2; cometMzMLOut3; cometMzMLOut4; cometMzMLOut5; cometMzMLOut6;cometMzMLOut7 }

process PeptideprophetComet {

    cache "lenient"
    cpus "$params.xinteract_threads"
    publishDir 'Results/Comet', mode: 'link'
    
    input:
    file pepxmls from cometsearchOut.map{ it[1] }.collect()
    file protein_db from file(params.protein_db)
    file mzML from Channel.fromPath("${params.dda_folder}/*.mzML").collect()
	    
    output:
    // Normal run
    set val ('merged'), file('comet_merged.pep.xml') into tppPepOutRaw
    file 'comet_merged.pep-MODELS.html'
    file 'comet_merged.pep.xml.index'
    file(protein_db) 
    file '*.log' into PeptideprophetCometlogOut
    
    script:
    """

    if [ "$params.instrument" = "Low" ]; then
        xinteract $params.tppl -p0 -THREADS=$params.xinteract_threads -d$params.decoy -D$params.protein_db -Ncomet_merged.pep.xml $pepxmls > xinteract-comet.log 2>&1
    else
        xinteract $params.tpp -p0 -THREADS=$params.xinteract_threads -d$params.decoy -D$params.protein_db -Ncomet_merged.pep.xml $pepxmls > xinteract-comet.log 2>&1

    fi
    """
}
  
process MSGFPlusSearch {
   
    cache "lenient"
    cpus params.msgf_threads
    // memory 64.GB
    publishDir 'Results/MSGF', mode: 'link'
    
    input:
    file mzMLs from Channel.fromPath("${params.dda_folder}/*.mzML")
    file protein_db from file(params.protein_db)
    file msgf_params from file(params.msgf_params)

    output:
    file ('./MZid')  into msgfSearchOut
    file mzMLs into msgfSearchmzMLOut

    script:
    """
    sed -i 's,DatabaseFile.*,DatabaseFile = $protein_db,' $msgf_params
    sed -i 's,NumThreads.*,NumThreads=$params.msgf_threads,' $msgf_params

    mkdir MZid
    for file in $mzMLs
    do
        java -jar -Xmx64G /opt/MSGFPlus_v20210322/MSGFPlus.jar -s \$file -conf $msgf_params > msgf.log
        mv *.mzid ./MZid
    done
  
    """
}

process MSGFmzid { 

    cache "lenient"
    publishDir 'Results/MSGF', mode: 'link'
	    
    input:
    file(mzidfolder) from msgfSearchOut
	    
    output:
    set val ('modified'), file('*.pepXML') into msgfpepXMLOut
	    
    script:
    """
    for file in ${mzidfolder}/*.mzid
    do
        fraction=`basename \${file} .mzid`
        cat \$file | sed 's/<Enzyme semiSpecific=.*/<Enzyme semiSpecific="true" missdCleavages="-1" id="Tryp">/g' | \
        sed 's/<cvParam.*name="unspecific cleavage"/<cvParam cvRef="PSI-MS" accession="MS:1001251" name="Trypsin"/g' > m.\$fraction.mzid
	idconvert m.\$fraction.mzid --pepXML;
        rm m.\$fraction.mzid
    done
    """
}

process PeptideprophetMSGF {

    cache "lenient"
    cpus "$params.xinteract_threads"
    publishDir 'Results/MSGF', mode: 'link'
	    
    input:
    file pepxmls from msgfpepXMLOut.map{ it[1] }.collect()
    file protein_db from file(params.protein_db)
    file mzML from Channel.fromPath("${params.dda_folder}/*.mzML").collect()
	    
    output:
    // Normal run
    set val ('merged'), file('msgf_merged.pep.xml') into msgfPepOut
    file 'msgf_merged.pep-MODELS.html'
    file 'msgf_merged.pep.xml.index'
    file(protein_db) into msgfdb
    file '*.log' into PeptideprophetMSGFlogOut
    
    script:
    """
    if [ "$params.instrument" = "Low" ]; then
        xinteract $params.tppl -p0 -THREADS=$params.xinteract_threads -d$params.decoy -D$params.protein_db -Nmsgf_merged.pep.xml $pepxmls > xinteract-msgf.log 2>&1
    else
        xinteract $params.tpp -p0 -THREADS=$params.xinteract_threads -d$params.decoy -D$params.protein_db -Nmsgf_merged.pep.xml $pepxmls > xinteract-msgf.log 2>&1
    fi
    """
}

process MsfraggerSearch {
   
    cache "lenient"
    cpus params.fragger_threads
    // memory 214.GB
    publishDir 'Results/Fragger', mode: 'link'
    
    input:

    file protein_db from file(params.protein_db)
    file mzML from cometMzMLOut1.collect()
    file fragger_params from file(params.fragger_params)

    output:
    set val ('fragger'), file ('*.pepXML')  into msfraggerSearchOut
    file mzML into msfraggerSearchmzMLOut
    
    script:
    """
    sed -i 's,num_threads.*,num_threads = $params.fragger_threads,' "$fragger_params"
    sed -i 's,database_name.*,database_name = $protein_db,' "$fragger_params"
    java -jar -Xmx214G /opt/MSFragger-3.4/MSFragger-3.4.jar $fragger_params $mzML

    """
}

process PeptideprophetFragger {
    
    cache "lenient"
    publishDir 'Results/Fragger', mode: 'link'
    
    input:
    file pepxml from msfraggerSearchOut.map{ it[1] }.collect()
    file protein_db from file(params.protein_db)
    file mzML from cometMzMLOut2.collect()

    output:
    set val ('merged'), file ('*pep.xml') into peptideProphetOut
    file '*.log' into peptideProphetlogOut

    script:
    """
    /home/biodocker/bin/philosopher workspace --init --nocheck

    if [ "$params.mode" = "closed" ]; then
        /home/biodocker/bin/philosopher peptideprophet --nonparam --expectscore --decoyprobs --nontt --decoy $params.decoy --database $params.protein_db --combine --output frag_merged ${pepxml} > peptideprophet-fragger.log 2>&1
    
    else
        /home/biodocker/bin/philosopher peptideprophet --nonparam --expectscore --decoyprobs --clevel -2 --nontt --decoy $params.decoy --database $params.protein_db --combine --output frag_merged ${pepxml} > peptideprophet-fragger.log 2>&1
    fi

    """
}

 peptideProphetOut.into{ peptideProphetOut1; peptideProphetOut2; peptideProphetOut3; peptideProphetOut4 }


process Interprophet {
    
    cache "lenient"
    cpus params.iprophet_threads

    publishDir 'Results/SampleLib', mode: 'link'

    input:
    tppPepOutRaw.join(msgfPepOut).join(peptideProphetOut4).set{joined_pepxml}
    set val ('merged'), file(comet_pepxml), file(msgf_pepxml), file(frag_pepxml) from joined_pepxml
    file protein_db from file(params.protein_db)
    file mzML from cometMzMLOut5.collect()

    output:
    file 'iprophet.pep.xml' into iProphetOut
    file mzML into iProphetmzMLOut
    file '*.log' into iprophetlogOut

    script:
    """
    InterProphetParser THREADS=$params.iprophet_threads DECOY=$params.decoy ${comet_pepxml} ${msgf_pepxml} ${frag_pepxml}  iprophet.pep.xml > iprophet.log 2>&1
    
    sed -ri 's|/tmp/nxf.{11}|${workflow.launchDir}/Results/SampleLib/|'g iprophet.pep.xml
    """
}

iProphetOut.into{ iProphetOut1; iProphetOut2 }

process FDRcal {
    
    cache "lenient"

    publishDir 'Results/Peplevel_FDR', mode: 'link'

    input:
    file pepxml from iProphetOut1.collect()

    output:
    file '*result*.csv' into fdrresultOutcsv
    file 'ionslist.csv' into fdrionsOutcsv
    file 'peplist.csv' into fdrpepOutcsv
    file 'iprophet.pep.csv' into iprophetOutcsv
    file 'prob.txt' into FDRout
    file 'probability_result.csv' into fdrprobaOutcsv
    file 'peptide_result.csv' into fdrpepresultOutcsv
    file 'psm_result.csv' into fdrpsmresultOutcsv

    script:
    """
    cat $pepxml | grep '<interprophet_result' | sed 's/.*probability=//g' | awk '{print \$1}' | sed 's/"//g' > iprobability.csv
    /usr/local/bin/pepxml2csv_modified.py $pepxml prophet.pep.csv

    if [ "$params.MHClass" = "MHC1" ]; then
        python3 $baseDir/bin/get_peptide_iprob.py prophet.pep.csv iprobability.csv $params.decoy > prob.txt
    
    else
        python3 $baseDir/bin/get_peptide_iprob_II.py prophet.pep.csv iprobability.csv $params.decoy > prob.txt
    fi
    """
}

process LibraryGen1 {
   
    cache "lenient"
    publishDir 'Results/SampleLib', mode: 'link'

    input:
    file prob from FDRout
    file pepxml from iProphetOut2.collect()
    file protein_db from file(params.protein_db)
    file mzML from cometMzMLOut6.collect()
    file mods from file(params.mods)

    output:
    file 'raw.splib' into LibraryrawlibOut
    file 'cons.splib' into LibraryconslibOut
    file 'raw.sptxt' into LibraryrawtxtOut
    file 'cons.sptxt' into LibraryconstxtOut
    file 'raw.pepidx' into LibraryrawpidxOut
    file 'cons.pepidx' into LibraryconspidxOut
    file '*.spidx' into LibraryidxOut
    file '*.log' into LibrarylogOut
    file 'peplist' into LibrarypeplistOut
    file 'peptides' into LibrarypeptidesOut
    file 'ions' into LibraryionsOut

    script:
    """
    spectrast -cNraw -cP\$(cat ${prob}) -cf "Protein! ~ $params.decoy" -cI$params.ionstype -M$mods ${pepxml}
    spectrast -cNcons -cAC -M$mods raw.splib

    if [ "$params.MHClass" = "MHC1" ]; then
        cat cons.pepidx | grep -v '#' | awk '{if (length(\$1)>=8 && length(\$1)<=14) print \$1}' | sort | uniq > peplist
    
    else
        cat cons.pepidx | grep -v '#' | awk '{if (length(\$1)>=9 && length(\$1)<=25) print \$1}' | sort | uniq > peplist
    fi

    cat cons.pepidx | grep -v '#' | awk '{print \$1}' > peptides
    cat cons.sptxt | grep 'Name' | grep -v 'FullName' | awk '{print \$2}' > ions
    """
}
LibraryrawlibOut.into{ LibraryrawlibOut1; LibraryrawlibOut2}
LibraryrawpidxOut.into{ LibraryrawpidxOut1; LibraryrawpidxOut2; LibraryrawpidxOut3}

if (params.alleles == 'noalleles') {
    
    if (params.neo != 'no') {
    
        process SVariant {
   
            cache "lenient"
            publishDir 'Results/SampleVariant', mode: 'link'

            input:

            file 'peps' from LibrarypeplistOut

            output:
    
            file '*.txt' into sampleVariantOut

            script:
            """
            $baseDir/bin/findneoantigen.sh $peps $baseDir/bin/new_targetDECOY.fasta
            """
        }

    }

    process SPTMslandscape {
   
        cache "lenient"
        publishDir 'Results/SamplePTMs', mode: 'link'

        input:

        file 'pepidx' from LibraryconspidxOut
        file 'pname' from file(params.ptmname)

        output:
    
        file '*.svg' into samplePTMfigOut
        file '*.csv' into samplePTMcsvOut
        file 'modify.csv' into samplePTMmcsvOut
        file '*.svg' into samplePTMpngOut

        script:
        """
        $baseDir/bin/get_modinfo.sh pepidx

        if [ ! -s modify.csv ]; then
            echo "no PTM" > out.csv
            touch none.svg
            touch none.png
        else
            $baseDir/bin/plot_ptmheatmap.py modify.csv $pname
        fi
        
        """
    }
    
    process SSummary {
   
        cache "lenient"
        publishDir 'Results/Summary', mode: 'link'

        input:
        file mzML from Channel.fromPath("${params.dda_folder}/*.mzML").collect()
        file pepidx from LibraryrawpidxOut3.collect()
        file modfile from samplePTMmcsvOut
        file probab from fdrprobaOutcsv

        output:
    
        file 'summary.txt' into SumtxtOut

        script:
        """
        cat $mzML | grep '"ms level" value="2"' | wc | awk '{print \$1}' > tot_spec.txt
        cat $probab | awk -F ',' '{print \$2}' | sed -n '2p' > psm.prob.txt
        cat $probab | awk -F ',' '{print \$2}' | sed -n '3p' > pep.prob.txt
        cat $pepidx | sed -n '3,9p' | sed 's/#//g' > sample.result.txt
        cat $modfile | awk '{print \$5}' | sort | uniq -c | sort -k 1 -g -r | head -5 | sed 's/#//g' > ptm.txt

        echo "### MS2 Spectra: " > summary.txt
        cat tot_spec.txt >> summary.txt
        echo -e "\n### The probability of PSM level 1% FDR:" >> summary.txt
        cat psm.prob.txt >> summary.txt
        echo -e "\n### The probability of peptide level 1% FDR:" >> summary.txt
        cat pep.prob.txt >> summary.txt
        echo -e "\n### Identified results at peptide level 1% FDR:" >> summary.txt
        cat sample.result.txt >> summary.txt
        echo -e "\n### PTM results:" >> summary.txt
        cat ptm.txt >> summary.txt
        """
    }
}

else {

    process NETMHCpan {

        cache "lenient"
        publishDir 'Results/MHCSpecificLib', mode: 'link'

        input:
    
        file peplist from LibrarypeplistOut
        file alleles from file(params.alleles)

        output:
    
        file netmhcout into speNetmhcpanOut
        file 'alles.txt' into spealles
        file 'newname' into spenamealles

        script:
        """
        if [ "$params.MHClass" = "MHC1" ]; then
            echo $alleles > alle
            sed -i 's/_/:/g' alle
            awk -F "," '{for(i=1;i<=NF;i++) a[i,NR]=\$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf a[i,j] " ";print ""}}' alle > alles.txt
    
            for name in \$(cat alles.txt)
            do
	        echo \$name > newname
                sed -i s,':','_', newname 
                /data/netMHCpan-4.1/netMHCpan -p $peplist -BA -xls -a \$name -xlsfile \$(cat newname).xls > \$(cat newname).log
     
            done
            mkdir netmhcout
            mv *.xls netmhcout
            mv *.log netmhcout
    
        else
            echo $alleles > alle
        
            awk -F "," '{for(i=1;i<=NF;i++) a[i,NR]=\$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf a[i,j] " ";print ""}}' alle > alles.txt
    
            for name in \$(cat alles.txt)
            do
	        echo \$name > newname 
                /data/netMHCIIpan-4.1/netMHCIIpan -f $peplist -inptype 1 -BA -xls -a \$name -xlsfile \$(cat newname).xls > \$(cat newname).log
     
            done
            mkdir netmhcout
            mv *.xls netmhcout
            mv *.log netmhcout
        fi
        
            """
    }
    speNetmhcpanOut.into{ speNetmhcpanOut1; speNetmhcpanOut2}
    spealles.into{ spealles1; spealles2}


    process specificpep {
   
        cache "lenient"
        publishDir 'Results/MHCSpecificLib', mode: 'link'

        input:
        file 'alles.txt' from spealles1.collect()
        file spefolder from speNetmhcpanOut1.collect()
        file pepcsv from fdrpepresultOutcsv

        output:
    
        file '*complete.csv' into specificcsvOut
        file '*annotation.csv' into annotationcsvOut
        file '*heatmap.svg' into heatmapOut

        script:
        """
        $baseDir/bin/getspecificpep.py ${spefolder}/ $pepcsv

        """
    }

    process LibraryGen2 {
   
        cache "lenient"
        publishDir 'Results/MHCSpecificLib', mode: 'link'

        input:
        file 'alles.txt' from spealles2.collect()
        file finaltsv from specificcsvOut
        file psmcsv from fdrpsmresultOutcsv

        file mods from file(params.mods)
        file splib from LibraryrawlibOut1.collect()
        file pepidx from LibraryrawpidxOut1.collect()

        output:
    
        file '*.csv' into Library2csvOut
        file '*.log' into Library2logOut
        file 'allbindions.txt' into Library2ionstxtOut
        file 'allbindpeps.txt' into Library2pepstxtOut
        file HLA into Library2Out
        file binderpeps into Library2Outpepfolder

        script:
        """
        mkdir HLA
        for name in \$(cat alles.txt)
        do
            echo \$name > newname
            Parse_binderions.py $psmcsv $finaltsv \$(cat newname)

            spectrast -cNraw_\${name} -M$mods -cI$params.ionstype -cT\${name}_bindions.csv $splib
            spectrast -cNcons_\${name} -M$mods -cAC raw_\${name}.splib
 
            mkdir \$name
            mv *_\$name* \$name
            mv \$name HLA
        done

        cat *bindions.csv | sort | uniq > allbindions.txt 
        cat *bindpeps.csv | sort | uniq > allbindpeps.txt 
    
        mkdir binderpeps
        mv *bindpeps.csv binderpeps
        """
    }

    process LibraryGen3 {
   
        cache "deep"
        publishDir 'Results/MHCLib', mode: 'link'

        input:

        file splib from LibraryrawlibOut2.collect()
        file pepidx from LibraryrawpidxOut2.collect()
        file mods from file(params.mods)
        file bindions from Library2ionstxtOut

        output:
    
        file '*.splib' into Library3libOut
        file 'con*.sptxt' into Library3contxtOut
        file 'con*.pepidx' into Library3conpepOut
        file 'raw*.sptxt' into Library3rawtxtOut
        file 'raw*.pepidx' into Library3rawpepOut
        file '*.spidx' into Library3idxOut
        file '*.log' into Library3logOut

        script:
        """
        spectrast -cNraw_binders -M$mods -cI$params.ionstype -cT$bindions $splib
        spectrast -cNcons_binders -cAC -M$mods raw_binders.splib
        """
    }
    
    process GibbsCluster {
   
        cache "lenient"
        publishDir 'Results/MHCmotifs', mode: 'link'

        input:

        file pepsfolder from Library2Outpepfolder

        output:
    
        file HLAspecificmotifs into LogospeOut
    
        script:
        """

        if [ "$params.MHClass" = "MHC1" ]; then
            mkdir HLAspecificmotifs
            for file in ${pepsfolder}/*.csv
            do
                echo \$file > name
                sed -i s,'_bindpeps.csv','', name 
                sed -i s,'binderpeps/','', name
                $baseDir/bin/get_motif.sh \$file \$(cat name)
             
                mkdir \$(cat name)
                mv \$(cat name).pep* \$(cat name)
                mv gibbs_* \$(cat name)
                mv \$(cat name) HLAspecificmotifs
            done

        else
            mkdir HLAspecificmotifs
            for file in ${pepsfolder}/*.csv
            do
                echo \$file > name
                sed -i s,'_bindpeps.csv','', name 
                sed -i s,'binderpeps/','', name
                $baseDir/bin/get_motif_II.sh \$file \$(cat name)
             
                mkdir \$(cat name)
                mv \$(cat name).pep* \$(cat name)
                mv gibbs_* \$(cat name)
                mv \$(cat name) HLAspecificmotifs
            done
        fi
        """
    }

    process SamPTMslandscape {
   
        cache "lenient"
        publishDir 'Results/SamplePTMs', mode: 'link'

        input:

        file 'pepidx' from LibraryconspidxOut
        file 'pname' from file(params.ptmname)

        output:
    
        file '*.svg' into samplePTMfigOut
        file '*.csv' into samplePTMcsvOut
        file 'modify.csv' into samplePTMmcsvOut
        file '*.png' into samplePTMpngOut

        script:
        """
        $baseDir/bin/get_modinfo.sh pepidx

        if [ ! -s modify.csv ]; then
            echo "no PTM" > out.csv
            touch none.svg
            touch none.png
        else
            $baseDir/bin/plot_ptmheatmap.py modify.csv $pname
        fi
        
        """
    }

    process PTMslandscape {
   
        cache "lenient"
        publishDir 'Results/MHCPTMs', mode: 'link'

        input:

        file 'pepidx' from Library3conpepOut
        file 'pname' from file(params.ptmname)

        output:
    
        file '*.svg' into PTMfigOut
        file '*.csv' into PTMcsvOut
        file 'modify.csv' into PTMmcsvOut
        file '*.png' into PTMpngOut
        script:
        """
        $baseDir/bin/get_modinfo.sh pepidx

        if [ ! -s modify.csv ]; then
            echo "no PTM" > out.csv
            touch none.svg
            touch none.png
        else
            $baseDir/bin/plot_ptmheatmap.py modify.csv $pname
        fi
       
        """
    }

    if (params.neo == 'y') {
    
        process SamVariant {
   
            cache "lenient"
            publishDir 'Results/SampleVariant', mode: 'link'

            input:

            file 'peps' from LibrarypeplistOut

            output:
    
            file '*.txt' into SamVariantOut

            script:
            """
            $baseDir/bin/findneoantigen.sh $peps $baseDir/bin/new_targetDECOY.fasta
            """
        }

        process Variant {
   
            cache "lenient"
            publishDir 'Results/Variant', mode: 'link'

            input:

            file 'peps' from Library2pepstxtOut

            output:
    
            file '*.txt' into VariantOut

            script:
            """
            $baseDir/bin/findneoantigen.sh $peps $baseDir/bin/new_targetDECOY.fasta
            """
        }

    }

    process Summary {
   
        cache "lenient"
        publishDir 'Results/Summary', mode: 'link'

        input:
        file mzML from Channel.fromPath("${params.dda_folder}/*.mzML").collect()
        file pepidx from LibraryrawpidxOut3.collect()
        file mhcpep from Library3rawpepOut
        file modfile from PTMmcsvOut
        file probab from fdrprobaOutcsv

        output:
    
        file 'summary.txt' into SumtxtOut

        script:
        """
        cat $mzML | grep '"ms level" value="2"' | wc | awk '{print \$1}' > tot_spec.txt
        cat $probab | awk -F ',' '{print \$2}' | sed -n '2p' > psm.prob.txt
        cat $probab | awk -F ',' '{print \$2}' | sed -n '3p' > pep.prob.txt
        cat $pepidx | sed -n '3,9p' | sed 's/#//g' > sample.result.txt
        cat $mhcpep | sed -n '3,9p' | sed 's/#//g' > mhc.result.txt 
        cat $modfile | awk '{print \$5}' | sort | uniq -c | sort -k 1 -g -r | head -5 | sed 's/#//g' > ptm.txt

        echo "### MS2 Spectra: " > summary.txt
        cat tot_spec.txt >> summary.txt
        echo -e "\n### The probability of PSM level 1% FDR:" >> summary.txt
        cat psm.prob.txt >> summary.txt
        echo -e "\n### The probability of peptide level 1% FDR:" >> summary.txt
        cat pep.prob.txt >> summary.txt
        echo -e "\n### Identified results at peptide level 1% FDR:" >> summary.txt
        cat sample.result.txt >> summary.txt
        echo -e "\n### MHC binder results:" >> summary.txt
        cat mhc.result.txt >> summary.txt
        echo -e "\n### PTM results:" >> summary.txt
        cat ptm.txt >> summary.txt
        """
    }

}

workflow.onComplete {
    // Make the results folder writable for the www-data group
    "chmod g+w Results/Comet".execute()
    "chmod g+w Results/Fragger".execute()
    "chmod g+w Results/MSGF".execute()
    "chmod g+w Results/SampleLib".execute()
    "chmod g+w Results/Peplevel_FDR".execute()
    "chmod g+w Results/MHCLib".execute()
    "chmod g+w Results/MHCSpecificLib".execute()
    "chmod g+w Results/MHCmotifs".execute()
    "chmod g+w Results/MHCPTMs".execute()
    "chmod g+w Results/Variant".execute()
    "chmod g+w Results/Summary".execute()
}
