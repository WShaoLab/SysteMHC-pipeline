#!/bin/sh

cat $1 | awk '{if (length($1)==8) print $1}' > $2.pep8.txt;
cat $1 | awk '{if (length($1)==9) print $1}' > $2.pep9.txt;
cat $1 | awk '{if (length($1)==10) print $1}' > $2.pep10.txt;
cat $1 | awk '{if (length($1)==11) print $1}' > $2.pep11.txt;
cat $1 | awk '{if (length($1)==12) print $1}' > $2.pep12.txt;
cat $1 | awk '{if (length($1)==13) print $1}' > $2.pep13.txt;
cat $1 | awk '{if (length($1)==14) print $1}' > $2.pep14.txt;

/gibbsCluster/gibbscluster-2.0/gibbscluster -f $1 -P gibbs_8 -l8 -g1-4 > $2.pep8;
        
/gibbsCluster/gibbscluster-2.0/gibbscluster -f $1 -P gibbs_9 -l9 -g1-4 > $2.pep9;
       
/gibbsCluster/gibbscluster-2.0/gibbscluster -f $1 -P gibbs_10 -l10 -g1-4 > $2.pep10;
              
/gibbsCluster/gibbscluster-2.0/gibbscluster -f $1 -P gibbs_11 -l11 -g1-4 > $2.pep11;
             
/gibbsCluster/gibbscluster-2.0/gibbscluster -f $1 -P gibbs_12 -l12 -g1-4 > $2.pep12;
             
/gibbsCluster/gibbscluster-2.0/gibbscluster -f $1 -P gibbs_13 -l13 -g1-4 > $2.pep13;
    
/gibbsCluster/gibbscluster-2.0/gibbscluster -f $1 -P gibbs_14 -l14 -g1-4 > $2.pep14;



 
