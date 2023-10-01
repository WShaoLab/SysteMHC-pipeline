#!/bin/sh


cat $1 | awk '{if (length($1)==9) print $1}' > $2.pep9.txt;
cat $1 | awk '{if (length($1)==10) print $1}' > $2.pep10.txt;
cat $1 | awk '{if (length($1)==11) print $1}' > $2.pep11.txt;
cat $1 | awk '{if (length($1)==12) print $1}' > $2.pep12.txt;
cat $1 | awk '{if (length($1)==13) print $1}' > $2.pep13.txt;
cat $1 | awk '{if (length($1)==14) print $1}' > $2.pep14.txt;
cat $1 | awk '{if (length($1)==15) print $1}' > $2.pep15.txt;
cat $1 | awk '{if (length($1)==16) print $1}' > $2.pep16.txt;
cat $1 | awk '{if (length($1)==17) print $1}' > $2.pep17.txt;
cat $1 | awk '{if (length($1)==18) print $1}' > $2.pep18.txt;
cat $1 | awk '{if (length($1)==19) print $1}' > $2.pep19.txt;
cat $1 | awk '{if (length($1)==20) print $1}' > $2.pep20.txt;
cat $1 | awk '{if (length($1)==21) print $1}' > $2.pep21.txt;
cat $1 | awk '{if (length($1)==22) print $1}' > $2.pep22.txt;
cat $1 | awk '{if (length($1)==23) print $1}' > $2.pep23.txt;
cat $1 | awk '{if (length($1)==24) print $1}' > $2.pep24.txt;
cat $1 | awk '{if (length($1)==25) print $1}' > $2.pep25.txt;
        
/gibbsCluster/gibbscluster-2.0/gibbscluster -f $1 -P gibbs_9 -l9 -g1-4 > $2.pep9;
       
/gibbsCluster/gibbscluster-2.0/gibbscluster -f $1 -P gibbs_10 -l10 -g1-4 > $2.pep10;
              
/gibbsCluster/gibbscluster-2.0/gibbscluster -f $1 -P gibbs_11 -l11 -g1-4 > $2.pep11;
             
/gibbsCluster/gibbscluster-2.0/gibbscluster -f $1 -P gibbs_12 -l12 -g1-4 > $2.pep12;
             
/gibbsCluster/gibbscluster-2.0/gibbscluster -f $1 -P gibbs_13 -l13 -g1-4 > $2.pep13;
    
/gibbsCluster/gibbscluster-2.0/gibbscluster -f $1 -P gibbs_14 -l14 -g1-4 > $2.pep14;

/gibbsCluster/gibbscluster-2.0/gibbscluster -f $1 -P gibbs_15 -l15 -g1-4 > $2.pep15;

/gibbsCluster/gibbscluster-2.0/gibbscluster -f $1 -P gibbs_16 -l16 -g1-4 > $2.pep16;
       
/gibbsCluster/gibbscluster-2.0/gibbscluster -f $1 -P gibbs_17 -l17 -g1-4 > $2.pep17;
              
/gibbsCluster/gibbscluster-2.0/gibbscluster -f $1 -P gibbs_18 -l18 -g1-4 > $2.pep18;
             
/gibbsCluster/gibbscluster-2.0/gibbscluster -f $1 -P gibbs_19 -l19 -g1-4 > $2.pep19;
             
/gibbsCluster/gibbscluster-2.0/gibbscluster -f $1 -P gibbs_20 -l20 -g1-4 > $2.pep20;
    
/gibbsCluster/gibbscluster-2.0/gibbscluster -f $1 -P gibbs_21 -l21 -g1-4 > $2.pep21;

/gibbsCluster/gibbscluster-2.0/gibbscluster -f $1 -P gibbs_22 -l22 -g1-4 > $2.pep22;

/gibbsCluster/gibbscluster-2.0/gibbscluster -f $1 -P gibbs_23 -l23 -g1-4 > $2.pep23;
    
/gibbsCluster/gibbscluster-2.0/gibbscluster -f $1 -P gibbs_24 -l24 -g1-4 > $2.pep24;

/gibbsCluster/gibbscluster-2.0/gibbscluster -f $1 -P gibbs_25 -l25 -g1-4 > $2.pep25;



 
