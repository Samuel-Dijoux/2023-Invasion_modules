FW='2spAC';

DATADIR=" ZZZZZ "
for gamma in {1,2,4,5,10,20,25,50,100}
do
    Gamma=$((gamma));
    FIL1=$"BMRs"_"$Gamma";
    DIR="$DATADIR/$FIL1";
    cd $DIR;
        
    chmod 744 *;
    for i in {1..40}
    do
        NUM=$((i));
        COMB="$FW"_"$Gamma"_"$NUM";
        SH="$COMB.sh";
        # Convert the file in the right format
        dos2unix $SH;
        # Submit the jobs
        qsub -m a $SH;    
    done
done

Gamma=2.5;
FIL1=$"BMRs"_"$Gamma";
cd $DIR;
chmod 744 *;
for i in {1..40}
do
    NUM=$((i));
    COMB="$FW"_"$Gamma"_"$NUM";
    SH="$COMB.sh";
    # Convert the file in the right format
    dos2unix $SH;
    # Submit the jobs
    qsub -m a $SH;    
done
