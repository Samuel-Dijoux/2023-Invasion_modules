FW='IGPC';

DATADIR="ZZZZZZZZZZ"

for alpha in {1,2,5,10}
do
    for beta in {1,2,5,10}
    do
        Alpha=$((alpha));
        Beta=$((beta));
        PPMR="$Alpha.$Beta";
        FIL1=$"BMRs"_"$PPMR";
        DIR="$DATADIR/$FIL1";
        cd $DIR;
        
        chmod 744 *;
        for i in {1..40}
        do
            NUM=$((i));
            COMB="$FW"_"$PPMR"_"$NUM";
            SH="$COMB.sh";
            # Convert the file in the right format
            dos2unix $SH;
            # Submit the jobs
            qsub -m a $SH;
        done
    
    done
done
