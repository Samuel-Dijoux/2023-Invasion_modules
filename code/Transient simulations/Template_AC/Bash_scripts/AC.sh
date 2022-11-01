FW='AC';
DATADIR="ZZZZZZZZZ";

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

declare -a tempcomb=(0.1.10 0.2.5 0.2.10 0.5.2 0.5.5 0.5.10)
for ((i=0; i<6; i++))
do
    FIL1=$"BMRs"_"${tempcomb[$i]}";
    DIR="$DATADIR/$FIL1";
    cd $DIR;
        
    chmod 744 *;
    for j in {1..40}
    do
        NUM=$((j));
        COMB="$FW"_"${tempcomb[$i]}"_"$NUM";
        SH="$COMB.sh";
        # Convert the file in the right format
        dos2unix $SH;
        # Submit the jobs
        qsub -m a $SH;
    done
done