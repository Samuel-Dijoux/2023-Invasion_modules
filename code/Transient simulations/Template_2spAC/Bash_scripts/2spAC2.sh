FW='2spAC';
DATADIR=" ZZZZZ ";
declare -a tempcomb=(1.5 3 3.75 7.5 15);

for ((i=0; i<4; i++))
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