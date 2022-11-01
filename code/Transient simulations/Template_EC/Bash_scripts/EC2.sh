FW='EC';
DATADIR="ZZZZZZZZZZZZZZZZZZZ";
declare -a tempcomb=(1.1.5 2.0.75 2.1.5 5.0.75 5.1.5 10.0.75 10.1.5);

for ((i=0; i<7; i++))
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