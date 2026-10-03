SP=genome_474
GENOME=00_genome/${SP}.fa
OUT=$(pwd)/03_repeat/${SP}.TRASH

mkdir -p $OUT

conda run -n trash_env path/to/home/sofeware/TRASH-main/TRASH_run.sh $GENOME \
  --o $OUT \
  --par 32
