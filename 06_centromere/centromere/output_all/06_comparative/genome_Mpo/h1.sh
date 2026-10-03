CHR=Chr1
CHR_SIZE=64347262

awk -v c=${CHR} -v s=${CHR_SIZE} '
$1==c {sum += $3-$2}
END{
    print "Chr peak length =", sum;
	    print "Chr peak ratio =", sum/s;
		
}' genome_Mpo.CENH3.q20_peaks.broadPeak.bed
