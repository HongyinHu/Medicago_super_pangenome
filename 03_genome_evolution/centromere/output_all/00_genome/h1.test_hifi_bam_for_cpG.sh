samtools view genome_Msa.hifi_reads.bam  | head -n 1000 | \                                                                                                                                                    
	awk '
/MM:Z:/ {mm++}
/ML:B:C/ {ml++}
/fi:B|fp:B|ri:B|rp:B|ip:B|pw:B/ {kin++}
END{
  print "MM tags:", mm+0
    print "ML tags:", ml+0
	  print "kinetics tags:", kin+0
	  
  }'
