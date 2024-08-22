test:
	#snap-aligner single . -fastq a.fq -x -o -sam - -d 15 -D 0 -f -C-- -A- -mrl 35 -t 4 >a.sam 2>a.log
	snap-aligner single . -fastq a.fq -x -o -sam - -d 5 -D 0 -f -F u -C-- -A- -mrl 35 -t 4 >a.sam 2>a.log

clean:
	rm -f a.sam a.log
