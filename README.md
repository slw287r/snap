# SNAP

Scalable Nucleotide Alignment Program - <https://www.microsoft.com/en-us/research/project/snap/>

## Overview

SNAP is a fast and accurate aligner for short DNA reads. It is optimized for
modern read lengths of 100 bases or higher, and takes advantage of these reads
to align data quickly through a hash-based indexing scheme.

It also includes support for sorting, marking duplicates and indexing its results, eliminating the 
need for several pipeline stages used by other aligners.

## Binaries

The SNAP executable
- [v2.0.2 For Linux](https://1drv.ms/u/s!AhuEg_0yZD86hcpYCkpLlDktZnVaow?e=QRDhs4)
- [v2.0.2 for Windows 10](https://1drv.ms/u/s!AhuEg_0yZD86hcpZQUgOEMrmA5qaLA?e=eUNeHZ)
- [v2.0.0 for OSX](https://1drv.ms/u/s!AhuEg_0yZD86hcphrIjwoeTjdSvgoA?e=coSU85)

The SNAPCommand tool
- [SNAPCommand for Linux](https://1drv.ms/u/s!AhuEg_0yZD86hcpdvv0ZBdB1BqF57g?e=IHVbq2>)
- [SNAPCommand for Windows 10](https://1drv.ms/u/s!AhuEg_0yZD86hcpaSLKPRGJ6dcvVgA?e=vXH8y6)
- [SNAPCommand for OSX](https://1drv.ms/u/s!AhuEg_0yZD86hcpgy-ONBaw0DjFpTQ?e=cMc6eE)


## Documentation

SNAP has a one page [Quick Start Guide](https://1drv.ms/b/s!AhuEg_0yZD86hcpcvhSwRyDwk1Ru0Q?e=4BvzLn) and a more extensive [Manual](https://1drv.ms/b/s!AhuEg_0yZD86hcpblUt-muHKYsG8fA?e=mbyUP5).

## Building

SNAP runs on Windows, Linux and OSX.

For Windows, we provide a Visual C++ solution, `snap.sln`. Requirements:
- Visual Studio 2022

When you build it, you will have to set it to build for x64, not "Any CPU" or 32 bit.

For Linux, simply type `make`. Requirements:
- g++ version 4.8.5 or later
- zlib 1.2.11 or later from http://zlib.net/

## Notes for this branch (dev)

This branch is tailored for CAPx's host depletion and others statistics include host rRNA, HSK and synthetic IC. It gathers extra statistics information while depleting the host reads to speed up the pipeline.

Coordinates presets are required in the `SNAPLib/T2TrRNA.h` file. The coordinate info can be obtained by simulating SE reads and map them back to the updated reference genome with `snap-alginer single` with the `--ppp` option, which is only available when compiled with the `-D_DEBUG_3P` flag in `Makefile`. After done with the `SNAPLib/T2TrRNA.h` updates, comment out the `CXX = g++ -D_DEBUG_3P` line and make again to get rid of the `printPaddedPosition` check which may slow the mapping step.



<details><summary><h3>Get contig offsets from the SNAP index file</h3></summary>

1. Get the total number of contigs 
```
$ head -n1 Genome | cut -d' ' -f2
74
```
 
2. Get their offsets (1-based coordinates)
```
$ head -n$(head -n1 Genome | cut -d' ' -f2) Genome | sed '1d' | awk '{print $3"\t"$8"\t"$1}'
```


```
0	chr1	2000
1	chr2	248391328
2	chr3	491090080
3	chr4	692198028
4	chr5	885774973
5	chr6	1067822412
6	chr7	1239951040
7	chr8	1400520468
8	chr9	1546781799
9	chr10	1697401046
10	chr11	1832161180
11	chr12	1967290949
12	chr13	2100617497
13	chr14	2214186183
14	chr15	2315349675
15	chr16	2415104870
16	chr17	2511437244
17	chr18	2595716141
18	chr19	2676260679
19	chr20	2737970043
20	chr21	2804182298
21	chr22	2849274980
22	chrX	2900601906
23	chrY	3054863472
24	chrM	3117325501
25	IC1	3117344070
26	IC2	3117346729
27	IC3	3117349449
28	IC4	3117352049
29	IC5	3117354797
30	IC6	3117357481
31	IC7	3117360166
32	IC8	3117362797
33	IC9	3117365479
34	IC10	3117368203
35	IC11	3117370971
36	IC12	3117373675
37	IC13	3117376284
38	IC14	3117378923
39	IC15	3117381590
40	IC16	3117384389
41	IC17	3117387022
42	IC18	3117389664
43	IC19	3117392434
44	IC20	3117395169
45	IC21	3117397924
46	IC22	3117400663
47	IC23	3117403451
48	IC24	3117406122
49	IC25	3117408852
50	IC26	3117411467
51	IC27	3117414211
52	IC28	3117417000
53	IC29	3117419638
54	IC30	3117422299
55	IC31	3117424964
56	IC32	3117427752
57	IC33	3117430365
58	IC34	3117433048
59	IC35	3117435675
60	IC36	3117438471
61	IC37	3117441131
62	IC38	3117443731
63	IC39	3117446368
64	IC40	3117449159
65	IC41	3117451884
66	IC42	3117454576
67	IC43	3117457228
68	IC44	3117459895
69	IC45	3117462659
70	IC46	3117465371
71	IC47	3117467974
72	IC48	3117470635
73	IC49	3117473263
```

This is enough for preparing intervals for internal controls, subtract one to get the start coordinate 3117473262 for IC49 (0-based in `T2TrRNA.h`), add IC length (such as 600) to get the end coordinate 3117473862.

</details>

<details><summary><h3>Get padded intervals for host ribosomal RNA, HSK and IC</h3></summary>

Get padded intervals for host ribosomal RNA genes or the HSK (Ribonuclease P protein subunit p30 and p38)

1. Simulate SE reads from the gene's reference sequences
2. Map the SE reads to snap reference and output padded mapping positions (via `--ppp`, uncomment `#CXX = g++ -D_DEBUG_3P` in Makefile and make again to enable)
3. Combine the padded mapping intervals

* Example
```
snap-aligner single /mnt/ONCOBOX/geneplus/pipeline/CNC/db/pubDB/human/gptk/tngs_max/snap -fastq rrna_SE120.fq -o -sam -ppp > rrna_SE120.sam 2>rrna_SE120.log
```
padded mapping positions go to log file with leading `[DEBUG]` symbol.
</details>