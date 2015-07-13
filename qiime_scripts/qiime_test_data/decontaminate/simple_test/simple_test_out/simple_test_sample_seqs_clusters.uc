# uclust --input /var/folders/x7/d7tds9qn0_74nljddc_n03q40000gn/T/uclust_fasta_sortaHyTHN.fasta --id 0.99 --tmpdir /var/folders/x7/d7tds9qn0_74nljddc_n03q40000gn/T --w 8 --stepwords 8 --maxaccepts 1 --libonly --maxrejects 8 --lib ./qiime_scripts/qiime_test_data/decontaminate/simple_test/simple_test_contam_ref_db.fna --uc ./qiime_scripts/qiime_test_data/decontaminate/simple_test/simple_test_out/simple_test_sample_seqs_clusters.uc
# version=1.2.22
# Tab-separated fields:
# 1=Type, 2=ClusterNr, 3=SeqLength or ClusterSize, 4=PctId, 5=Strand, 6=QueryStart, 7=SeedStart, 8=Alignment, 9=QueryLabel, 10=TargetLabel
# Record types (field 1): L=LibSeed, S=NewSeed, H=Hit, R=Reject, D=LibCluster, C=NewCluster, N=NoHit
# For C and D types, PctId is average id with seed.
# QueryStart and SeedStart are zero-based relative to start of sequence.
# If minus strand, SeedStart is relative to reverse-complemented seed.
N	*	400	*	*	*	*	*	contam1 denovo35 SHOG_15	*
L	0	400	*	*	*	*	*	denovo0 SHNZ_90	*
H	0	400	99.2	+	0	0	400M	contam2 denovo39 SHNO_87	denovo0 SHNZ_90
L	3	400	*	*	*	*	*	denovo3 SHNO_86	*
H	3	400	100.0	+	0	0	400M	contam3 denovo40 SHNO_86	denovo3 SHNO_86
L	4	400	*	*	*	*	*	denovo4 SHNO_85	*
H	4	400	100.0	+	0	0	400M	contam4 denovo41 SHNO_85	denovo4 SHNO_85
N	*	400	*	*	*	*	*	otu1 denovo43 SHNX_68	*
N	*	400	*	*	*	*	*	otu2 denovo44 SHOG_9	*
N	*	400	*	*	*	*	*	otu3 denovo45 SHOG_5	*
N	*	400	*	*	*	*	*	otu4 denovo47 SHOG_6	*
D	0	2	*	*	*	*	99.3	denovo0 SHNZ_90	*
D	3	2	*	*	*	*	100.0	denovo3 SHNO_86	*
D	4	2	*	*	*	*	100.0	denovo4 SHNO_85	*
