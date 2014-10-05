# uclust --input /var/folders/x7/d7tds9qn0_74nljddc_n03q40000gn/T/UclustExactMatchFilterUgp_B6.fasta --id 0.97 --tmpdir /var/folders/x7/d7tds9qn0_74nljddc_n03q40000gn/T --w 8 --stepwords 8 --usersort --maxaccepts 1 --stable_sort --maxrejects 8 --uc blank_clustered/test_seqs.blanks_clusters.uc
# version=1.2.22
# Tab-separated fields:
# 1=Type, 2=ClusterNr, 3=SeqLength or ClusterSize, 4=PctId, 5=Strand, 6=QueryStart, 7=SeedStart, 8=Alignment, 9=QueryLabel, 10=TargetLabel
# Record types (field 1): L=LibSeed, S=NewSeed, H=Hit, R=Reject, D=LibCluster, C=NewCluster, N=NoHit
# For C and D types, PctId is average id with seed.
# QueryStart and SeedStart are zero-based relative to start of sequence.
# If minus strand, SeedStart is relative to reverse-complemented seed.
S	0	400	*	*	*	*	*	QiimeExactMatch.SHNZ_505	*
S	1	400	*	*	*	*	*	QiimeExactMatch.SHNO_726	*
S	2	400	*	*	*	*	*	QiimeExactMatch.SHNO_0011	*
H	1	400	99.5	+	0	0	400M	QiimeExactMatch.SHNO_0010	QiimeExactMatch.SHNO_726
H	1	400	99.8	+	0	0	400M	QiimeExactMatch.SHNO_587	QiimeExactMatch.SHNO_726
S	3	400	*	*	*	*	*	QiimeExactMatch.SHNZ_0011	*
H	0	400	99.5	+	0	0	400M	QiimeExactMatch.SHNZ_0010	QiimeExactMatch.SHNZ_505
H	1	400	99.0	+	0	0	2I398M2D	QiimeExactMatch.SHNO_467	QiimeExactMatch.SHNO_726
H	0	400	99.8	+	0	0	400M	QiimeExactMatch.SHNO_401	QiimeExactMatch.SHNZ_505
S	4	400	*	*	*	*	*	QiimeExactMatch.SHNO_0001	*
C	0	3	99.6	*	*	*	*	QiimeExactMatch.SHNZ_505	*
C	1	4	99.4	*	*	*	*	QiimeExactMatch.SHNO_726	*
C	2	1	*	*	*	*	*	QiimeExactMatch.SHNO_0011	*
C	3	1	*	*	*	*	*	QiimeExactMatch.SHNZ_0011	*
C	4	1	*	*	*	*	*	QiimeExactMatch.SHNO_0001	*
