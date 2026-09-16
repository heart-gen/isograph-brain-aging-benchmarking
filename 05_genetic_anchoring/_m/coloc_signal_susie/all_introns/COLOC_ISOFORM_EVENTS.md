# Resolved isoform events — signal-level colocalization (all-introns arm)

Every tissue where a `coloc.susie` nomination's sQTL reaches PP4 >= 0.8, the intron that call used, and the IsoGraph switch it lands on in the same GTEx tissue. Junction matching, switch evidence and the BrainSeq check are the same code as the CLPP layer (`coloc/coloc_isoform_events_combined.parquet`), which this table sits beside and does not replace.

Three differences from the CLPP layer, each limiting what a row can say:

- **No direction.** The risk-allele-signed effect is resolved only for the CLPP layer (`coloc_direction`); `risk_qtl_effect` is empty here rather than borrowed.
- **`estimator = abf` rows name GTEx's representative intron**, the only one the fallback scored. They resolve an event but never tested the gene's other introns.
- **One intron per tissue:** the intron behind that tissue's call. Other introns of the gene that also colocalize are in `signal_pairs.parquet`, not here.

- events (nomination x tissue): **164** over **41** gene x trait nominations
- signal-level (`susie`): **30**; abf fallback: **134**
- junction mapped to a GENCODE v47 transcript: **156/164**
- junction on an IsoGraph switch-pair isoform in the same tissue: **61** events, **21** gene x trait
- also a switch-pair isoform in the matching BrainSeq region: **9**
- cross-tissue exceptions (switch pair taken from another region; flagged, never counted as tissue-concordant, scored apart in the long-read arm): **2** events (UNC13A)

| gene | trait | tissue | junction | PP4 | estimator | prior | transcripts | switch pair | BrainSeq | structural event |
|---|---|---|---|---|---|---|---|---|---|---|
| NDUFS3 | AD | cerebellar_hemisphere | chr11:47565986-47571471(+) | 0.864 | abf | primary_prior | ENST00000533507.5 | no | — | — |
| SIRPA | AD | cerebellar_hemisphere | chr20:1915455-1921395(+) | 0.965 | abf | intermediate | ENST00000356025.7,ENST00000358771.5,ENST00000400068.7,ENST00000622179.4 | yes | — | no annotated structural change |
| SIRPA | AD | cerebellum | chr20:1915455-1921395(+) | 0.965 | abf | intermediate | ENST00000356025.7,ENST00000358771.5,ENST00000400068.7,ENST00000622179.4 | no | — | — |
| SIRPA | AD | anterior_cingulate_cortex_ba24 | chr20:1915455-1921395(+) | 0.965 | abf | intermediate | ENST00000356025.7,ENST00000358771.5,ENST00000400068.7,ENST00000622179.4 | yes | — | no annotated structural change |
| SIRPA | AD | putamen_basal_ganglia | chr20:1915455-1921395(+) | 0.965 | abf | intermediate | ENST00000356025.7,ENST00000358771.5,ENST00000400068.7,ENST00000622179.4 | no | — | — |
| SIRPA | AD | spinal_cord_cervical_c_1 | chr20:1915455-1921395(+) | 0.965 | abf | intermediate | ENST00000356025.7,ENST00000358771.5,ENST00000400068.7,ENST00000622179.4 | no | — | — |
| SIRPA | AD | frontal_cortex_ba9 | chr20:1915455-1921395(+) | 0.964 | abf | intermediate | ENST00000356025.7,ENST00000358771.5,ENST00000400068.7,ENST00000622179.4 | yes | — | no annotated structural change |
| SIRPA | AD | caudate_basal_ganglia | chr20:1915455-1921395(+) | 0.963 | abf | intermediate | ENST00000356025.7,ENST00000358771.5,ENST00000400068.7,ENST00000622179.4 | no | — | — |
| SIRPA | AD | hypothalamus | chr20:1924877-1927875(+) | 0.961 | abf | intermediate | ENST00000356025.7,ENST00000358771.5,ENST00000400068.7,ENST00000622179.4 | no | — | — |
| SIRPA | AD | amygdala | chr20:1915455-1921395(+) | 0.960 | abf | intermediate | ENST00000356025.7,ENST00000358771.5,ENST00000400068.7,ENST00000622179.4 | yes | — | no annotated structural change |
| SIRPA | AD | cortex | chr20:1915455-1921395(+) | 0.957 | abf | intermediate | ENST00000356025.7,ENST00000358771.5,ENST00000400068.7,ENST00000622179.4 | yes | — | no annotated structural change |
| SIRPA | AD | substantia_nigra | chr20:1915455-1921395(+) | 0.955 | abf | intermediate | ENST00000356025.7,ENST00000358771.5,ENST00000400068.7,ENST00000622179.4 | yes | — | no annotated structural change |
| SIRPA | AD | hippocampus | chr20:1915455-1921395(+) | 0.950 | abf | intermediate | ENST00000356025.7,ENST00000358771.5,ENST00000400068.7,ENST00000622179.4 | yes | — | no annotated structural change |
| SIRPA | AD | nucleus_accumbens_basal_ganglia | chr20:1915455-1921395(+) | 0.936 | abf | intermediate | ENST00000356025.7,ENST00000358771.5,ENST00000400068.7,ENST00000622179.4 | yes | — | no annotated structural change |
| TPCN1 | AD | cerebellar_hemisphere | chr12:113280195-113283693(+) | 0.861 | susie | primary_prior | — | — | — | — |
| TPCN1 | AD | cerebellum | chr12:113280195-113284581(+) | 0.835 | susie | primary_prior | ENST00000335509.11,ENST00000392569.8,ENST00000428632.7,ENST00000541517.5,ENST00000546781.5,ENST00000550785.5,ENST00000551127.5,ENST00000552077.5 | yes | — | no annotated structural change |
| ZNF232 | AD | spinal_cord_cervical_c_1 | chr17:5111219-5111800(-) | 0.899 | abf | intermediate | ENST00000570486.5,ENST00000571076.1,ENST00000573015.6,ENST00000575538.1,ENST00000696408.1,ENST00000696538.1 | no | — | — |
| GGNBP2 | ALS | cerebellum | chr17:36560871-36567663(+) | 0.912 | susie | intermediate | ENST00000613102.5,ENST00000617860.4,ENST00000618837.4 | no | — | — |
| PRDM2 | ALS | cortex | chr1:13816570-13823159(+) | 0.964 | abf | intermediate | ENST00000235372.11,ENST00000311066.10,ENST00000376048.9,ENST00000503842.5,ENST00000505823.5 | yes | — | no annotated structural change |
| PRDM2 | ALS | substantia_nigra | chr1:13816570-13823159(+) | 0.947 | abf | intermediate | ENST00000235372.11,ENST00000311066.10,ENST00000376048.9,ENST00000503842.5,ENST00000505823.5 | no | — | — |
| PRDM2 | ALS | nucleus_accumbens_basal_ganglia | chr1:13816570-13823159(+) | 0.945 | abf | intermediate | ENST00000235372.11,ENST00000311066.10,ENST00000376048.9,ENST00000503842.5,ENST00000505823.5 | yes | — | no annotated structural change |
| PRDM2 | ALS | hippocampus | chr1:13816570-13823159(+) | 0.901 | abf | intermediate | ENST00000235372.11,ENST00000311066.10,ENST00000376048.9,ENST00000503842.5,ENST00000505823.5 | no | — | — |
| PRDM2 | ALS | putamen_basal_ganglia | chr1:13816570-13823159(+) | 0.857 | abf | primary_prior | ENST00000235372.11,ENST00000311066.10,ENST00000376048.9,ENST00000503842.5,ENST00000505823.5 | no | — | — |
| PRDM2 | ALS | anterior_cingulate_cortex_ba24 | chr1:13816570-13823159(+) | 0.826 | abf | primary_prior | ENST00000235372.11,ENST00000311066.10,ENST00000376048.9,ENST00000503842.5,ENST00000505823.5 | yes | — | no annotated structural change |
| PTPRN | ALS | nucleus_accumbens_basal_ganglia | chr2:219295141-219295897(-) | 0.918 | abf | intermediate | — | — | — | — |
| PTPRN | ALS | cerebellum | chr2:219295141-219296226(-) | 0.890 | susie | intermediate | ENST00000295718.7,ENST00000409251.7,ENST00000423636.6,ENST00000443981.5 | yes | — | no annotated structural change |
| PTPRN | ALS | cerebellar_hemisphere | chr2:219295141-219296226(-) | 0.890 | susie | intermediate | ENST00000295718.7,ENST00000409251.7,ENST00000423636.6,ENST00000443981.5 | no | — | — |
| PTPRN | ALS | cortex | chr2:219295141-219296226(-) | 0.889 | susie | intermediate | ENST00000295718.7,ENST00000409251.7,ENST00000423636.6,ENST00000443981.5 | no | — | — |
| PTPRN | ALS | frontal_cortex_ba9 | chr2:219295141-219296226(-) | 0.832 | abf | primary_prior | ENST00000295718.7,ENST00000409251.7,ENST00000423636.6,ENST00000443981.5 | no | — | — |
| TXNDC15 | ALS | caudate_basal_ganglia | chr5:134874530-134893492(+) | 0.969 | abf | intermediate | ENST00000511070.5 | no | — | — |
| UNC13A | ALS | cerebellum | chr19:17630750-17632782(-) | 0.959 | susie | intermediate | ENST00000519716.7,ENST00000550896.1,ENST00000551649.5,ENST00000552293.5 | cross-tissue exception (amygdala;anterior_cingulate_cortex_ba24;cortex;frontal_cortex_ba9;hippocampus;nucleus_accumbens_basal_ganglia;substantia_nigra) | — | — |
| UNC13A | ALS | cerebellar_hemisphere | chr19:17630750-17632782(-) | 0.936 | susie | intermediate | ENST00000519716.7,ENST00000550896.1,ENST00000551649.5,ENST00000552293.5 | cross-tissue exception (amygdala;anterior_cingulate_cortex_ba24;cortex;frontal_cortex_ba9;hippocampus;nucleus_accumbens_basal_ganglia;substantia_nigra) | — | — |
| WIPI2 | ALS | nucleus_accumbens_basal_ganglia | chr7:5229705-5230835(+) | 0.893 | abf | intermediate | ENST00000382384.6,ENST00000404704.7,ENST00000484262.1 | yes | — | no annotated structural change |
| CTSB | PD | amygdala | chr8:11844241-11845082(-) | 0.868 | susie | primary_prior | ENST00000530640.7,ENST00000531089.6 | no | — | — |
| DDRGK1 | PD | anterior_cingulate_cortex_ba24 | chr20:3195353-3200001(-) | 0.884 | abf | primary_prior | ENST00000354488.8,ENST00000380201.2 | yes | — | no annotated structural change |
| NCOR1 | PD | cerebellar_hemisphere | chr17:16171995-16186554(-) | 0.895 | susie | intermediate | ENST00000268712.8,ENST00000395851.5,ENST00000411510.5,ENST00000430577.2,ENST00000436068.2,ENST00000436828.5,ENST00000582357.5,ENST00000585296.1,ENST00000704744.1,ENST00000704745.1 | yes | — | no annotated structural change |
| NCOR1 | PD | nucleus_accumbens_basal_ganglia | chr17:16171995-16186554(-) | 0.805 | abf | primary_prior | ENST00000268712.8,ENST00000395851.5,ENST00000411510.5,ENST00000430577.2,ENST00000436068.2,ENST00000436828.5,ENST00000582357.5,ENST00000585296.1,ENST00000704744.1,ENST00000704745.1 | yes | — | no annotated structural change |
| SH3GL2 | PD | substantia_nigra | chr9:17579287-17747066(+) | 0.834 | abf | primary_prior | ENST00000380607.5 | yes | — | no annotated structural change |
| SH3GL2 | PD | cortex | chr9:17579287-17747066(+) | 0.834 | abf | primary_prior | ENST00000380607.5 | yes | — | no annotated structural change |
| SH3GL2 | PD | frontal_cortex_ba9 | chr9:17579287-17747066(+) | 0.834 | abf | primary_prior | ENST00000380607.5 | yes | — | no annotated structural change |
| SH3GL2 | PD | anterior_cingulate_cortex_ba24 | chr9:17579287-17747066(+) | 0.834 | abf | primary_prior | ENST00000380607.5 | yes | — | no annotated structural change |
| SH3GL2 | PD | caudate_basal_ganglia | chr9:17579287-17747066(+) | 0.834 | abf | primary_prior | ENST00000380607.5 | no | — | — |
| SH3GL2 | PD | nucleus_accumbens_basal_ganglia | chr9:17579287-17747066(+) | 0.834 | abf | primary_prior | ENST00000380607.5 | yes | — | no annotated structural change |
| SH3GL2 | PD | hypothalamus | chr9:17579287-17747066(+) | 0.834 | abf | primary_prior | ENST00000380607.5 | no | — | — |
| SH3GL2 | PD | spinal_cord_cervical_c_1 | chr9:17579287-17747066(+) | 0.834 | abf | primary_prior | ENST00000380607.5 | no | — | — |
| SH3GL2 | PD | amygdala | chr9:17579287-17747066(+) | 0.834 | abf | primary_prior | ENST00000380607.5 | yes | — | no annotated structural change |
| TTC19 | PD | putamen_basal_ganglia | chr17:16000245-16001915(+) | 0.898 | susie | intermediate | ENST00000261647.10,ENST00000470399.1,ENST00000497842.6 | no | — | — |
| TTC19 | PD | caudate_basal_ganglia | chr17:16006568-16025017(+) | 0.805 | susie | primary_prior | ENST00000261647.10,ENST00000475723.5 | no | — | — |
| ACTR1B | SCZ | caudate_basal_ganglia | chr2:97658316-97658427(-) | 0.997 | abf | robust | ENST00000289228.7,ENST00000451664.1 | no | yes | — |
| ACTR1B | SCZ | caudate_basal_ganglia | chr2:97658316-97658427(-) | 0.997 | abf | robust | ENST00000289228.7,ENST00000451664.1 | no | yes | — |
| ACTR1B | SCZ | hypothalamus | chr2:97658316-97658427(-) | 0.996 | abf | robust | ENST00000289228.7,ENST00000451664.1 | no | — | — |
| ACTR1B | SCZ | hypothalamus | chr2:97658316-97658427(-) | 0.996 | abf | robust | ENST00000289228.7,ENST00000451664.1 | no | — | — |
| ACTR1B | SCZ | frontal_cortex_ba9 | chr2:97658316-97658427(-) | 0.996 | abf | robust | ENST00000289228.7,ENST00000451664.1 | yes | — | no annotated structural change |
| ACTR1B | SCZ | frontal_cortex_ba9 | chr2:97658316-97658427(-) | 0.996 | abf | robust | ENST00000289228.7,ENST00000451664.1 | yes | — | no annotated structural change |
| ACTR1B | SCZ | nucleus_accumbens_basal_ganglia | chr2:97658316-97658427(-) | 0.995 | abf | robust | ENST00000289228.7,ENST00000451664.1 | yes | — | no annotated structural change |
| ACTR1B | SCZ | nucleus_accumbens_basal_ganglia | chr2:97658316-97658427(-) | 0.995 | abf | robust | ENST00000289228.7,ENST00000451664.1 | yes | — | no annotated structural change |
| ACTR1B | SCZ | cortex | chr2:97658316-97658427(-) | 0.995 | abf | robust | ENST00000289228.7,ENST00000451664.1 | yes | — | no annotated structural change |
| ACTR1B | SCZ | cortex | chr2:97658316-97658427(-) | 0.995 | abf | robust | ENST00000289228.7,ENST00000451664.1 | yes | — | no annotated structural change |
| ACTR1B | SCZ | cerebellar_hemisphere | chr2:97658316-97658427(-) | 0.995 | abf | robust | ENST00000289228.7,ENST00000451664.1 | no | — | — |
| ACTR1B | SCZ | cerebellar_hemisphere | chr2:97658316-97658427(-) | 0.995 | abf | robust | ENST00000289228.7,ENST00000451664.1 | no | — | — |
| ACTR1B | SCZ | cerebellum | chr2:97658316-97658427(-) | 0.995 | abf | robust | ENST00000289228.7,ENST00000451664.1 | no | — | — |
| ACTR1B | SCZ | cerebellum | chr2:97658316-97658427(-) | 0.995 | abf | robust | ENST00000289228.7,ENST00000451664.1 | no | — | — |
| ACTR1B | SCZ | spinal_cord_cervical_c_1 | chr2:97658316-97658427(-) | 0.989 | abf | robust | ENST00000289228.7,ENST00000451664.1 | no | — | — |
| ACTR1B | SCZ | spinal_cord_cervical_c_1 | chr2:97658316-97658427(-) | 0.989 | abf | robust | ENST00000289228.7,ENST00000451664.1 | no | — | — |
| ACTR1B | SCZ | hippocampus | chr2:97658316-97658427(-) | 0.978 | abf | robust | ENST00000289228.7,ENST00000451664.1 | yes | — | no annotated structural change |
| ACTR1B | SCZ | hippocampus | chr2:97658316-97658427(-) | 0.978 | abf | robust | ENST00000289228.7,ENST00000451664.1 | yes | — | no annotated structural change |
| ACTR1B | SCZ | anterior_cingulate_cortex_ba24 | chr2:97658316-97658427(-) | 0.848 | abf | primary_prior | ENST00000289228.7,ENST00000451664.1 | yes | — | no annotated structural change |
| ACTR1B | SCZ | anterior_cingulate_cortex_ba24 | chr2:97658316-97658427(-) | 0.848 | abf | primary_prior | ENST00000289228.7,ENST00000451664.1 | yes | — | no annotated structural change |
| ACTR1B | SCZ | amygdala | chr2:97658316-97658427(-) | 0.811 | abf | primary_prior | ENST00000289228.7,ENST00000451664.1 | yes | — | no annotated structural change |
| ACTR1B | SCZ | amygdala | chr2:97658316-97658427(-) | 0.811 | abf | primary_prior | ENST00000289228.7,ENST00000451664.1 | yes | — | no annotated structural change |
| CDIP1 | SCZ | caudate_basal_ganglia | chr16:4514664-4538325(-) | 0.951 | susie | intermediate | ENST00000563332.6,ENST00000588381.1 | no | yes | — |
| CDIP1 | SCZ | hypothalamus | chr16:4514664-4538325(-) | 0.936 | susie | intermediate | ENST00000563332.6,ENST00000588381.1 | no | — | — |
| CDIP1 | SCZ | amygdala | chr16:4514664-4538702(-) | 0.929 | susie | intermediate | ENST00000562334.5,ENST00000562579.5,ENST00000563507.5,ENST00000564828.5,ENST00000567695.6 | yes | — | no annotated structural change |
| CDIP1 | SCZ | substantia_nigra | chr16:4514664-4538702(-) | 0.928 | abf | intermediate | ENST00000562334.5,ENST00000562579.5,ENST00000563507.5,ENST00000564828.5,ENST00000567695.6 | no | — | — |
| CDIP1 | SCZ | cortex | chr16:4514664-4538702(-) | 0.831 | susie | primary_prior | ENST00000562334.5,ENST00000562579.5,ENST00000563507.5,ENST00000564828.5,ENST00000567695.6 | yes | — | no annotated structural change |
| CDIP1 | SCZ | hippocampus | chr16:4514664-4538702(-) | 0.802 | susie | primary_prior | ENST00000562334.5,ENST00000562579.5,ENST00000563507.5,ENST00000564828.5,ENST00000567695.6 | yes | — | no annotated structural change |
| COPA | SCZ | cerebellum | chr1:160332557-160335242(-) | 0.830 | abf | primary_prior | ENST00000647799.1,ENST00000696207.1 | no | — | — |
| DGKZ | SCZ | frontal_cortex_ba9 | chr11:46369550-46371313(+) | 0.953 | susie | intermediate | ENST00000318201.12,ENST00000531879.5 | no | yes | — |
| DGKZ | SCZ | cortex | chr11:46369550-46371313(+) | 0.951 | susie | intermediate | ENST00000318201.12,ENST00000531879.5 | no | yes | — |
| DGKZ | SCZ | caudate_basal_ganglia | chr11:46369550-46371313(+) | 0.935 | susie | intermediate | ENST00000318201.12,ENST00000531879.5 | no | yes | — |
| DGKZ | SCZ | nucleus_accumbens_basal_ganglia | chr11:46369550-46371313(+) | 0.928 | abf | intermediate | ENST00000318201.12,ENST00000531879.5 | no | — | — |
| DGKZ | SCZ | cerebellar_hemisphere | chr11:46369550-46371313(+) | 0.924 | abf | intermediate | ENST00000318201.12,ENST00000531879.5 | no | — | — |
| DGKZ | SCZ | spinal_cord_cervical_c_1 | chr11:46369550-46371313(+) | 0.892 | abf | intermediate | ENST00000318201.12,ENST00000531879.5 | no | — | — |
| DGKZ | SCZ | anterior_cingulate_cortex_ba24 | chr11:46369550-46371313(+) | 0.882 | susie | primary_prior | ENST00000318201.12,ENST00000531879.5 | yes | — | no annotated structural change |
| DGKZ | SCZ | amygdala | chr11:46369550-46371313(+) | 0.878 | abf | primary_prior | ENST00000318201.12,ENST00000531879.5 | yes | — | no annotated structural change |
| DGKZ | SCZ | hippocampus | chr11:46369550-46371313(+) | 0.838 | abf | primary_prior | ENST00000318201.12,ENST00000531879.5 | yes | yes | no annotated structural change |
| FGFR1 | SCZ | cerebellar_hemisphere | chr8:38400788-38410894(-) | 0.957 | susie | intermediate | — | — | — | — |
| FGFR1 | SCZ | cerebellum | chr8:38400788-38402894(-) | 0.937 | susie | intermediate | — | — | — | — |
| FGFR1 | SCZ | caudate_basal_ganglia | chr8:38429390-38429682(-) | 0.908 | susie | intermediate | ENST00000484370.5 | no | — | — |
| GABBR2 | SCZ | cerebellum | chr9:98299353-98303241(-) | 0.976 | susie | robust | ENST00000259455.4,ENST00000637410.1 | no | — | — |
| GLYCTK | SCZ | anterior_cingulate_cortex_ba24 | chr3:52293192-52293307(+) | 0.812 | abf | primary_prior | ENST00000471180.5,ENST00000473032.5 | no | — | — |
| IRF3 | SCZ | cortex | chr19:49664846-49665672(-) | 0.994 | abf | robust | ENST00000593337.5,ENST00000594387.1,ENST00000601809.5 | no | — | — |
| IRF3 | SCZ | anterior_cingulate_cortex_ba24 | chr19:49664846-49665672(-) | 0.992 | abf | robust | ENST00000593337.5,ENST00000594387.1,ENST00000601809.5 | no | — | — |
| IRF3 | SCZ | spinal_cord_cervical_c_1 | chr19:49664846-49665672(-) | 0.991 | abf | robust | ENST00000593337.5,ENST00000594387.1,ENST00000601809.5 | no | — | — |
| IRF3 | SCZ | putamen_basal_ganglia | chr19:49664846-49665672(-) | 0.991 | abf | robust | ENST00000593337.5,ENST00000594387.1,ENST00000601809.5 | no | — | — |
| IRF3 | SCZ | hypothalamus | chr19:49664846-49665672(-) | 0.990 | abf | robust | ENST00000593337.5,ENST00000594387.1,ENST00000601809.5 | no | — | — |
| IRF3 | SCZ | substantia_nigra | chr19:49664846-49665672(-) | 0.990 | abf | robust | ENST00000593337.5,ENST00000594387.1,ENST00000601809.5 | no | — | — |
| IRF3 | SCZ | nucleus_accumbens_basal_ganglia | chr19:49664846-49665672(-) | 0.990 | abf | robust | ENST00000593337.5,ENST00000594387.1,ENST00000601809.5 | no | — | — |
| IRF3 | SCZ | caudate_basal_ganglia | chr19:49664846-49665672(-) | 0.988 | abf | robust | ENST00000593337.5,ENST00000594387.1,ENST00000601809.5 | no | no | — |
| IRF3 | SCZ | frontal_cortex_ba9 | chr19:49664846-49665672(-) | 0.983 | abf | robust | ENST00000593337.5,ENST00000594387.1,ENST00000601809.5 | no | — | — |
| IRF3 | SCZ | hippocampus | chr19:49664846-49665672(-) | 0.976 | abf | intermediate | ENST00000593337.5,ENST00000594387.1,ENST00000601809.5 | no | — | — |
| IRF3 | SCZ | amygdala | chr19:49664846-49665672(-) | 0.971 | abf | intermediate | ENST00000593337.5,ENST00000594387.1,ENST00000601809.5 | no | — | — |
| KLC1 | SCZ | cortex | chr14:103679518-103684986(+) | 0.890 | susie | intermediate | ENST00000380038.7,ENST00000445352.8,ENST00000553325.5,ENST00000555856.1 | yes | — | no annotated structural change |
| LPCAT4 | SCZ | putamen_basal_ganglia | chr15:34360209-34362196(-) | 0.884 | abf | primary_prior | ENST00000567507.1 | no | — | — |
| MED19 | SCZ | hippocampus | chr11:57704396-57704719(-) | 0.810 | abf | primary_prior | ENST00000431606.5 | yes | — | no annotated structural change |
| MRPS33 | SCZ | hypothalamus | chr7:141006535-141014551(-) | 0.994 | abf | robust | ENST00000484502.1 | no | — | — |
| MRPS33 | SCZ | cerebellum | chr7:141006535-141014551(-) | 0.859 | abf | primary_prior | ENST00000484502.1 | no | — | — |
| MRPS33 | SCZ | nucleus_accumbens_basal_ganglia | chr7:141006535-141014551(-) | 0.848 | abf | primary_prior | ENST00000484502.1 | no | — | — |
| NDUFAF7 | SCZ | putamen_basal_ganglia | chr2:37248426-37253086(+) | 0.909 | susie | intermediate | ENST00000441905.1 | no | — | — |
| NDUFAF7 | SCZ | amygdala | chr2:37248426-37253086(+) | 0.859 | abf | primary_prior | ENST00000441905.1 | yes | — | no annotated structural change |
| NEK4 | SCZ | frontal_cortex_ba9 | chr3:52711869-52737586(-) | 0.820 | abf | primary_prior | ENST00000233027.10,ENST00000535191.5 | yes | — | no annotated structural change |
| NUP50 | SCZ | caudate_basal_ganglia | chr22:45183520-45184453(+) | 0.960 | abf | intermediate | ENST00000347635.9,ENST00000396096.6,ENST00000407019.6,ENST00000493456.5 | no | — | — |
| NUP50 | SCZ | nucleus_accumbens_basal_ganglia | chr22:45183520-45184453(+) | 0.956 | abf | intermediate | ENST00000347635.9,ENST00000396096.6,ENST00000407019.6,ENST00000493456.5 | yes | — | no annotated structural change |
| NUP50 | SCZ | cerebellum | chr22:45183520-45184453(+) | 0.953 | abf | intermediate | ENST00000347635.9,ENST00000396096.6,ENST00000407019.6,ENST00000493456.5 | no | — | — |
| NUP50 | SCZ | cerebellar_hemisphere | chr22:45183520-45184453(+) | 0.945 | abf | intermediate | ENST00000347635.9,ENST00000396096.6,ENST00000407019.6,ENST00000493456.5 | no | — | — |
| NUP50 | SCZ | cortex | chr22:45164296-45168168(+) | 0.935 | abf | intermediate | ENST00000347635.9,ENST00000407019.6,ENST00000417702.5,ENST00000430547.5,ENST00000434760.5,ENST00000484186.5,ENST00000497960.5 | no | — | — |
| NUP50 | SCZ | frontal_cortex_ba9 | chr22:45164296-45168168(+) | 0.927 | abf | intermediate | ENST00000347635.9,ENST00000407019.6,ENST00000417702.5,ENST00000430547.5,ENST00000434760.5,ENST00000484186.5,ENST00000497960.5 | yes | — | no annotated structural change |
| NUP50 | SCZ | hypothalamus | chr22:45164296-45168168(+) | 0.918 | abf | intermediate | ENST00000347635.9,ENST00000407019.6,ENST00000417702.5,ENST00000430547.5,ENST00000434760.5,ENST00000484186.5,ENST00000497960.5 | no | — | — |
| NUP50 | SCZ | anterior_cingulate_cortex_ba24 | chr22:45164296-45168168(+) | 0.918 | abf | intermediate | ENST00000347635.9,ENST00000407019.6,ENST00000417702.5,ENST00000430547.5,ENST00000434760.5,ENST00000484186.5,ENST00000497960.5 | yes | — | no annotated structural change |
| NUP50 | SCZ | hippocampus | chr22:45164296-45168168(+) | 0.910 | abf | intermediate | ENST00000347635.9,ENST00000407019.6,ENST00000417702.5,ENST00000430547.5,ENST00000434760.5,ENST00000484186.5,ENST00000497960.5 | no | — | — |
| NUP50 | SCZ | amygdala | chr22:45183520-45184453(+) | 0.835 | abf | primary_prior | ENST00000347635.9,ENST00000396096.6,ENST00000407019.6,ENST00000493456.5 | yes | — | no annotated structural change |
| NUP50 | SCZ | putamen_basal_ganglia | chr22:45164296-45168168(+) | 0.805 | abf | primary_prior | ENST00000347635.9,ENST00000407019.6,ENST00000417702.5,ENST00000430547.5,ENST00000434760.5,ENST00000484186.5,ENST00000497960.5 | no | — | — |
| POLG | SCZ | cerebellar_hemisphere | chr15:89321012-89321736(-) | 0.982 | abf | robust | — | — | — | — |
| POLG | SCZ | cerebellum | chr15:89321012-89321736(-) | 0.972 | abf | intermediate | — | — | — | — |
| PPIL2 | SCZ | amygdala | chr22:21695496-21696747(+) | 0.943 | abf | intermediate | ENST00000335025.12,ENST00000498109.2,ENST00000679477.1,ENST00000679479.1,ENST00000679564.1,ENST00000679580.1,ENST00000679692.1,ENST00000679827.1,ENST00000680109.1,ENST00000680183.1,ENST00000680434.1,ENST00000680860.1,ENST00000681137.1,ENST00000681286.1,ENST00000681338.1,ENST00000681791.1 | no | — | — |
| PPIL2 | SCZ | cerebellum | chr22:21695496-21696747(+) | 0.941 | abf | intermediate | ENST00000335025.12,ENST00000498109.2,ENST00000679477.1,ENST00000679479.1,ENST00000679564.1,ENST00000679580.1,ENST00000679692.1,ENST00000679827.1,ENST00000680109.1,ENST00000680183.1,ENST00000680434.1,ENST00000680860.1,ENST00000681137.1,ENST00000681286.1,ENST00000681338.1,ENST00000681791.1 | yes | — | no annotated structural change |
| PPIL2 | SCZ | frontal_cortex_ba9 | chr22:21695496-21696747(+) | 0.939 | abf | intermediate | ENST00000335025.12,ENST00000498109.2,ENST00000679477.1,ENST00000679479.1,ENST00000679564.1,ENST00000679580.1,ENST00000679692.1,ENST00000679827.1,ENST00000680109.1,ENST00000680183.1,ENST00000680434.1,ENST00000680860.1,ENST00000681137.1,ENST00000681286.1,ENST00000681338.1,ENST00000681791.1 | no | — | — |
| PPIL2 | SCZ | cerebellar_hemisphere | chr22:21695496-21696747(+) | 0.938 | abf | intermediate | ENST00000335025.12,ENST00000498109.2,ENST00000679477.1,ENST00000679479.1,ENST00000679564.1,ENST00000679580.1,ENST00000679692.1,ENST00000679827.1,ENST00000680109.1,ENST00000680183.1,ENST00000680434.1,ENST00000680860.1,ENST00000681137.1,ENST00000681286.1,ENST00000681338.1,ENST00000681791.1 | no | — | — |
| PPIL2 | SCZ | substantia_nigra | chr22:21695496-21696747(+) | 0.936 | abf | intermediate | ENST00000335025.12,ENST00000498109.2,ENST00000679477.1,ENST00000679479.1,ENST00000679564.1,ENST00000679580.1,ENST00000679692.1,ENST00000679827.1,ENST00000680109.1,ENST00000680183.1,ENST00000680434.1,ENST00000680860.1,ENST00000681137.1,ENST00000681286.1,ENST00000681338.1,ENST00000681791.1 | no | — | — |
| PPIL2 | SCZ | hypothalamus | chr22:21695496-21696747(+) | 0.935 | abf | intermediate | ENST00000335025.12,ENST00000498109.2,ENST00000679477.1,ENST00000679479.1,ENST00000679564.1,ENST00000679580.1,ENST00000679692.1,ENST00000679827.1,ENST00000680109.1,ENST00000680183.1,ENST00000680434.1,ENST00000680860.1,ENST00000681137.1,ENST00000681286.1,ENST00000681338.1,ENST00000681791.1 | no | — | — |
| PPIL2 | SCZ | caudate_basal_ganglia | chr22:21695496-21696747(+) | 0.933 | abf | intermediate | ENST00000335025.12,ENST00000498109.2,ENST00000679477.1,ENST00000679479.1,ENST00000679564.1,ENST00000679580.1,ENST00000679692.1,ENST00000679827.1,ENST00000680109.1,ENST00000680183.1,ENST00000680434.1,ENST00000680860.1,ENST00000681137.1,ENST00000681286.1,ENST00000681338.1,ENST00000681791.1 | no | — | — |
| PPIL2 | SCZ | putamen_basal_ganglia | chr22:21695496-21696747(+) | 0.855 | abf | primary_prior | ENST00000335025.12,ENST00000498109.2,ENST00000679477.1,ENST00000679479.1,ENST00000679564.1,ENST00000679580.1,ENST00000679692.1,ENST00000679827.1,ENST00000680109.1,ENST00000680183.1,ENST00000680434.1,ENST00000680860.1,ENST00000681137.1,ENST00000681286.1,ENST00000681338.1,ENST00000681791.1 | yes | — | no annotated structural change |
| PSMD6 | SCZ | hippocampus | chr3:64010955-64018599(-) | 0.976 | abf | robust | ENST00000480205.5 | no | no | — |
| RAI1 | SCZ | cerebellum | chr17:17793288-17803756(+) | 0.971 | abf | intermediate | — | — | — | — |
| RAI1 | SCZ | cerebellar_hemisphere | chr17:17798513-17803756(+) | 0.945 | susie | intermediate | ENST00000353383.6,ENST00000583166.1 | no | — | — |
| RASA1 | SCZ | cerebellum | chr5:87372195-87374840(+) | 0.946 | abf | intermediate | — | — | — | — |
| RERE | SCZ | nucleus_accumbens_basal_ganglia | chr1:8438300-8465925(-) | 0.856 | abf | primary_prior | ENST00000659924.1 | yes | — | no annotated structural change |
| RERE | SCZ | cerebellar_hemisphere | chr1:8438300-8465925(-) | 0.818 | abf | primary_prior | ENST00000659924.1 | no | — | — |
| SNAP91 | SCZ | cerebellum | chr6:83607808-83610650(-) | 0.974 | susie | intermediate | ENST00000195649.10,ENST00000369694.7,ENST00000439399.6,ENST00000518312.5,ENST00000520213.5,ENST00000520302.5,ENST00000521485.5,ENST00000521616.5,ENST00000521743.5,ENST00000521931.5 | yes | — | no annotated structural change |
| SNAP91 | SCZ | nucleus_accumbens_basal_ganglia | chr6:83607808-83610650(-) | 0.964 | abf | intermediate | ENST00000195649.10,ENST00000369694.7,ENST00000439399.6,ENST00000518312.5,ENST00000520213.5,ENST00000520302.5,ENST00000521485.5,ENST00000521616.5,ENST00000521743.5,ENST00000521931.5 | no | — | — |
| SNAP91 | SCZ | spinal_cord_cervical_c_1 | chr6:83607808-83610650(-) | 0.962 | abf | intermediate | ENST00000195649.10,ENST00000369694.7,ENST00000439399.6,ENST00000518312.5,ENST00000520213.5,ENST00000520302.5,ENST00000521485.5,ENST00000521616.5,ENST00000521743.5,ENST00000521931.5 | no | — | — |
| SNAP91 | SCZ | hypothalamus | chr6:83607808-83610650(-) | 0.961 | abf | intermediate | ENST00000195649.10,ENST00000369694.7,ENST00000439399.6,ENST00000518312.5,ENST00000520213.5,ENST00000520302.5,ENST00000521485.5,ENST00000521616.5,ENST00000521743.5,ENST00000521931.5 | yes | — | no annotated structural change |
| SNAP91 | SCZ | cortex | chr6:83607808-83610650(-) | 0.958 | abf | intermediate | ENST00000195649.10,ENST00000369694.7,ENST00000439399.6,ENST00000518312.5,ENST00000520213.5,ENST00000520302.5,ENST00000521485.5,ENST00000521616.5,ENST00000521743.5,ENST00000521931.5 | no | — | — |
| SNAP91 | SCZ | cerebellar_hemisphere | chr6:83593017-83593478(-) | 0.955 | susie | intermediate | ENST00000518312.5 | no | — | — |
| SNAP91 | SCZ | caudate_basal_ganglia | chr6:83607808-83610650(-) | 0.954 | abf | intermediate | ENST00000195649.10,ENST00000369694.7,ENST00000439399.6,ENST00000518312.5,ENST00000520213.5,ENST00000520302.5,ENST00000521485.5,ENST00000521616.5,ENST00000521743.5,ENST00000521931.5 | no | yes | — |
| SYT5 | SCZ | frontal_cortex_ba9 | chr19:55179086-55179964(-) | 0.904 | abf | intermediate | ENST00000589172.5 | no | — | — |
| TAOK2 | SCZ | cerebellum | chr16:29986504-29986844(+) | 0.829 | abf | primary_prior | ENST00000543033.5 | no | — | — |
| YPEL1 | SCZ | cerebellar_hemisphere | chr22:21701218-21703370(-) | 0.985 | abf | robust | ENST00000339468.8,ENST00000477675.1,ENST00000672036.2 | no | — | — |
| YPEL1 | SCZ | cerebellum | chr22:21701218-21703370(-) | 0.985 | abf | robust | ENST00000339468.8,ENST00000477675.1,ENST00000672036.2 | no | — | — |
| YPEL1 | SCZ | caudate_basal_ganglia | chr22:21701218-21703370(-) | 0.984 | abf | robust | ENST00000339468.8,ENST00000477675.1,ENST00000672036.2 | no | — | — |
| YPEL1 | SCZ | putamen_basal_ganglia | chr22:21701218-21703370(-) | 0.982 | abf | robust | ENST00000339468.8,ENST00000477675.1,ENST00000672036.2 | no | — | — |
| YPEL1 | SCZ | frontal_cortex_ba9 | chr22:21701218-21703370(-) | 0.982 | abf | robust | ENST00000339468.8,ENST00000477675.1,ENST00000672036.2 | yes | — | no annotated structural change |
| YPEL1 | SCZ | nucleus_accumbens_basal_ganglia | chr22:21701218-21703370(-) | 0.904 | abf | intermediate | ENST00000339468.8,ENST00000477675.1,ENST00000672036.2 | yes | — | no annotated structural change |
| YPEL1 | SCZ | anterior_cingulate_cortex_ba24 | chr22:21701218-21703370(-) | 0.899 | abf | intermediate | ENST00000339468.8,ENST00000477675.1,ENST00000672036.2 | yes | — | no annotated structural change |
| YPEL1 | SCZ | cortex | chr22:21701218-21703370(-) | 0.866 | abf | primary_prior | ENST00000339468.8,ENST00000477675.1,ENST00000672036.2 | no | — | — |
| YWHAB | SCZ | frontal_cortex_ba9 | chr20:44885886-44901531(+) | 0.891 | abf | intermediate | ENST00000353703.9,ENST00000479421.5 | yes | — | no annotated structural change |
| YWHAB | SCZ | nucleus_accumbens_basal_ganglia | chr20:44885886-44901531(+) | 0.881 | abf | primary_prior | ENST00000353703.9,ENST00000479421.5 | yes | — | no annotated structural change |
| YWHAB | SCZ | spinal_cord_cervical_c_1 | chr20:44885886-44901531(+) | 0.861 | abf | primary_prior | ENST00000353703.9,ENST00000479421.5 | no | — | — |
| YWHAB | SCZ | cerebellar_hemisphere | chr20:44885886-44901531(+) | 0.860 | abf | primary_prior | ENST00000353703.9,ENST00000479421.5 | no | — | — |
| YWHAB | SCZ | cortex | chr20:44885886-44901531(+) | 0.856 | abf | primary_prior | ENST00000353703.9,ENST00000479421.5 | yes | — | no annotated structural change |
| YWHAB | SCZ | cerebellum | chr20:44885886-44901531(+) | 0.834 | abf | primary_prior | ENST00000353703.9,ENST00000479421.5 | no | — | — |
| YWHAB | SCZ | caudate_basal_ganglia | chr20:44885886-44901531(+) | 0.832 | abf | primary_prior | ENST00000353703.9,ENST00000479421.5 | no | — | — |
| YWHAB | SCZ | hippocampus | chr20:44885886-44901531(+) | 0.803 | abf | primary_prior | ENST00000353703.9,ENST00000479421.5 | yes | yes | no annotated structural change |
| YWHAB | SCZ | hypothalamus | chr20:44885886-44901531(+) | 0.802 | abf | primary_prior | ENST00000353703.9,ENST00000479421.5 | no | — | — |
