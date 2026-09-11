# Resolved isoform events — signal-level colocalization (all-introns arm)

Every tissue where a `coloc.susie` nomination's sQTL reaches PP4 >= 0.8, the intron that call used, and the IsoGraph switch it lands on in the same GTEx tissue. Junction matching, switch evidence and the BrainSeq check are the same code as the CLPP layer (`coloc/coloc_isoform_events_combined.parquet`), which this table sits beside and does not replace.

Three differences from the CLPP layer, each limiting what a row can say:

- **No direction.** The risk-allele-signed effect is resolved only for the CLPP layer (`coloc_direction`); `risk_qtl_effect` is empty here rather than borrowed.
- **`estimator = abf` rows name GTEx's representative intron**, the only one the fallback scored. They resolve an event but never tested the gene's other introns.
- **One intron per tissue:** the intron behind that tissue's call. Other introns of the gene that also colocalize are in `signal_pairs.parquet`, not here.

- events (nomination x tissue): **161** over **42** gene x trait nominations
- signal-level (`susie`): **67**; abf fallback: **94**
- junction mapped to a GENCODE v47 transcript: **147/161**
- junction on an IsoGraph switch-pair isoform in the same tissue: **64** events, **26** gene x trait
- also a switch-pair isoform in the matching BrainSeq region: **17**
- cross-tissue exceptions (switch pair taken from another region; flagged, never counted as tissue-concordant, scored apart in the long-read arm): **2** events (UNC13A)

| gene | trait | tissue | junction | PP4 | estimator | prior | transcripts | switch pair | BrainSeq | structural event |
|---|---|---|---|---|---|---|---|---|---|---|
| CTSH | AD | hippocampus | chr15:78937423-78939140(-) | 0.989 | abf | robust | ENST00000220166.10,ENST00000525807.6,ENST00000527715.6,ENST00000528741.6,ENST00000529612.6,ENST00000529861.6,ENST00000530010.6,ENST00000530929.5,ENST00000534038.6,ENST00000534268.6,ENST00000615999.5,ENST00000649928.2,ENST00000676510.1,ENST00000676596.1,ENST00000676639.1,ENST00000676671.1,ENST00000676850.1,ENST00000676865.1,ENST00000676880.1,ENST00000677238.1,ENST00000677254.1,ENST00000677316.1,ENST00000677320.1,ENST00000677367.1,ENST00000677448.1,ENST00000677534.1,ENST00000677789.1,ENST00000677810.1,ENST00000677874.1,ENST00000677936.1,ENST00000678031.1,ENST00000678281.1,ENST00000678283.1,ENST00000678415.1,ENST00000678727.1,ENST00000678799.1,ENST00000678841.1,ENST00000678886.1,ENST00000678940.1,ENST00000679017.1,ENST00000679047.1,ENST00000679211.1,ENST00000679334.1 | yes | — | no annotated structural change |
| CTSH | AD | nucleus_accumbens_basal_ganglia | chr15:78939171-78944778(-) | 0.988 | abf | robust | ENST00000529861.6,ENST00000677320.1,ENST00000677921.1,ENST00000678799.1,ENST00000678886.1,ENST00000679211.1 | no | — | — |
| CTSH | AD | substantia_nigra | chr15:78937423-78939140(-) | 0.985 | abf | robust | ENST00000220166.10,ENST00000525807.6,ENST00000527715.6,ENST00000528741.6,ENST00000529612.6,ENST00000529861.6,ENST00000530010.6,ENST00000530929.5,ENST00000534038.6,ENST00000534268.6,ENST00000615999.5,ENST00000649928.2,ENST00000676510.1,ENST00000676596.1,ENST00000676639.1,ENST00000676671.1,ENST00000676850.1,ENST00000676865.1,ENST00000676880.1,ENST00000677238.1,ENST00000677254.1,ENST00000677316.1,ENST00000677320.1,ENST00000677367.1,ENST00000677448.1,ENST00000677534.1,ENST00000677789.1,ENST00000677810.1,ENST00000677874.1,ENST00000677936.1,ENST00000678031.1,ENST00000678281.1,ENST00000678283.1,ENST00000678415.1,ENST00000678727.1,ENST00000678799.1,ENST00000678841.1,ENST00000678886.1,ENST00000678940.1,ENST00000679017.1,ENST00000679047.1,ENST00000679211.1,ENST00000679334.1 | no | — | — |
| CTSH | AD | hypothalamus | chr15:78939171-78944778(-) | 0.985 | abf | robust | ENST00000529861.6,ENST00000677320.1,ENST00000677921.1,ENST00000678799.1,ENST00000678886.1,ENST00000679211.1 | no | — | — |
| CTSH | AD | putamen_basal_ganglia | chr15:78939171-78944778(-) | 0.977 | abf | robust | ENST00000529861.6,ENST00000677320.1,ENST00000677921.1,ENST00000678799.1,ENST00000678886.1,ENST00000679211.1 | no | — | — |
| CTSH | AD | spinal_cord_cervical_c_1 | chr15:78939171-78944778(-) | 0.899 | abf | intermediate | ENST00000529861.6,ENST00000677320.1,ENST00000677921.1,ENST00000678799.1,ENST00000678886.1,ENST00000679211.1 | no | — | — |
| DOC2A | AD | amygdala | chr16:30007090-30007173(-) | 0.881 | abf | primary_prior | ENST00000350119.9,ENST00000561671.5,ENST00000564944.5,ENST00000564979.5,ENST00000566310.5,ENST00000616445.4 | yes | — | no annotated structural change |
| NDUFS3 | AD | cerebellar_hemisphere | chr11:47565986-47571471(+) | 0.864 | abf | primary_prior | ENST00000533507.5 | no | — | — |
| PICALM | AD | cortex | chr11:85974812-85981129(-) | 0.818 | abf | primary_prior | ENST00000356360.9,ENST00000532603.5 | no | — | — |
| SIRPA | AD | cerebellar_hemisphere | chr20:1915455-1921395(+) | 0.965 | abf | intermediate | ENST00000356025.7,ENST00000358771.5,ENST00000400068.7,ENST00000622179.4 | no | — | — |
| SIRPA | AD | cerebellum | chr20:1915455-1921395(+) | 0.965 | abf | intermediate | ENST00000356025.7,ENST00000358771.5,ENST00000400068.7,ENST00000622179.4 | yes | — | no annotated structural change |
| SIRPA | AD | anterior_cingulate_cortex_ba24 | chr20:1915455-1921395(+) | 0.965 | abf | intermediate | ENST00000356025.7,ENST00000358771.5,ENST00000400068.7,ENST00000622179.4 | yes | — | no annotated structural change |
| SIRPA | AD | putamen_basal_ganglia | chr20:1915455-1921395(+) | 0.965 | abf | intermediate | ENST00000356025.7,ENST00000358771.5,ENST00000400068.7,ENST00000622179.4 | no | — | — |
| SIRPA | AD | spinal_cord_cervical_c_1 | chr20:1915455-1921395(+) | 0.965 | abf | intermediate | ENST00000356025.7,ENST00000358771.5,ENST00000400068.7,ENST00000622179.4 | no | — | — |
| SIRPA | AD | frontal_cortex_ba9 | chr20:1915455-1921395(+) | 0.964 | abf | intermediate | ENST00000356025.7,ENST00000358771.5,ENST00000400068.7,ENST00000622179.4 | yes | — | no annotated structural change |
| SIRPA | AD | caudate_basal_ganglia | chr20:1915455-1921395(+) | 0.963 | abf | intermediate | ENST00000356025.7,ENST00000358771.5,ENST00000400068.7,ENST00000622179.4 | no | yes | — |
| SIRPA | AD | hypothalamus | chr20:1924877-1927875(+) | 0.961 | abf | intermediate | ENST00000356025.7,ENST00000358771.5,ENST00000400068.7,ENST00000622179.4 | no | — | — |
| SIRPA | AD | amygdala | chr20:1915455-1921395(+) | 0.960 | abf | intermediate | ENST00000356025.7,ENST00000358771.5,ENST00000400068.7,ENST00000622179.4 | yes | — | no annotated structural change |
| SIRPA | AD | cortex | chr20:1915455-1921395(+) | 0.957 | abf | intermediate | ENST00000356025.7,ENST00000358771.5,ENST00000400068.7,ENST00000622179.4 | yes | — | no annotated structural change |
| SIRPA | AD | substantia_nigra | chr20:1915455-1921395(+) | 0.955 | abf | intermediate | ENST00000356025.7,ENST00000358771.5,ENST00000400068.7,ENST00000622179.4 | no | — | — |
| SIRPA | AD | hippocampus | chr20:1915455-1921395(+) | 0.950 | abf | intermediate | ENST00000356025.7,ENST00000358771.5,ENST00000400068.7,ENST00000622179.4 | yes | yes | no annotated structural change |
| SIRPA | AD | nucleus_accumbens_basal_ganglia | chr20:1915455-1921395(+) | 0.936 | abf | intermediate | ENST00000356025.7,ENST00000358771.5,ENST00000400068.7,ENST00000622179.4 | no | — | — |
| SPAG9 | AD | hypothalamus | chr17:50974538-50974771(-) | 0.839 | abf | primary_prior | — | — | — | — |
| SPAG9 | AD | putamen_basal_ganglia | chr17:50974538-50974771(-) | 0.827 | abf | primary_prior | — | — | — | — |
| SPAG9 | AD | hippocampus | chr17:50974538-50974771(-) | 0.823 | abf | primary_prior | — | — | — | — |
| SPAG9 | AD | anterior_cingulate_cortex_ba24 | chr17:50974538-50974771(-) | 0.823 | abf | primary_prior | — | — | — | — |
| SPAG9 | AD | substantia_nigra | chr17:50974538-50974771(-) | 0.821 | abf | primary_prior | — | — | — | — |
| SPAG9 | AD | amygdala | chr17:50974538-50974771(-) | 0.811 | abf | primary_prior | — | — | — | — |
| SPAG9 | AD | caudate_basal_ganglia | chr17:50974538-50974771(-) | 0.809 | abf | primary_prior | — | — | — | — |
| SPAG9 | AD | cortex | chr17:50970856-50974771(-) | 0.802 | abf | primary_prior | ENST00000262013.12,ENST00000357122.8,ENST00000505279.5,ENST00000506500.1,ENST00000510283.5 | no | yes | — |
| SPI1 | AD | cerebellar_hemisphere | chr11:47383832-47408952(-) | 0.926 | abf | intermediate | ENST00000713543.1 | no | — | — |
| TPCN1 | AD | cerebellar_hemisphere | chr12:113280195-113283693(+) | 0.861 | susie | primary_prior | — | — | — | — |
| TPCN1 | AD | cerebellum | chr12:113280195-113284581(+) | 0.835 | susie | primary_prior | ENST00000335509.11,ENST00000392569.8,ENST00000428632.7,ENST00000541517.5,ENST00000546781.5,ENST00000550785.5,ENST00000551127.5,ENST00000552077.5 | yes | — | no annotated structural change |
| ZNF232 | AD | spinal_cord_cervical_c_1 | chr17:5111219-5111800(-) | 0.899 | abf | intermediate | ENST00000570486.5,ENST00000571076.1,ENST00000573015.6,ENST00000575538.1,ENST00000696408.1,ENST00000696538.1 | no | — | — |
| G2E3 | ALS | nucleus_accumbens_basal_ganglia | chr14:30581116-30586718(+) | 0.839 | abf | primary_prior | ENST00000206595.11,ENST00000547532.5,ENST00000550944.5,ENST00000552488.1,ENST00000553504.5,ENST00000554714.5,ENST00000555429.1,ENST00000648008.1 | yes | — | no annotated structural change |
| GGNBP2 | ALS | cerebellum | chr17:36560871-36567663(+) | 0.911 | susie | intermediate | ENST00000613102.5,ENST00000617860.4,ENST00000618837.4 | yes | — | no annotated structural change |
| PGS1 | ALS | cortex | chr17:78378808-78392476(+) | 0.976 | abf | robust | ENST00000262764.11,ENST00000585521.2,ENST00000589426.5,ENST00000589689.5,ENST00000592043.5 | yes | — | no annotated structural change |
| PGS1 | ALS | hippocampus | chr17:78378808-78392476(+) | 0.974 | abf | intermediate | ENST00000262764.11,ENST00000585521.2,ENST00000589426.5,ENST00000589689.5,ENST00000592043.5 | no | — | — |
| PGS1 | ALS | amygdala | chr17:78378808-78392476(+) | 0.973 | abf | intermediate | ENST00000262764.11,ENST00000585521.2,ENST00000589426.5,ENST00000589689.5,ENST00000592043.5 | no | — | — |
| PGS1 | ALS | cerebellar_hemisphere | chr17:78378808-78392476(+) | 0.972 | abf | intermediate | ENST00000262764.11,ENST00000585521.2,ENST00000589426.5,ENST00000589689.5,ENST00000592043.5 | yes | — | no annotated structural change |
| PGS1 | ALS | substantia_nigra | chr17:78378808-78392476(+) | 0.970 | abf | intermediate | ENST00000262764.11,ENST00000585521.2,ENST00000589426.5,ENST00000589689.5,ENST00000592043.5 | no | — | — |
| PGS1 | ALS | putamen_basal_ganglia | chr17:78378808-78392476(+) | 0.970 | abf | intermediate | ENST00000262764.11,ENST00000585521.2,ENST00000589426.5,ENST00000589689.5,ENST00000592043.5 | no | — | — |
| PGS1 | ALS | cerebellum | chr17:78378808-78392476(+) | 0.969 | abf | intermediate | ENST00000262764.11,ENST00000585521.2,ENST00000589426.5,ENST00000589689.5,ENST00000592043.5 | yes | — | no annotated structural change |
| PGS1 | ALS | spinal_cord_cervical_c_1 | chr17:78378808-78392476(+) | 0.968 | abf | intermediate | ENST00000262764.11,ENST00000585521.2,ENST00000589426.5,ENST00000589689.5,ENST00000592043.5 | no | — | — |
| PGS1 | ALS | caudate_basal_ganglia | chr17:78378808-78392476(+) | 0.967 | abf | intermediate | ENST00000262764.11,ENST00000585521.2,ENST00000589426.5,ENST00000589689.5,ENST00000592043.5 | no | — | — |
| PGS1 | ALS | hypothalamus | chr17:78378808-78392476(+) | 0.954 | abf | intermediate | ENST00000262764.11,ENST00000585521.2,ENST00000589426.5,ENST00000589689.5,ENST00000592043.5 | no | — | — |
| PGS1 | ALS | nucleus_accumbens_basal_ganglia | chr17:78378808-78392476(+) | 0.945 | abf | intermediate | ENST00000262764.11,ENST00000585521.2,ENST00000589426.5,ENST00000589689.5,ENST00000592043.5 | yes | — | no annotated structural change |
| PGS1 | ALS | frontal_cortex_ba9 | chr17:78378808-78392476(+) | 0.942 | abf | intermediate | ENST00000262764.11,ENST00000585521.2,ENST00000589426.5,ENST00000589689.5,ENST00000592043.5 | yes | — | no annotated structural change |
| PGS1 | ALS | anterior_cingulate_cortex_ba24 | chr17:78378808-78392476(+) | 0.906 | abf | intermediate | ENST00000262764.11,ENST00000585521.2,ENST00000589426.5,ENST00000589689.5,ENST00000592043.5 | no | — | — |
| PTPRN | ALS | nucleus_accumbens_basal_ganglia | chr2:219295141-219295897(-) | 0.918 | abf | intermediate | — | — | — | — |
| PTPRN | ALS | cerebellar_hemisphere | chr2:219295966-219296226(-) | 0.876 | susie | primary_prior | — | — | — | — |
| PTPRN | ALS | cerebellum | chr2:219295141-219296226(-) | 0.875 | susie | primary_prior | ENST00000295718.7,ENST00000409251.7,ENST00000423636.6,ENST00000443981.5 | no | — | — |
| PTPRN | ALS | cortex | chr2:219295141-219296226(-) | 0.874 | susie | primary_prior | ENST00000295718.7,ENST00000409251.7,ENST00000423636.6,ENST00000443981.5 | no | yes | — |
| PTPRN | ALS | frontal_cortex_ba9 | chr2:219295141-219296226(-) | 0.832 | abf | primary_prior | ENST00000295718.7,ENST00000409251.7,ENST00000423636.6,ENST00000443981.5 | yes | yes | no annotated structural change |
| SCFD1 | ALS | hypothalamus | chr14:30622399-30630477(+) | 0.856 | abf | primary_prior | ENST00000484733.6,ENST00000544052.6,ENST00000553693.6,ENST00000554437.6,ENST00000676465.1,ENST00000676473.1,ENST00000676520.1,ENST00000676674.1,ENST00000676834.1,ENST00000676914.1,ENST00000676954.1,ENST00000677413.1,ENST00000677456.1,ENST00000677690.1,ENST00000678579.1,ENST00000678637.1,ENST00000678858.1,ENST00000679342.1 | yes | — | no annotated structural change |
| TPP1 | ALS | cerebellum | chr11:6611818-6611958(-) | 0.996 | abf | robust | — | — | — | — |
| TXNDC15 | ALS | caudate_basal_ganglia | chr5:134874530-134893492(+) | 0.969 | abf | intermediate | ENST00000511070.5 | no | — | — |
| UNC13A | ALS | cerebellum | chr19:17630750-17632782(-) | 0.961 | susie | intermediate | ENST00000519716.7,ENST00000550896.1,ENST00000551649.5,ENST00000552293.5 | cross-tissue exception (anterior_cingulate_cortex_ba24;caudate_basal_ganglia;frontal_cortex_ba9;hippocampus;nucleus_accumbens_basal_ganglia) | — | — |
| UNC13A | ALS | cerebellar_hemisphere | chr19:17630750-17632782(-) | 0.939 | susie | intermediate | ENST00000519716.7,ENST00000550896.1,ENST00000551649.5,ENST00000552293.5 | cross-tissue exception (anterior_cingulate_cortex_ba24;caudate_basal_ganglia;frontal_cortex_ba9;hippocampus;nucleus_accumbens_basal_ganglia) | — | — |
| WIPI2 | ALS | nucleus_accumbens_basal_ganglia | chr7:5229705-5230835(+) | 0.893 | abf | intermediate | ENST00000382384.6,ENST00000404704.7,ENST00000484262.1 | no | — | — |
| SNCA | LBD | cortex | chr4:89835692-89836127(-) | 0.974 | susie | intermediate | ENST00000508895.5,ENST00000618500.4 | yes | — | no annotated structural change |
| SNCA | LBD | substantia_nigra | chr4:89835692-89836127(-) | 0.971 | susie | intermediate | ENST00000508895.5,ENST00000618500.4 | no | — | — |
| SNCA | LBD | hypothalamus | chr4:89835692-89836127(-) | 0.971 | susie | intermediate | ENST00000508895.5,ENST00000618500.4 | no | — | — |
| SNCA | LBD | spinal_cord_cervical_c_1 | chr4:89835692-89836127(-) | 0.962 | susie | intermediate | ENST00000508895.5,ENST00000618500.4 | no | — | — |
| SNCA | LBD | frontal_cortex_ba9 | chr4:89835692-89836127(-) | 0.956 | susie | intermediate | ENST00000508895.5,ENST00000618500.4 | yes | — | no annotated structural change |
| SNCA | LBD | nucleus_accumbens_basal_ganglia | chr4:89835692-89836127(-) | 0.956 | susie | intermediate | ENST00000508895.5,ENST00000618500.4 | no | — | — |
| SNCA | LBD | amygdala | chr4:89835692-89836127(-) | 0.955 | susie | intermediate | ENST00000508895.5,ENST00000618500.4 | yes | — | no annotated structural change |
| SNCA | LBD | hippocampus | chr4:89835692-89836127(-) | 0.954 | susie | intermediate | ENST00000508895.5,ENST00000618500.4 | yes | yes | no annotated structural change |
| SNCA | LBD | cerebellar_hemisphere | chr4:89835692-89836127(-) | 0.952 | susie | intermediate | ENST00000508895.5,ENST00000618500.4 | no | — | — |
| SNCA | LBD | anterior_cingulate_cortex_ba24 | chr4:89835692-89836127(-) | 0.952 | susie | intermediate | ENST00000508895.5,ENST00000618500.4 | no | — | — |
| SNCA | LBD | cerebellum | chr4:89835692-89836127(-) | 0.941 | susie | intermediate | ENST00000508895.5,ENST00000618500.4 | yes | — | no annotated structural change |
| AZI2 | PD | cerebellar_hemisphere | chr3:28331954-28332369(-) | 0.803 | abf | primary_prior | ENST00000420543.6,ENST00000457172.5,ENST00000488978.1 | yes | — | no annotated structural change |
| NCOR1 | PD | cerebellar_hemisphere | chr17:16171995-16186554(-) | 0.895 | susie | intermediate | ENST00000268712.8,ENST00000395851.5,ENST00000411510.5,ENST00000430577.2,ENST00000436068.2,ENST00000436828.5,ENST00000582357.5,ENST00000585296.1,ENST00000704744.1,ENST00000704745.1 | no | — | — |
| NCOR1 | PD | nucleus_accumbens_basal_ganglia | chr17:16171995-16186554(-) | 0.805 | abf | primary_prior | ENST00000268712.8,ENST00000395851.5,ENST00000411510.5,ENST00000430577.2,ENST00000436068.2,ENST00000436828.5,ENST00000582357.5,ENST00000585296.1,ENST00000704744.1,ENST00000704745.1 | yes | — | no annotated structural change |
| PITPNM2 | PD | cerebellar_hemisphere | chr12:122996910-122997325(-) | 0.989 | abf | robust | ENST00000280562.9,ENST00000320201.10 | no | — | — |
| SH3GL2 | PD | substantia_nigra | chr9:17579287-17747066(+) | 0.834 | abf | primary_prior | ENST00000380607.5 | no | — | — |
| SH3GL2 | PD | cortex | chr9:17579287-17747066(+) | 0.834 | abf | primary_prior | ENST00000380607.5 | no | yes | — |
| SH3GL2 | PD | frontal_cortex_ba9 | chr9:17579287-17747066(+) | 0.834 | abf | primary_prior | ENST00000380607.5 | yes | yes | no annotated structural change |
| SH3GL2 | PD | anterior_cingulate_cortex_ba24 | chr9:17579287-17747066(+) | 0.834 | abf | primary_prior | ENST00000380607.5 | yes | — | no annotated structural change |
| SH3GL2 | PD | caudate_basal_ganglia | chr9:17579287-17747066(+) | 0.834 | abf | primary_prior | ENST00000380607.5 | yes | — | no annotated structural change |
| SH3GL2 | PD | nucleus_accumbens_basal_ganglia | chr9:17579287-17747066(+) | 0.834 | abf | primary_prior | ENST00000380607.5 | no | — | — |
| SH3GL2 | PD | hypothalamus | chr9:17579287-17747066(+) | 0.834 | abf | primary_prior | ENST00000380607.5 | no | — | — |
| SH3GL2 | PD | spinal_cord_cervical_c_1 | chr9:17579287-17747066(+) | 0.834 | abf | primary_prior | ENST00000380607.5 | no | — | — |
| SH3GL2 | PD | amygdala | chr9:17579287-17747066(+) | 0.834 | abf | primary_prior | ENST00000380607.5 | no | — | — |
| SNCA | PD | cerebellum | chr4:89835692-89836127(-) | 0.974 | susie | intermediate | ENST00000508895.5,ENST00000618500.4 | yes | — | no annotated structural change |
| SNCA | PD | substantia_nigra | chr4:89835692-89836127(-) | 0.970 | susie | intermediate | ENST00000508895.5,ENST00000618500.4 | no | — | — |
| SNCA | PD | frontal_cortex_ba9 | chr4:89835692-89836127(-) | 0.968 | susie | intermediate | ENST00000508895.5,ENST00000618500.4 | yes | — | no annotated structural change |
| SNCA | PD | cortex | chr4:89835692-89836127(-) | 0.959 | susie | intermediate | ENST00000508895.5,ENST00000618500.4 | yes | — | no annotated structural change |
| SNCA | PD | nucleus_accumbens_basal_ganglia | chr4:89835692-89836127(-) | 0.950 | susie | intermediate | ENST00000508895.5,ENST00000618500.4 | no | — | — |
| SNCA | PD | hypothalamus | chr4:89835692-89836127(-) | 0.903 | susie | intermediate | ENST00000508895.5,ENST00000618500.4 | no | — | — |
| SNCA | PD | hippocampus | chr4:89835692-89836127(-) | 0.894 | susie | intermediate | ENST00000508895.5,ENST00000618500.4 | yes | yes | no annotated structural change |
| SNCA | PD | spinal_cord_cervical_c_1 | chr4:89835692-89836127(-) | 0.894 | susie | intermediate | ENST00000508895.5,ENST00000618500.4 | no | — | — |
| TTC19 | PD | putamen_basal_ganglia | chr17:16000245-16001915(+) | 0.898 | susie | intermediate | ENST00000261647.10,ENST00000470399.1,ENST00000497842.6 | no | — | — |
| TTC19 | PD | caudate_basal_ganglia | chr17:16006568-16025017(+) | 0.805 | susie | primary_prior | ENST00000261647.10,ENST00000475723.5 | no | — | — |
| ASB3 | SCZ | cerebellar_hemisphere | chr2:53729570-53750783(-) | 0.859 | abf | primary_prior | ENST00000263634.8,ENST00000394717.3,ENST00000406053.5,ENST00000406625.6,ENST00000406687.5,ENST00000482829.5,ENST00000489508.5 | yes | — | no annotated structural change |
| CDIP1 | SCZ | caudate_basal_ganglia | chr16:4514664-4538325(-) | 0.949 | susie | intermediate | ENST00000563332.6,ENST00000588381.1 | yes | — | no annotated structural change |
| CDIP1 | SCZ | hypothalamus | chr16:4514664-4538325(-) | 0.934 | susie | intermediate | ENST00000563332.6,ENST00000588381.1 | yes | — | no annotated structural change |
| CDIP1 | SCZ | substantia_nigra | chr16:4514664-4538702(-) | 0.928 | abf | intermediate | ENST00000562334.5,ENST00000562579.5,ENST00000563507.5,ENST00000564828.5,ENST00000567695.6 | yes | — | no annotated structural change |
| CDIP1 | SCZ | amygdala | chr16:4514664-4538702(-) | 0.928 | susie | intermediate | ENST00000562334.5,ENST00000562579.5,ENST00000563507.5,ENST00000564828.5,ENST00000567695.6 | yes | — | no annotated structural change |
| CDIP1 | SCZ | cortex | chr16:4514664-4538702(-) | 0.842 | susie | primary_prior | ENST00000562334.5,ENST00000562579.5,ENST00000563507.5,ENST00000564828.5,ENST00000567695.6 | yes | — | no annotated structural change |
| CDIP1 | SCZ | hippocampus | chr16:4514664-4538702(-) | 0.804 | susie | primary_prior | ENST00000562334.5,ENST00000562579.5,ENST00000563507.5,ENST00000564828.5,ENST00000567695.6 | yes | — | no annotated structural change |
| COPA | SCZ | cerebellum | chr1:160332557-160335242(-) | 0.830 | abf | primary_prior | ENST00000647799.1,ENST00000696207.1 | no | — | — |
| DOC2A | SCZ | frontal_cortex_ba9 | chr16:30007090-30007179(-) | 0.970 | susie | intermediate | — | — | — | — |
| DOC2A | SCZ | cerebellar_hemisphere | chr16:30007299-30008996(-) | 0.969 | susie | intermediate | ENST00000350119.9,ENST00000563378.5,ENST00000564944.5,ENST00000564979.5,ENST00000565273.5,ENST00000566310.5,ENST00000567332.6,ENST00000616445.4 | yes | — | no annotated structural change |
| DOC2A | SCZ | cortex | chr16:30007090-30007179(-) | 0.967 | susie | intermediate | — | — | — | — |
| DOC2A | SCZ | hippocampus | chr16:30007090-30007173(-) | 0.964 | abf | intermediate | ENST00000350119.9,ENST00000561671.5,ENST00000564944.5,ENST00000564979.5,ENST00000566310.5,ENST00000616445.4 | yes | yes | no annotated structural change |
| DOC2A | SCZ | anterior_cingulate_cortex_ba24 | chr16:30007090-30007173(-) | 0.959 | abf | intermediate | ENST00000350119.9,ENST00000561671.5,ENST00000564944.5,ENST00000564979.5,ENST00000566310.5,ENST00000616445.4 | no | — | — |
| DOC2A | SCZ | hypothalamus | chr16:30007090-30007173(-) | 0.949 | abf | intermediate | ENST00000350119.9,ENST00000561671.5,ENST00000564944.5,ENST00000564979.5,ENST00000566310.5,ENST00000616445.4 | yes | — | no annotated structural change |
| DOC2A | SCZ | amygdala | chr16:30007090-30007173(-) | 0.947 | abf | intermediate | ENST00000350119.9,ENST00000561671.5,ENST00000564944.5,ENST00000564979.5,ENST00000566310.5,ENST00000616445.4 | yes | — | no annotated structural change |
| DOC2A | SCZ | nucleus_accumbens_basal_ganglia | chr16:30007090-30007173(-) | 0.901 | abf | intermediate | ENST00000350119.9,ENST00000561671.5,ENST00000564944.5,ENST00000564979.5,ENST00000566310.5,ENST00000616445.4 | no | — | — |
| DOC2A | SCZ | caudate_basal_ganglia | chr16:30007090-30007173(-) | 0.871 | abf | primary_prior | ENST00000350119.9,ENST00000561671.5,ENST00000564944.5,ENST00000564979.5,ENST00000566310.5,ENST00000616445.4 | yes | — | no annotated structural change |
| FAM221A | SCZ | anterior_cingulate_cortex_ba24 | chr7:23691596-23698192(+) | 0.817 | susie | primary_prior | ENST00000344962.9,ENST00000409653.5 | yes | — | no annotated structural change |
| FAM221A | SCZ | nucleus_accumbens_basal_ganglia | chr7:23698299-23700789(+) | 0.812 | susie | primary_prior | — | — | — | — |
| FAM221A | SCZ | cortex | chr7:23691596-23700786(+) | 0.803 | susie | primary_prior | ENST00000409192.7,ENST00000409994.3 | yes | — | no annotated structural change |
| GABBR2 | SCZ | cerebellum | chr9:98299353-98303241(-) | 0.976 | susie | robust | ENST00000259455.4,ENST00000637410.1 | no | — | — |
| GPM6A | SCZ | caudate_basal_ganglia | chr4:175701767-176002309(-) | 0.991 | susie | robust | ENST00000506894.5 | no | — | — |
| GPM6A | SCZ | anterior_cingulate_cortex_ba24 | chr4:175701767-176002309(-) | 0.990 | susie | robust | ENST00000506894.5 | yes | — | no annotated structural change |
| GPM6A | SCZ | putamen_basal_ganglia | chr4:175701767-175812191(-) | 0.908 | abf | intermediate | ENST00000280187.11,ENST00000393658.7,ENST00000513365.1 | no | — | — |
| GPM6A | SCZ | amygdala | chr4:175701767-175812191(-) | 0.841 | abf | primary_prior | ENST00000280187.11,ENST00000393658.7,ENST00000513365.1 | yes | — | no annotated structural change |
| HMOX2 | SCZ | frontal_cortex_ba9 | chr16:4483754-4505484(+) | 0.831 | susie | primary_prior | ENST00000458134.7,ENST00000570445.5 | yes | — | no annotated structural change |
| KLC1 | SCZ | cortex | chr14:103679518-103684986(+) | 0.882 | susie | primary_prior | ENST00000380038.7,ENST00000445352.8,ENST00000553325.5,ENST00000555856.1 | no | — | — |
| NEK4 | SCZ | frontal_cortex_ba9 | chr3:52711869-52737586(-) | 0.820 | abf | primary_prior | ENST00000233027.10,ENST00000535191.5 | yes | yes | no annotated structural change |
| NT5C2 | SCZ | cerebellum | chr10:103174982-103181185(-) | 0.958 | susie | intermediate | ENST00000404739.8,ENST00000467380.1,ENST00000674860.1,ENST00000675020.1,ENST00000675326.1,ENST00000675436.1 | no | — | — |
| NT5C2 | SCZ | cerebellar_hemisphere | chr10:103174982-103181185(-) | 0.943 | susie | intermediate | ENST00000404739.8,ENST00000467380.1,ENST00000674860.1,ENST00000675020.1,ENST00000675326.1,ENST00000675436.1 | yes | — | no annotated structural change |
| NT5C2 | SCZ | putamen_basal_ganglia | chr10:103174982-103181185(-) | 0.907 | abf | intermediate | ENST00000404739.8,ENST00000467380.1,ENST00000674860.1,ENST00000675020.1,ENST00000675326.1,ENST00000675436.1 | no | — | — |
| NT5C2 | SCZ | cortex | chr10:103174982-103181185(-) | 0.906 | susie | intermediate | ENST00000404739.8,ENST00000467380.1,ENST00000674860.1,ENST00000675020.1,ENST00000675326.1,ENST00000675436.1 | yes | yes | no annotated structural change |
| NT5C2 | SCZ | caudate_basal_ganglia | chr10:103174982-103181185(-) | 0.901 | susie | intermediate | ENST00000404739.8,ENST00000467380.1,ENST00000674860.1,ENST00000675020.1,ENST00000675326.1,ENST00000675436.1 | no | — | — |
| NT5C2 | SCZ | amygdala | chr10:103174982-103181185(-) | 0.823 | abf | primary_prior | ENST00000404739.8,ENST00000467380.1,ENST00000674860.1,ENST00000675020.1,ENST00000675326.1,ENST00000675436.1 | no | — | — |
| NT5C2 | SCZ | frontal_cortex_ba9 | chr10:103174982-103181185(-) | 0.822 | abf | primary_prior | ENST00000404739.8,ENST00000467380.1,ENST00000674860.1,ENST00000675020.1,ENST00000675326.1,ENST00000675436.1 | yes | yes | no annotated structural change |
| PLCB2 | SCZ | hypothalamus | chr15:40288642-40289272(-) | 0.975 | abf | intermediate | ENST00000558505.5 | yes | — | no annotated structural change |
| PLCB2 | SCZ | nucleus_accumbens_basal_ganglia | chr15:40288642-40289272(-) | 0.826 | susie | primary_prior | ENST00000558505.5 | no | — | — |
| PPIL2 | SCZ | amygdala | chr22:21695496-21696747(+) | 0.943 | abf | intermediate | ENST00000335025.12,ENST00000498109.2,ENST00000679477.1,ENST00000679479.1,ENST00000679564.1,ENST00000679580.1,ENST00000679692.1,ENST00000679827.1,ENST00000680109.1,ENST00000680183.1,ENST00000680434.1,ENST00000680860.1,ENST00000681137.1,ENST00000681286.1,ENST00000681338.1,ENST00000681791.1 | no | — | — |
| PPIL2 | SCZ | cerebellum | chr22:21695496-21696747(+) | 0.941 | abf | intermediate | ENST00000335025.12,ENST00000498109.2,ENST00000679477.1,ENST00000679479.1,ENST00000679564.1,ENST00000679580.1,ENST00000679692.1,ENST00000679827.1,ENST00000680109.1,ENST00000680183.1,ENST00000680434.1,ENST00000680860.1,ENST00000681137.1,ENST00000681286.1,ENST00000681338.1,ENST00000681791.1 | yes | — | no annotated structural change |
| PPIL2 | SCZ | frontal_cortex_ba9 | chr22:21695496-21696747(+) | 0.939 | abf | intermediate | ENST00000335025.12,ENST00000498109.2,ENST00000679477.1,ENST00000679479.1,ENST00000679564.1,ENST00000679580.1,ENST00000679692.1,ENST00000679827.1,ENST00000680109.1,ENST00000680183.1,ENST00000680434.1,ENST00000680860.1,ENST00000681137.1,ENST00000681286.1,ENST00000681338.1,ENST00000681791.1 | yes | — | no annotated structural change |
| PPIL2 | SCZ | cerebellar_hemisphere | chr22:21695496-21696747(+) | 0.938 | abf | intermediate | ENST00000335025.12,ENST00000498109.2,ENST00000679477.1,ENST00000679479.1,ENST00000679564.1,ENST00000679580.1,ENST00000679692.1,ENST00000679827.1,ENST00000680109.1,ENST00000680183.1,ENST00000680434.1,ENST00000680860.1,ENST00000681137.1,ENST00000681286.1,ENST00000681338.1,ENST00000681791.1 | yes | — | no annotated structural change |
| PPIL2 | SCZ | substantia_nigra | chr22:21695496-21696747(+) | 0.936 | abf | intermediate | ENST00000335025.12,ENST00000498109.2,ENST00000679477.1,ENST00000679479.1,ENST00000679564.1,ENST00000679580.1,ENST00000679692.1,ENST00000679827.1,ENST00000680109.1,ENST00000680183.1,ENST00000680434.1,ENST00000680860.1,ENST00000681137.1,ENST00000681286.1,ENST00000681338.1,ENST00000681791.1 | yes | — | no annotated structural change |
| PPIL2 | SCZ | hypothalamus | chr22:21695496-21696747(+) | 0.935 | abf | intermediate | ENST00000335025.12,ENST00000498109.2,ENST00000679477.1,ENST00000679479.1,ENST00000679564.1,ENST00000679580.1,ENST00000679692.1,ENST00000679827.1,ENST00000680109.1,ENST00000680183.1,ENST00000680434.1,ENST00000680860.1,ENST00000681137.1,ENST00000681286.1,ENST00000681338.1,ENST00000681791.1 | no | — | — |
| PPIL2 | SCZ | caudate_basal_ganglia | chr22:21695496-21696747(+) | 0.933 | abf | intermediate | ENST00000335025.12,ENST00000498109.2,ENST00000679477.1,ENST00000679479.1,ENST00000679564.1,ENST00000679580.1,ENST00000679692.1,ENST00000679827.1,ENST00000680109.1,ENST00000680183.1,ENST00000680434.1,ENST00000680860.1,ENST00000681137.1,ENST00000681286.1,ENST00000681338.1,ENST00000681791.1 | no | — | — |
| PPIL2 | SCZ | putamen_basal_ganglia | chr22:21695496-21696747(+) | 0.855 | abf | primary_prior | ENST00000335025.12,ENST00000498109.2,ENST00000679477.1,ENST00000679479.1,ENST00000679564.1,ENST00000679580.1,ENST00000679692.1,ENST00000679827.1,ENST00000680109.1,ENST00000680183.1,ENST00000680434.1,ENST00000680860.1,ENST00000681137.1,ENST00000681286.1,ENST00000681338.1,ENST00000681791.1 | yes | — | no annotated structural change |
| PPIP5K1 | SCZ | cerebellar_hemisphere | chr15:43558932-43560413(-) | 0.926 | abf | intermediate | ENST00000420765.6,ENST00000439195.5,ENST00000644537.1 | no | — | — |
| PPIP5K1 | SCZ | cerebellum | chr15:43558932-43564103(-) | 0.869 | abf | primary_prior | ENST00000381879.8,ENST00000381885.5,ENST00000396923.7 | no | — | — |
| RNASEH2C | SCZ | anterior_cingulate_cortex_ba24 | chr11:65715317-65719681(-) | 0.981 | susie | robust | ENST00000644198.1 | no | — | — |
| RNASEH2C | SCZ | caudate_basal_ganglia | chr11:65715317-65719681(-) | 0.967 | susie | intermediate | ENST00000644198.1 | no | — | — |
| RNASEH2C | SCZ | putamen_basal_ganglia | chr11:65715317-65719681(-) | 0.961 | susie | intermediate | ENST00000644198.1 | no | — | — |
| RNASEH2C | SCZ | cerebellum | chr11:65719809-65720045(-) | 0.946 | susie | intermediate | ENST00000308418.10,ENST00000528220.2,ENST00000531596.6,ENST00000533698.5,ENST00000534482.6,ENST00000642430.1,ENST00000643214.1,ENST00000644142.1,ENST00000644198.1,ENST00000646597.1 | yes | — | no annotated structural change |
| RNASEH2C | SCZ | cerebellar_hemisphere | chr11:65719809-65720045(-) | 0.895 | susie | intermediate | ENST00000308418.10,ENST00000528220.2,ENST00000531596.6,ENST00000533698.5,ENST00000534482.6,ENST00000642430.1,ENST00000643214.1,ENST00000644142.1,ENST00000644198.1,ENST00000646597.1 | no | — | — |
| RNASEH2C | SCZ | nucleus_accumbens_basal_ganglia | chr11:65715317-65719681(-) | 0.866 | abf | primary_prior | ENST00000644198.1 | no | — | — |
| SYT5 | SCZ | frontal_cortex_ba9 | chr19:55179086-55179964(-) | 0.904 | abf | intermediate | ENST00000589172.5 | no | — | — |
| TMED4 | SCZ | hypothalamus | chr7:44579628-44581449(-) | 0.965 | susie | intermediate | ENST00000289577.10 | no | — | — |
| TMED4 | SCZ | cerebellar_hemisphere | chr7:44579628-44581093(-) | 0.963 | susie | intermediate | ENST00000457408.7 | no | — | — |
| TMED4 | SCZ | cerebellum | chr7:44579628-44581449(-) | 0.961 | susie | intermediate | ENST00000289577.10 | yes | — | no annotated structural change |
| TMED4 | SCZ | putamen_basal_ganglia | chr7:44579628-44581093(-) | 0.956 | susie | intermediate | ENST00000457408.7 | no | — | — |
| TMED4 | SCZ | nucleus_accumbens_basal_ganglia | chr7:44579628-44581449(-) | 0.955 | susie | intermediate | ENST00000289577.10 | no | — | — |
| TMED4 | SCZ | caudate_basal_ganglia | chr7:44579628-44581093(-) | 0.955 | susie | intermediate | ENST00000457408.7 | no | yes | — |
| TMED4 | SCZ | hippocampus | chr7:44579628-44581093(-) | 0.954 | susie | intermediate | ENST00000457408.7 | yes | yes | no annotated structural change |
| TMED4 | SCZ | amygdala | chr7:44579628-44581449(-) | 0.954 | susie | intermediate | ENST00000289577.10 | yes | — | no annotated structural change |
| TMED4 | SCZ | frontal_cortex_ba9 | chr7:44579628-44581449(-) | 0.954 | susie | intermediate | ENST00000289577.10 | yes | yes | no annotated structural change |
| TMED4 | SCZ | cortex | chr7:44579628-44581449(-) | 0.952 | susie | intermediate | ENST00000289577.10 | no | yes | — |
| TMED4 | SCZ | spinal_cord_cervical_c_1 | chr7:44579628-44581093(-) | 0.922 | abf | intermediate | ENST00000457408.7 | no | — | — |
| TMED4 | SCZ | anterior_cingulate_cortex_ba24 | chr7:44579628-44581093(-) | 0.914 | susie | intermediate | ENST00000457408.7 | no | — | — |
| TMED4 | SCZ | substantia_nigra | chr7:44579628-44581449(-) | 0.914 | abf | intermediate | ENST00000289577.10 | no | — | — |
