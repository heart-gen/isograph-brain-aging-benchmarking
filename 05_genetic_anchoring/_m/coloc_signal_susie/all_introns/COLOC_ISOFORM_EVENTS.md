# Resolved isoform events — signal-level colocalization (all-introns arm)

Every tissue where a `coloc.susie` nomination's sQTL reaches PP4 >= 0.8, the intron that call used, and the IsoGraph switch it lands on in the same GTEx tissue. Junction matching, switch evidence and the BrainSeq check are the same code as the CLPP layer (`coloc/coloc_isoform_events_combined.parquet`), which this table sits beside and does not replace.

Three differences from the CLPP layer, each limiting what a row can say:

- **No direction.** The risk-allele-signed effect is resolved only for the CLPP layer (`coloc_direction`); `risk_qtl_effect` is empty here rather than borrowed.
- **`estimator = abf` rows name GTEx's representative intron**, the only one the fallback scored. They resolve an event but never tested the gene's other introns.
- **One intron per tissue:** the intron behind that tissue's call. Other introns of the gene that also colocalize are in `signal_pairs.parquet`, not here.

- events (nomination x tissue): **410** over **119** gene x trait nominations
- signal-level (`susie`): **111**; abf fallback: **299**
- junction mapped to a GENCODE v47 transcript: **362/410**
- junction on an IsoGraph switch-pair isoform in the same tissue: **115** events, **56** gene x trait
- also a switch-pair isoform in the matching BrainSeq region: **29**
- cross-tissue exceptions (switch pair taken from another region; flagged, never counted as tissue-concordant, scored apart in the long-read arm): **2** events (UNC13A)

| gene | trait | tissue | junction | PP4 | estimator | prior | transcripts | switch pair | BrainSeq | structural event |
|---|---|---|---|---|---|---|---|---|---|---|
| AKT1 | AD | caudate_basal_ganglia | chr14:104770847-104771722(-) | 0.840 | abf | primary_prior | ENST00000554585.6 | no | — | — |
| BCKDK | AD | spinal_cord_cervical_c_1 | chr16:31115869-31117463(+) | 0.995 | susie | robust | — | — | — | — |
| BCKDK | AD | cerebellum | chr16:31115869-31117062(+) | 0.995 | susie | robust | — | — | — | — |
| BCKDK | AD | cerebellar_hemisphere | chr16:31115869-31117062(+) | 0.994 | susie | robust | — | — | — | — |
| BCKDK | AD | frontal_cortex_ba9 | chr16:31115869-31117463(+) | 0.993 | susie | robust | — | — | — | — |
| COG7 | AD | amygdala | chr16:23389086-23392380(-) | 0.837 | abf | primary_prior | ENST00000307149.10,ENST00000566364.1 | no | — | — |
| DOC2A | AD | amygdala | chr16:30007090-30007173(-) | 0.881 | abf | primary_prior | ENST00000350119.9,ENST00000561671.5,ENST00000564944.5,ENST00000564979.5,ENST00000566310.5,ENST00000616445.4 | no | — | — |
| IFNAR2 | AD | nucleus_accumbens_basal_ganglia | chr21:33245074-33246718(+) | 0.870 | abf | primary_prior | ENST00000342101.7,ENST00000342136.9,ENST00000382238.6,ENST00000382264.7,ENST00000404220.7,ENST00000413881.5,ENST00000443073.5,ENST00000447980.1,ENST00000683941.1,ENST00000700429.1 | no | — | — |
| IFNAR2 | AD | spinal_cord_cervical_c_1 | chr21:33245074-33246718(+) | 0.865 | abf | primary_prior | ENST00000342101.7,ENST00000342136.9,ENST00000382238.6,ENST00000382264.7,ENST00000404220.7,ENST00000413881.5,ENST00000443073.5,ENST00000447980.1,ENST00000683941.1,ENST00000700429.1 | no | — | — |
| IFNAR2 | AD | cerebellum | chr21:33260727-33262793(+) | 0.864 | abf | primary_prior | ENST00000342136.9,ENST00000382238.6,ENST00000443073.5,ENST00000683941.1 | yes | — | no annotated structural change |
| IFNAR2 | AD | cerebellar_hemisphere | chr21:33245074-33246718(+) | 0.861 | abf | primary_prior | ENST00000342101.7,ENST00000342136.9,ENST00000382238.6,ENST00000382264.7,ENST00000404220.7,ENST00000413881.5,ENST00000443073.5,ENST00000447980.1,ENST00000683941.1,ENST00000700429.1 | no | — | — |
| IFNAR2 | AD | caudate_basal_ganglia | chr21:33246482-33246718(+) | 0.855 | abf | primary_prior | ENST00000682044.1,ENST00000700427.1 | no | — | — |
| IFNAR2 | AD | putamen_basal_ganglia | chr21:33245074-33246718(+) | 0.843 | abf | primary_prior | ENST00000342101.7,ENST00000342136.9,ENST00000382238.6,ENST00000382264.7,ENST00000404220.7,ENST00000413881.5,ENST00000443073.5,ENST00000447980.1,ENST00000683941.1,ENST00000700429.1 | no | — | — |
| IFNAR2 | AD | hippocampus | chr21:33245074-33246718(+) | 0.840 | abf | primary_prior | ENST00000342101.7,ENST00000342136.9,ENST00000382238.6,ENST00000382264.7,ENST00000404220.7,ENST00000413881.5,ENST00000443073.5,ENST00000447980.1,ENST00000683941.1,ENST00000700429.1 | no | — | — |
| INO80E | AD | cerebellum | chr16:30004657-30005221(+) | 0.903 | susie | intermediate | ENST00000540562.1,ENST00000562441.5,ENST00000567987.5,ENST00000569957.5 | yes | — | no annotated structural change |
| INO80E | AD | nucleus_accumbens_basal_ganglia | chr16:30004657-30005221(+) | 0.832 | abf | primary_prior | ENST00000540562.1,ENST00000562441.5,ENST00000567987.5,ENST00000569957.5 | no | — | — |
| INO80E | AD | cerebellar_hemisphere | chr16:30001040-30001212(+) | 0.815 | susie | primary_prior | ENST00000540562.1,ENST00000562441.5,ENST00000567065.5,ENST00000569957.5,ENST00000620599.4 | yes | — | no annotated structural change |
| INTS8 | AD | nucleus_accumbens_basal_ganglia | chr8:94873477-94874552(+) | 0.940 | susie | intermediate | ENST00000343161.8,ENST00000519736.5,ENST00000521155.5,ENST00000523206.5,ENST00000523731.6,ENST00000524333.5,ENST00000715987.1 | no | — | — |
| INTS8 | AD | cerebellum | chr8:94873477-94874552(+) | 0.931 | susie | intermediate | ENST00000343161.8,ENST00000519736.5,ENST00000521155.5,ENST00000523206.5,ENST00000523731.6,ENST00000524333.5,ENST00000715987.1 | yes | — | no annotated structural change |
| INTS8 | AD | hippocampus | chr8:94874602-94876074(+) | 0.925 | abf | intermediate | ENST00000343161.8,ENST00000521155.5,ENST00000523206.5,ENST00000523731.6,ENST00000524333.5,ENST00000715987.1 | yes | — | no annotated structural change |
| INTS8 | AD | hypothalamus | chr8:94873477-94874552(+) | 0.828 | susie | primary_prior | ENST00000343161.8,ENST00000519736.5,ENST00000521155.5,ENST00000523206.5,ENST00000523731.6,ENST00000524333.5,ENST00000715987.1 | no | — | — |
| ITGB1BP1 | AD | cerebellar_hemisphere | chr2:9418732-9418818(-) | 0.979 | susie | robust | — | — | — | — |
| ITGB1BP1 | AD | cerebellum | chr2:9418732-9418818(-) | 0.977 | susie | robust | — | — | — | — |
| ITGB1BP1 | AD | cortex | chr2:9418732-9418818(-) | 0.955 | susie | intermediate | — | — | — | — |
| ITGB1BP1 | AD | frontal_cortex_ba9 | chr2:9418732-9419990(-) | 0.951 | susie | intermediate | ENST00000360635.7,ENST00000488451.5,ENST00000492079.5,ENST00000494563.5 | yes | — | no annotated structural change |
| ITGB1BP1 | AD | nucleus_accumbens_basal_ganglia | chr2:9418732-9419990(-) | 0.947 | susie | intermediate | ENST00000360635.7,ENST00000488451.5,ENST00000492079.5,ENST00000494563.5 | no | — | — |
| NDUFS3 | AD | cerebellar_hemisphere | chr11:47565986-47571471(+) | 0.864 | abf | primary_prior | ENST00000533507.5 | no | — | — |
| PICALM | AD | cortex | chr11:85974812-85981129(-) | 0.818 | abf | primary_prior | ENST00000356360.9,ENST00000532603.5 | no | — | — |
| PILRB | AD | putamen_basal_ganglia | chr7:100358366-100358933(+) | 0.999 | susie | robust | ENST00000438028.5 | no | — | — |
| RAD51C | AD | cortex | chr17:58692959-58694931(+) | 0.992 | susie | robust | ENST00000697692.1 | no | — | — |
| SERPINB1 | AD | spinal_cord_cervical_c_1 | chr6:2840594-2841216(-) | 0.917 | abf | intermediate | ENST00000476896.5 | no | — | — |
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
| SIRPA | AD | substantia_nigra | chr20:1915455-1921395(+) | 0.955 | abf | intermediate | ENST00000356025.7,ENST00000358771.5,ENST00000400068.7,ENST00000622179.4 | yes | — | no annotated structural change |
| SIRPA | AD | hippocampus | chr20:1915455-1921395(+) | 0.950 | abf | intermediate | ENST00000356025.7,ENST00000358771.5,ENST00000400068.7,ENST00000622179.4 | yes | — | no annotated structural change |
| SIRPA | AD | nucleus_accumbens_basal_ganglia | chr20:1915455-1921395(+) | 0.936 | abf | intermediate | ENST00000356025.7,ENST00000358771.5,ENST00000400068.7,ENST00000622179.4 | no | — | — |
| SLC39A13 | AD | spinal_cord_cervical_c_1 | chr11:47412467-47413597(+) | 0.987 | abf | robust | ENST00000531419.5,ENST00000531865.5 | yes | — | no annotated structural change |
| SLC39A13 | AD | hippocampus | chr11:47412467-47413597(+) | 0.983 | abf | robust | ENST00000531419.5,ENST00000531865.5 | no | — | — |
| SPAG9 | AD | hypothalamus | chr17:50974538-50974771(-) | 0.839 | abf | primary_prior | — | — | — | — |
| SPAG9 | AD | putamen_basal_ganglia | chr17:50974538-50974771(-) | 0.827 | abf | primary_prior | — | — | — | — |
| SPAG9 | AD | hippocampus | chr17:50974538-50974771(-) | 0.823 | abf | primary_prior | — | — | — | — |
| SPAG9 | AD | anterior_cingulate_cortex_ba24 | chr17:50974538-50974771(-) | 0.823 | abf | primary_prior | — | — | — | — |
| SPAG9 | AD | substantia_nigra | chr17:50974538-50974771(-) | 0.821 | abf | primary_prior | — | — | — | — |
| SPAG9 | AD | amygdala | chr17:50974538-50974771(-) | 0.811 | abf | primary_prior | — | — | — | — |
| SPAG9 | AD | caudate_basal_ganglia | chr17:50974538-50974771(-) | 0.809 | abf | primary_prior | — | — | — | — |
| SPAG9 | AD | cortex | chr17:50970856-50974771(-) | 0.802 | abf | primary_prior | ENST00000262013.12,ENST00000357122.8,ENST00000505279.5,ENST00000506500.1,ENST00000510283.5 | no | yes | — |
| SPI1 | AD | cerebellar_hemisphere | chr11:47383832-47408952(-) | 0.926 | abf | intermediate | ENST00000713543.1 | yes | — | no annotated structural change |
| TPCN1 | AD | cerebellar_hemisphere | chr12:113280195-113283693(+) | 0.861 | susie | primary_prior | — | — | — | — |
| TPCN1 | AD | cerebellum | chr12:113280195-113284581(+) | 0.836 | susie | primary_prior | ENST00000335509.11,ENST00000392569.8,ENST00000428632.7,ENST00000541517.5,ENST00000546781.5,ENST00000550785.5,ENST00000551127.5,ENST00000552077.5 | yes | — | no annotated structural change |
| VWA5B2 | AD | amygdala | chr3:184233409-184233576(+) | 0.864 | abf | primary_prior | — | — | — | — |
| YPEL3 | AD | substantia_nigra | chr16:30094897-30095252(-) | 0.935 | abf | intermediate | — | — | — | — |
| YPEL3 | AD | putamen_basal_ganglia | chr16:30094897-30095252(-) | 0.862 | abf | primary_prior | — | — | — | — |
| YPEL3 | AD | frontal_cortex_ba9 | chr16:30095335-30096097(-) | 0.845 | abf | primary_prior | ENST00000565479.5,ENST00000568674.1 | no | yes | — |
| YPEL3 | AD | caudate_basal_ganglia | chr16:30094897-30095252(-) | 0.822 | abf | primary_prior | — | — | — | — |
| ZNF232 | AD | spinal_cord_cervical_c_1 | chr17:5111219-5111800(-) | 0.899 | abf | intermediate | ENST00000570486.5,ENST00000571076.1,ENST00000573015.6,ENST00000575538.1,ENST00000696408.1,ENST00000696538.1 | no | — | — |
| C9orf72 | ALS | cerebellum | chr9:27556796-27558491(-) | 0.991 | abf | robust | ENST00000380003.8,ENST00000488117.5,ENST00000619707.5,ENST00000644136.1,ENST00000647196.1,ENST00000673600.1 | yes | — | no annotated structural change |
| C9orf72 | ALS | cerebellar_hemisphere | chr9:27567164-27573431(-) | 0.988 | abf | robust | ENST00000380003.8,ENST00000488117.5,ENST00000644136.1,ENST00000673600.1 | no | — | — |
| C9orf72 | ALS | hippocampus | chr9:27567164-27573431(-) | 0.984 | abf | robust | ENST00000380003.8,ENST00000488117.5,ENST00000644136.1,ENST00000673600.1 | no | — | — |
| C9orf72 | ALS | nucleus_accumbens_basal_ganglia | chr9:27567164-27573431(-) | 0.984 | abf | robust | ENST00000380003.8,ENST00000488117.5,ENST00000644136.1,ENST00000673600.1 | no | — | — |
| C9orf72 | ALS | caudate_basal_ganglia | chr9:27567164-27573431(-) | 0.966 | abf | intermediate | ENST00000380003.8,ENST00000488117.5,ENST00000644136.1,ENST00000673600.1 | no | — | — |
| C9orf72 | ALS | cortex | chr9:27567164-27573431(-) | 0.929 | abf | intermediate | ENST00000380003.8,ENST00000488117.5,ENST00000644136.1,ENST00000673600.1 | no | — | — |
| C9orf72 | ALS | anterior_cingulate_cortex_ba24 | chr9:27567164-27573431(-) | 0.910 | abf | intermediate | ENST00000380003.8,ENST00000488117.5,ENST00000644136.1,ENST00000673600.1 | no | — | — |
| C9orf72 | ALS | frontal_cortex_ba9 | chr9:27567164-27573431(-) | 0.831 | abf | primary_prior | ENST00000380003.8,ENST00000488117.5,ENST00000644136.1,ENST00000673600.1 | no | — | — |
| FNBP1 | ALS | hippocampus | chr9:129908999-129923844(-) | 0.929 | abf | intermediate | ENST00000449089.7,ENST00000703532.1 | no | — | — |
| FNBP1 | ALS | nucleus_accumbens_basal_ganglia | chr9:129957464-129958491(-) | 0.882 | abf | primary_prior | ENST00000355681.3,ENST00000446176.7,ENST00000449089.7,ENST00000699492.1,ENST00000703532.1,ENST00000703533.1,ENST00000703558.1,ENST00000703559.1,ENST00000703560.1 | no | — | — |
| FNBP1 | ALS | cortex | chr9:129957464-129958491(-) | 0.882 | abf | primary_prior | ENST00000355681.3,ENST00000446176.7,ENST00000449089.7,ENST00000699492.1,ENST00000703532.1,ENST00000703533.1,ENST00000703558.1,ENST00000703559.1,ENST00000703560.1 | yes | — | no annotated structural change |
| FNBP1 | ALS | anterior_cingulate_cortex_ba24 | chr9:129957464-129958491(-) | 0.844 | abf | primary_prior | ENST00000355681.3,ENST00000446176.7,ENST00000449089.7,ENST00000699492.1,ENST00000703532.1,ENST00000703533.1,ENST00000703558.1,ENST00000703559.1,ENST00000703560.1 | no | — | — |
| FNBP1 | ALS | hypothalamus | chr9:129957464-129958491(-) | 0.821 | abf | primary_prior | ENST00000355681.3,ENST00000446176.7,ENST00000449089.7,ENST00000699492.1,ENST00000703532.1,ENST00000703533.1,ENST00000703558.1,ENST00000703559.1,ENST00000703560.1 | yes | — | no annotated structural change |
| G2E3 | ALS | nucleus_accumbens_basal_ganglia | chr14:30581116-30586718(+) | 0.839 | abf | primary_prior | ENST00000206595.11,ENST00000547532.5,ENST00000550944.5,ENST00000552488.1,ENST00000553504.5,ENST00000554714.5,ENST00000555429.1,ENST00000648008.1 | no | — | — |
| GGNBP2 | ALS | cerebellum | chr17:36560871-36567663(+) | 0.912 | susie | intermediate | ENST00000613102.5,ENST00000617860.4,ENST00000618837.4 | no | — | — |
| NME4 | ALS | cortex | chr16:399739-400219(+) | 0.866 | abf | primary_prior | ENST00000219479.7,ENST00000382940.8,ENST00000397722.5,ENST00000444498.5,ENST00000448828.5,ENST00000450036.1,ENST00000460297.1,ENST00000468031.1,ENST00000620944.4,ENST00000621774.4 | no | — | — |
| NSMAF | ALS | cerebellum | chr8:58635546-58642984(-) | 0.940 | abf | intermediate | ENST00000038176.8,ENST00000427130.7,ENST00000519871.1,ENST00000522645.5 | no | — | — |
| NSMAF | ALS | cerebellar_hemisphere | chr8:58635546-58642984(-) | 0.867 | abf | primary_prior | ENST00000038176.8,ENST00000427130.7,ENST00000519871.1,ENST00000522645.5 | yes | — | no annotated structural change |
| PRDM2 | ALS | cortex | chr1:13816570-13823159(+) | 0.964 | abf | intermediate | ENST00000235372.11,ENST00000311066.10,ENST00000376048.9,ENST00000503842.5,ENST00000505823.5 | yes | yes | no annotated structural change |
| PRDM2 | ALS | substantia_nigra | chr1:13816570-13823159(+) | 0.947 | abf | intermediate | ENST00000235372.11,ENST00000311066.10,ENST00000376048.9,ENST00000503842.5,ENST00000505823.5 | no | — | — |
| PRDM2 | ALS | nucleus_accumbens_basal_ganglia | chr1:13816570-13823159(+) | 0.945 | abf | intermediate | ENST00000235372.11,ENST00000311066.10,ENST00000376048.9,ENST00000503842.5,ENST00000505823.5 | no | — | — |
| PRDM2 | ALS | hippocampus | chr1:13816570-13823159(+) | 0.900 | abf | intermediate | ENST00000235372.11,ENST00000311066.10,ENST00000376048.9,ENST00000503842.5,ENST00000505823.5 | no | — | — |
| PRDM2 | ALS | putamen_basal_ganglia | chr1:13816570-13823159(+) | 0.856 | abf | primary_prior | ENST00000235372.11,ENST00000311066.10,ENST00000376048.9,ENST00000503842.5,ENST00000505823.5 | no | — | — |
| PRDM2 | ALS | anterior_cingulate_cortex_ba24 | chr1:13816570-13823159(+) | 0.826 | abf | primary_prior | ENST00000235372.11,ENST00000311066.10,ENST00000376048.9,ENST00000503842.5,ENST00000505823.5 | yes | — | no annotated structural change |
| PTPRN | ALS | nucleus_accumbens_basal_ganglia | chr2:219295141-219295897(-) | 0.918 | abf | intermediate | — | — | — | — |
| PTPRN | ALS | cerebellum | chr2:219295141-219296226(-) | 0.881 | susie | primary_prior | ENST00000295718.7,ENST00000409251.7,ENST00000423636.6,ENST00000443981.5 | yes | — | no annotated structural change |
| PTPRN | ALS | cerebellar_hemisphere | chr2:219295966-219296226(-) | 0.880 | susie | primary_prior | — | — | — | — |
| PTPRN | ALS | cortex | chr2:219295141-219296226(-) | 0.879 | susie | primary_prior | ENST00000295718.7,ENST00000409251.7,ENST00000423636.6,ENST00000443981.5 | yes | — | no annotated structural change |
| PTPRN | ALS | frontal_cortex_ba9 | chr2:219295141-219296226(-) | 0.832 | abf | primary_prior | ENST00000295718.7,ENST00000409251.7,ENST00000423636.6,ENST00000443981.5 | no | — | — |
| RCSD1 | ALS | spinal_cord_cervical_c_1 | chr1:167685510-167690049(+) | 0.873 | abf | primary_prior | ENST00000367854.8 | no | — | — |
| RPS6KL1 | ALS | cerebellum | chr14:74921561-74923180(-) | 0.912 | abf | intermediate | ENST00000555834.5 | yes | — | no annotated structural change |
| RPS6KL1 | ALS | cerebellar_hemisphere | chr14:74921561-74923180(-) | 0.858 | abf | primary_prior | ENST00000555834.5 | no | — | — |
| RPS6KL1 | ALS | nucleus_accumbens_basal_ganglia | chr14:74921561-74923180(-) | 0.829 | abf | primary_prior | ENST00000555834.5 | no | — | — |
| SCFD1 | ALS | hypothalamus | chr14:30622399-30630477(+) | 0.856 | abf | primary_prior | ENST00000484733.6,ENST00000544052.6,ENST00000553693.6,ENST00000554437.6,ENST00000676465.1,ENST00000676473.1,ENST00000676520.1,ENST00000676674.1,ENST00000676834.1,ENST00000676914.1,ENST00000676954.1,ENST00000677413.1,ENST00000677456.1,ENST00000677690.1,ENST00000678579.1,ENST00000678637.1,ENST00000678858.1,ENST00000679342.1 | yes | — | no annotated structural change |
| TMEM175 | ALS | cerebellar_hemisphere | chr4:948615-950424(+) | 0.999 | susie | robust | ENST00000504180.5,ENST00000504744.5 | no | — | — |
| TMEM175 | ALS | cortex | chr4:932540-947709(+) | 0.998 | susie | robust | ENST00000264771.9,ENST00000438836.6,ENST00000452360.6,ENST00000504505.1,ENST00000504744.5,ENST00000505734.1,ENST00000507319.5,ENST00000513682.5,ENST00000513952.5,ENST00000514453.5,ENST00000515876.5,ENST00000622959.3 | yes | — | no annotated structural change |
| TMEM175 | ALS | cerebellum | chr4:947892-948116(+) | 0.997 | susie | robust | ENST00000264771.9,ENST00000438836.6,ENST00000452360.6,ENST00000504744.5,ENST00000507319.5,ENST00000513682.5,ENST00000513952.5,ENST00000514546.5,ENST00000515876.5,ENST00000622959.3 | yes | — | no annotated structural change |
| TNFSF13 | ALS | putamen_basal_ganglia | chr17:7559297-7559846(+) | 0.911 | abf | intermediate | ENST00000380535.8 | no | — | — |
| TNFSF13 | ALS | anterior_cingulate_cortex_ba24 | chr17:7559297-7559846(+) | 0.906 | abf | intermediate | ENST00000380535.8 | no | — | — |
| TPP1 | ALS | cerebellum | chr11:6611818-6611958(-) | 0.996 | abf | robust | — | — | — | — |
| TXNDC15 | ALS | caudate_basal_ganglia | chr5:134874530-134893492(+) | 0.969 | abf | intermediate | ENST00000511070.5 | no | — | — |
| UNC13A | ALS | cerebellum | chr19:17630750-17632782(-) | 0.959 | susie | intermediate | ENST00000519716.7,ENST00000550896.1,ENST00000551649.5,ENST00000552293.5 | cross-tissue exception (amygdala;anterior_cingulate_cortex_ba24;cortex;frontal_cortex_ba9;hippocampus;substantia_nigra) | — | — |
| UNC13A | ALS | cerebellar_hemisphere | chr19:17630750-17632782(-) | 0.936 | susie | intermediate | ENST00000519716.7,ENST00000550896.1,ENST00000551649.5,ENST00000552293.5 | cross-tissue exception (amygdala;anterior_cingulate_cortex_ba24;cortex;frontal_cortex_ba9;hippocampus;substantia_nigra) | — | — |
| WHAMM | ALS | hippocampus | chr15:82823287-82826410(+) | 0.916 | abf | intermediate | ENST00000286760.5 | no | — | — |
| WIPI2 | ALS | nucleus_accumbens_basal_ganglia | chr7:5229705-5230835(+) | 0.893 | abf | intermediate | ENST00000382384.6,ENST00000404704.7,ENST00000484262.1 | no | — | — |
| SNCA | LBD | cortex | chr4:89835692-89836127(-) | 0.974 | susie | intermediate | ENST00000508895.5,ENST00000618500.4 | no | — | — |
| SNCA | LBD | substantia_nigra | chr4:89835692-89836127(-) | 0.971 | susie | intermediate | ENST00000508895.5,ENST00000618500.4 | no | — | — |
| SNCA | LBD | hypothalamus | chr4:89835692-89836127(-) | 0.970 | susie | intermediate | ENST00000508895.5,ENST00000618500.4 | yes | — | no annotated structural change |
| SNCA | LBD | spinal_cord_cervical_c_1 | chr4:89835692-89836127(-) | 0.962 | susie | intermediate | ENST00000508895.5,ENST00000618500.4 | no | — | — |
| SNCA | LBD | frontal_cortex_ba9 | chr4:89835692-89836127(-) | 0.956 | susie | intermediate | ENST00000508895.5,ENST00000618500.4 | no | — | — |
| SNCA | LBD | nucleus_accumbens_basal_ganglia | chr4:89835692-89836127(-) | 0.956 | susie | intermediate | ENST00000508895.5,ENST00000618500.4 | no | — | — |
| SNCA | LBD | amygdala | chr4:89835692-89836127(-) | 0.955 | susie | intermediate | ENST00000508895.5,ENST00000618500.4 | no | — | — |
| SNCA | LBD | hippocampus | chr4:89835692-89836127(-) | 0.954 | susie | intermediate | ENST00000508895.5,ENST00000618500.4 | no | yes | — |
| SNCA | LBD | cerebellar_hemisphere | chr4:89835692-89836127(-) | 0.952 | susie | intermediate | ENST00000508895.5,ENST00000618500.4 | no | — | — |
| SNCA | LBD | anterior_cingulate_cortex_ba24 | chr4:89835692-89836127(-) | 0.951 | susie | intermediate | ENST00000508895.5,ENST00000618500.4 | no | — | — |
| SNCA | LBD | cerebellum | chr4:89835692-89836127(-) | 0.940 | susie | intermediate | ENST00000508895.5,ENST00000618500.4 | yes | — | no annotated structural change |
| AZI2 | PD | cerebellar_hemisphere | chr3:28331954-28332369(-) | 0.803 | abf | primary_prior | ENST00000420543.6,ENST00000457172.5,ENST00000488978.1 | no | — | — |
| CCDC62 | PD | cerebellar_hemisphere | chr12:122798200-122801124(+) | 0.946 | abf | intermediate | ENST00000253079.11,ENST00000341952.8,ENST00000392440.3,ENST00000392441.8,ENST00000537566.5 | yes | — | no annotated structural change |
| CDHR3 | PD | cerebellar_hemisphere | chr7:106018072-106020373(+) | 0.995 | abf | robust | ENST00000317716.14,ENST00000478080.5 | yes | — | no annotated structural change |
| CDHR3 | PD | cerebellum | chr7:106018072-106020373(+) | 0.922 | abf | intermediate | ENST00000317716.14,ENST00000478080.5 | yes | — | no annotated structural change |
| CTSB | PD | amygdala | chr8:11844241-11845082(-) | 0.868 | susie | primary_prior | ENST00000530640.7,ENST00000531089.6 | no | — | — |
| DDRGK1 | PD | anterior_cingulate_cortex_ba24 | chr20:3195353-3200001(-) | 0.884 | abf | primary_prior | ENST00000354488.8,ENST00000380201.2 | no | — | — |
| NCOR1 | PD | cerebellar_hemisphere | chr17:16171995-16186554(-) | 0.896 | susie | intermediate | ENST00000268712.8,ENST00000395851.5,ENST00000411510.5,ENST00000430577.2,ENST00000436068.2,ENST00000436828.5,ENST00000582357.5,ENST00000585296.1,ENST00000704744.1,ENST00000704745.1 | yes | — | no annotated structural change |
| NCOR1 | PD | nucleus_accumbens_basal_ganglia | chr17:16171995-16186554(-) | 0.805 | abf | primary_prior | ENST00000268712.8,ENST00000395851.5,ENST00000411510.5,ENST00000430577.2,ENST00000436068.2,ENST00000436828.5,ENST00000582357.5,ENST00000585296.1,ENST00000704744.1,ENST00000704745.1 | no | — | — |
| SH3GL2 | PD | substantia_nigra | chr9:17579287-17747066(+) | 0.834 | abf | primary_prior | ENST00000380607.5 | yes | — | no annotated structural change |
| SH3GL2 | PD | cortex | chr9:17579287-17747066(+) | 0.834 | abf | primary_prior | ENST00000380607.5 | yes | yes | no annotated structural change |
| SH3GL2 | PD | frontal_cortex_ba9 | chr9:17579287-17747066(+) | 0.834 | abf | primary_prior | ENST00000380607.5 | yes | yes | no annotated structural change |
| SH3GL2 | PD | anterior_cingulate_cortex_ba24 | chr9:17579287-17747066(+) | 0.834 | abf | primary_prior | ENST00000380607.5 | yes | — | no annotated structural change |
| SH3GL2 | PD | caudate_basal_ganglia | chr9:17579287-17747066(+) | 0.834 | abf | primary_prior | ENST00000380607.5 | no | — | — |
| SH3GL2 | PD | nucleus_accumbens_basal_ganglia | chr9:17579287-17747066(+) | 0.834 | abf | primary_prior | ENST00000380607.5 | no | — | — |
| SH3GL2 | PD | hypothalamus | chr9:17579287-17747066(+) | 0.834 | abf | primary_prior | ENST00000380607.5 | no | — | — |
| SH3GL2 | PD | spinal_cord_cervical_c_1 | chr9:17579287-17747066(+) | 0.834 | abf | primary_prior | ENST00000380607.5 | no | — | — |
| SH3GL2 | PD | amygdala | chr9:17579287-17747066(+) | 0.834 | abf | primary_prior | ENST00000380607.5 | yes | — | no annotated structural change |
| TTC19 | PD | putamen_basal_ganglia | chr17:16000245-16001915(+) | 0.899 | susie | intermediate | ENST00000261647.10,ENST00000470399.1,ENST00000497842.6 | no | — | — |
| TTC19 | PD | caudate_basal_ganglia | chr17:16006568-16025017(+) | 0.804 | susie | primary_prior | ENST00000261647.10,ENST00000475723.5 | no | yes | — |
| ZSWIM7 | PD | nucleus_accumbens_basal_ganglia | chr17:15977092-15978044(-) | 0.906 | susie | intermediate | ENST00000476496.5 | no | — | — |
| ZSWIM7 | PD | putamen_basal_ganglia | chr17:15977092-15978044(-) | 0.904 | susie | intermediate | ENST00000476496.5 | no | — | — |
| ZSWIM7 | PD | caudate_basal_ganglia | chr17:15993778-15999660(-) | 0.902 | susie | intermediate | ENST00000399280.6,ENST00000495825.6,ENST00000497719.5 | no | — | — |
| ZSWIM7 | PD | cortex | chr17:15993778-15999519(-) | 0.885 | susie | primary_prior | ENST00000399277.6,ENST00000472495.5,ENST00000486655.5 | yes | — | no annotated structural change |
| ZSWIM7 | PD | substantia_nigra | chr17:15977092-15978044(-) | 0.878 | susie | primary_prior | ENST00000476496.5 | no | — | — |
| ZSWIM7 | PD | frontal_cortex_ba9 | chr17:15977092-15978044(-) | 0.846 | susie | primary_prior | ENST00000476496.5 | no | — | — |
| ZSWIM7 | PD | hippocampus | chr17:15977092-15978044(-) | 0.844 | susie | primary_prior | ENST00000476496.5 | no | — | — |
| ZSWIM7 | PD | cerebellar_hemisphere | chr17:15977092-15978044(-) | 0.827 | susie | primary_prior | ENST00000476496.5 | no | — | — |
| ZSWIM7 | PD | spinal_cord_cervical_c_1 | chr17:15978163-15981040(-) | 0.823 | susie | primary_prior | ENST00000399277.6,ENST00000399280.6,ENST00000460252.5,ENST00000460315.5,ENST00000472495.5,ENST00000474716.5,ENST00000475498.5,ENST00000476496.5,ENST00000486706.6,ENST00000490395.5,ENST00000491631.5,ENST00000495825.6,ENST00000497719.5,ENST00000579955.1,ENST00000584519.5,ENST00000585208.5 | no | — | — |
| ZSWIM7 | PD | anterior_cingulate_cortex_ba24 | chr17:15977092-15978044(-) | 0.822 | susie | primary_prior | ENST00000476496.5 | no | — | — |
| ZSWIM7 | PD | amygdala | chr17:15977092-15978044(-) | 0.816 | susie | primary_prior | ENST00000476496.5 | no | — | — |
| ZSWIM7 | PD | hypothalamus | chr17:15977092-15978044(-) | 0.803 | susie | primary_prior | ENST00000476496.5 | no | — | — |
| ACTR1B | SCZ | caudate_basal_ganglia | chr2:97658316-97658427(-) | 0.997 | abf | robust | ENST00000289228.7,ENST00000451664.1 | no | yes | — |
| ACTR1B | SCZ | caudate_basal_ganglia | chr2:97658316-97658427(-) | 0.997 | abf | robust | ENST00000289228.7,ENST00000451664.1 | no | yes | — |
| ACTR1B | SCZ | hypothalamus | chr2:97658316-97658427(-) | 0.996 | abf | robust | ENST00000289228.7,ENST00000451664.1 | no | — | — |
| ACTR1B | SCZ | hypothalamus | chr2:97658316-97658427(-) | 0.996 | abf | robust | ENST00000289228.7,ENST00000451664.1 | no | — | — |
| ACTR1B | SCZ | frontal_cortex_ba9 | chr2:97658316-97658427(-) | 0.996 | abf | robust | ENST00000289228.7,ENST00000451664.1 | yes | — | no annotated structural change |
| ACTR1B | SCZ | frontal_cortex_ba9 | chr2:97658316-97658427(-) | 0.996 | abf | robust | ENST00000289228.7,ENST00000451664.1 | yes | — | no annotated structural change |
| ACTR1B | SCZ | nucleus_accumbens_basal_ganglia | chr2:97658316-97658427(-) | 0.995 | abf | robust | ENST00000289228.7,ENST00000451664.1 | no | — | — |
| ACTR1B | SCZ | nucleus_accumbens_basal_ganglia | chr2:97658316-97658427(-) | 0.995 | abf | robust | ENST00000289228.7,ENST00000451664.1 | no | — | — |
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
| ASB3 | SCZ | cerebellar_hemisphere | chr2:53729570-53750783(-) | 0.859 | abf | primary_prior | ENST00000263634.8,ENST00000394717.3,ENST00000406053.5,ENST00000406625.6,ENST00000406687.5,ENST00000482829.5,ENST00000489508.5 | yes | — | no annotated structural change |
| CATSPER2 | SCZ | hippocampus | chr15:43632363-43632717(-) | 0.894 | abf | intermediate | ENST00000321596.6,ENST00000381761.6,ENST00000396879.8,ENST00000419262.1,ENST00000450810.5,ENST00000472960.1 | no | yes | — |
| CATSPER2 | SCZ | hippocampus | chr15:43632363-43632717(-) | 0.894 | abf | intermediate | ENST00000321596.6,ENST00000381761.6,ENST00000396879.8,ENST00000419262.1,ENST00000450810.5,ENST00000472960.1 | no | yes | — |
| CCDC122 | SCZ | spinal_cord_cervical_c_1 | chr13:43874927-43879418(-) | 0.944 | abf | intermediate | — | — | — | — |
| CCDC122 | SCZ | putamen_basal_ganglia | chr13:43874927-43879418(-) | 0.932 | abf | intermediate | — | — | — | — |
| CCDC122 | SCZ | hippocampus | chr13:43874927-43879418(-) | 0.913 | abf | intermediate | — | — | — | — |
| CCDC122 | SCZ | hypothalamus | chr13:43874927-43879418(-) | 0.856 | abf | primary_prior | — | — | — | — |
| CCDC122 | SCZ | substantia_nigra | chr13:43874927-43879418(-) | 0.834 | abf | primary_prior | — | — | — | — |
| CCS | SCZ | anterior_cingulate_cortex_ba24 | chr11:66593714-66599116(+) | 0.937 | abf | intermediate | ENST00000310190.8,ENST00000526066.5,ENST00000530384.5,ENST00000530961.5,ENST00000531990.1,ENST00000533244.6 | no | — | — |
| CD46 | SCZ | cerebellum | chr1:207767195-207770321(+) | 0.804 | abf | primary_prior | ENST00000322918.9,ENST00000357714.5,ENST00000367041.5,ENST00000695778.1,ENST00000695780.1,ENST00000695781.1,ENST00000695782.1 | yes | — | no annotated structural change |
| CDIP1 | SCZ | caudate_basal_ganglia | chr16:4514664-4538325(-) | 0.951 | susie | intermediate | ENST00000563332.6,ENST00000588381.1 | no | yes | — |
| CDIP1 | SCZ | hypothalamus | chr16:4514664-4538325(-) | 0.935 | susie | intermediate | ENST00000563332.6,ENST00000588381.1 | no | — | — |
| CDIP1 | SCZ | amygdala | chr16:4514664-4538702(-) | 0.929 | susie | intermediate | ENST00000562334.5,ENST00000562579.5,ENST00000563507.5,ENST00000564828.5,ENST00000567695.6 | yes | — | no annotated structural change |
| CDIP1 | SCZ | substantia_nigra | chr16:4514664-4538702(-) | 0.928 | abf | intermediate | ENST00000562334.5,ENST00000562579.5,ENST00000563507.5,ENST00000564828.5,ENST00000567695.6 | yes | — | no annotated structural change |
| CDIP1 | SCZ | cortex | chr16:4514664-4538702(-) | 0.832 | susie | primary_prior | ENST00000562334.5,ENST00000562579.5,ENST00000563507.5,ENST00000564828.5,ENST00000567695.6 | yes | — | no annotated structural change |
| CDIP1 | SCZ | hippocampus | chr16:4514664-4538702(-) | 0.801 | susie | primary_prior | ENST00000562334.5,ENST00000562579.5,ENST00000563507.5,ENST00000564828.5,ENST00000567695.6 | yes | — | no annotated structural change |
| COPA | SCZ | cerebellum | chr1:160332557-160335242(-) | 0.830 | abf | primary_prior | ENST00000647799.1,ENST00000696207.1 | no | — | — |
| CRELD2 | SCZ | cerebellum | chr22:49921761-49922612(+) | 0.949 | susie | intermediate | ENST00000328268.9,ENST00000403427.3,ENST00000407217.7,ENST00000444954.1,ENST00000450207.5,ENST00000462253.5,ENST00000482956.5,ENST00000483652.5 | yes | — | no annotated structural change |
| CRELD2 | SCZ | cortex | chr22:49921761-49922612(+) | 0.947 | susie | intermediate | ENST00000328268.9,ENST00000403427.3,ENST00000407217.7,ENST00000444954.1,ENST00000450207.5,ENST00000462253.5,ENST00000482956.5,ENST00000483652.5 | no | — | — |
| CRELD2 | SCZ | putamen_basal_ganglia | chr22:49921761-49922612(+) | 0.945 | susie | intermediate | ENST00000328268.9,ENST00000403427.3,ENST00000407217.7,ENST00000444954.1,ENST00000450207.5,ENST00000462253.5,ENST00000482956.5,ENST00000483652.5 | no | — | — |
| CRELD2 | SCZ | cerebellar_hemisphere | chr22:49922443-49922612(+) | 0.943 | susie | intermediate | ENST00000404488.7 | yes | — | no annotated structural change |
| CRELD2 | SCZ | frontal_cortex_ba9 | chr22:49922443-49922612(+) | 0.943 | susie | intermediate | ENST00000404488.7 | no | — | — |
| CRELD2 | SCZ | nucleus_accumbens_basal_ganglia | chr22:49922443-49922612(+) | 0.939 | susie | intermediate | ENST00000404488.7 | no | — | — |
| CRELD2 | SCZ | anterior_cingulate_cortex_ba24 | chr22:49921761-49922612(+) | 0.934 | susie | intermediate | ENST00000328268.9,ENST00000403427.3,ENST00000407217.7,ENST00000444954.1,ENST00000450207.5,ENST00000462253.5,ENST00000482956.5,ENST00000483652.5 | no | — | — |
| CRELD2 | SCZ | caudate_basal_ganglia | chr22:49921761-49922294(+) | 0.933 | susie | intermediate | ENST00000404488.7 | no | — | — |
| CRELD2 | SCZ | hypothalamus | chr22:49921761-49922612(+) | 0.932 | susie | intermediate | ENST00000328268.9,ENST00000403427.3,ENST00000407217.7,ENST00000444954.1,ENST00000450207.5,ENST00000462253.5,ENST00000482956.5,ENST00000483652.5 | no | — | — |
| CRELD2 | SCZ | substantia_nigra | chr22:49921761-49922612(+) | 0.927 | susie | intermediate | ENST00000328268.9,ENST00000403427.3,ENST00000407217.7,ENST00000444954.1,ENST00000450207.5,ENST00000462253.5,ENST00000482956.5,ENST00000483652.5 | yes | — | no annotated structural change |
| CRELD2 | SCZ | spinal_cord_cervical_c_1 | chr22:49921761-49922612(+) | 0.919 | susie | intermediate | ENST00000328268.9,ENST00000403427.3,ENST00000407217.7,ENST00000444954.1,ENST00000450207.5,ENST00000462253.5,ENST00000482956.5,ENST00000483652.5 | no | — | — |
| CRELD2 | SCZ | hippocampus | chr22:49922443-49922612(+) | 0.919 | susie | intermediate | ENST00000404488.7 | yes | — | no annotated structural change |
| CRELD2 | SCZ | amygdala | chr22:49921761-49922294(+) | 0.918 | susie | intermediate | ENST00000404488.7 | no | — | — |
| DGKZ | SCZ | frontal_cortex_ba9 | chr11:46369550-46371313(+) | 0.949 | abf | intermediate | ENST00000318201.12,ENST00000531879.5 | no | yes | — |
| DGKZ | SCZ | cortex | chr11:46369550-46371313(+) | 0.939 | abf | intermediate | ENST00000318201.12,ENST00000531879.5 | no | yes | — |
| DGKZ | SCZ | nucleus_accumbens_basal_ganglia | chr11:46369550-46371313(+) | 0.928 | abf | intermediate | ENST00000318201.12,ENST00000531879.5 | no | — | — |
| DGKZ | SCZ | cerebellar_hemisphere | chr11:46369550-46371313(+) | 0.924 | abf | intermediate | ENST00000318201.12,ENST00000531879.5 | no | — | — |
| DGKZ | SCZ | caudate_basal_ganglia | chr11:46369550-46371313(+) | 0.917 | abf | intermediate | ENST00000318201.12,ENST00000531879.5 | no | yes | — |
| DGKZ | SCZ | spinal_cord_cervical_c_1 | chr11:46369550-46371313(+) | 0.892 | abf | intermediate | ENST00000318201.12,ENST00000531879.5 | no | — | — |
| DGKZ | SCZ | amygdala | chr11:46369550-46371313(+) | 0.878 | abf | primary_prior | ENST00000318201.12,ENST00000531879.5 | yes | — | no annotated structural change |
| DGKZ | SCZ | anterior_cingulate_cortex_ba24 | chr11:46369550-46371313(+) | 0.849 | abf | primary_prior | ENST00000318201.12,ENST00000531879.5 | yes | — | no annotated structural change |
| DGKZ | SCZ | hippocampus | chr11:46369550-46371313(+) | 0.838 | abf | primary_prior | ENST00000318201.12,ENST00000531879.5 | yes | yes | no annotated structural change |
| DNAJA3 | SCZ | cerebellar_hemisphere | chr16:4434517-4437402(+) | 0.963 | susie | intermediate | ENST00000262375.11,ENST00000355296.8,ENST00000572009.1,ENST00000572139.5,ENST00000573120.5,ENST00000574895.1,ENST00000575106.5 | yes | — | no annotated structural change |
| DNAJA3 | SCZ | cerebellum | chr16:4434517-4437402(+) | 0.902 | susie | intermediate | ENST00000262375.11,ENST00000355296.8,ENST00000572009.1,ENST00000572139.5,ENST00000573120.5,ENST00000574895.1,ENST00000575106.5 | yes | — | no annotated structural change |
| DNAJA3 | SCZ | hypothalamus | chr16:4450497-4455546(+) | 0.900 | abf | intermediate | ENST00000355296.8,ENST00000431375.6,ENST00000576911.5 | yes | — | no annotated structural change |
| DOC2A | SCZ | frontal_cortex_ba9 | chr16:30007090-30007179(-) | 0.970 | susie | intermediate | — | — | — | — |
| DOC2A | SCZ | cerebellar_hemisphere | chr16:30007299-30008996(-) | 0.969 | susie | intermediate | ENST00000350119.9,ENST00000563378.5,ENST00000564944.5,ENST00000564979.5,ENST00000565273.5,ENST00000566310.5,ENST00000567332.6,ENST00000616445.4 | yes | — | no annotated structural change |
| DOC2A | SCZ | cortex | chr16:30007090-30007179(-) | 0.967 | susie | intermediate | — | — | — | — |
| DOC2A | SCZ | hippocampus | chr16:30007090-30007173(-) | 0.964 | abf | intermediate | ENST00000350119.9,ENST00000561671.5,ENST00000564944.5,ENST00000564979.5,ENST00000566310.5,ENST00000616445.4 | no | — | — |
| DOC2A | SCZ | anterior_cingulate_cortex_ba24 | chr16:30007090-30007173(-) | 0.959 | abf | intermediate | ENST00000350119.9,ENST00000561671.5,ENST00000564944.5,ENST00000564979.5,ENST00000566310.5,ENST00000616445.4 | no | — | — |
| DOC2A | SCZ | hypothalamus | chr16:30007090-30007173(-) | 0.949 | abf | intermediate | ENST00000350119.9,ENST00000561671.5,ENST00000564944.5,ENST00000564979.5,ENST00000566310.5,ENST00000616445.4 | yes | — | no annotated structural change |
| DOC2A | SCZ | amygdala | chr16:30007090-30007173(-) | 0.947 | abf | intermediate | ENST00000350119.9,ENST00000561671.5,ENST00000564944.5,ENST00000564979.5,ENST00000566310.5,ENST00000616445.4 | no | — | — |
| DOC2A | SCZ | nucleus_accumbens_basal_ganglia | chr16:30007090-30007173(-) | 0.901 | abf | intermediate | ENST00000350119.9,ENST00000561671.5,ENST00000564944.5,ENST00000564979.5,ENST00000566310.5,ENST00000616445.4 | no | — | — |
| DOC2A | SCZ | caudate_basal_ganglia | chr16:30007090-30007173(-) | 0.871 | abf | primary_prior | ENST00000350119.9,ENST00000561671.5,ENST00000564944.5,ENST00000564979.5,ENST00000566310.5,ENST00000616445.4 | no | — | — |
| EFHB | SCZ | frontal_cortex_ba9 | chr3:19946254-19946919(-) | 0.965 | abf | intermediate | ENST00000474780.1 | no | — | — |
| FAM120AOS | SCZ | cerebellum | chr9:93437866-93438627(-) | 0.919 | abf | intermediate | — | — | — | — |
| FAM120AOS | SCZ | cerebellar_hemisphere | chr9:93435922-93437022(-) | 0.889 | abf | primary_prior | ENST00000445280.1 | no | — | — |
| FAM184A | SCZ | cerebellar_hemisphere | chr6:118966952-118975024(-) | 0.832 | abf | primary_prior | ENST00000352896.9,ENST00000368475.8,ENST00000521531.5 | no | — | — |
| FANCI | SCZ | nucleus_accumbens_basal_ganglia | chr15:89307672-89312904(+) | 0.879 | abf | primary_prior | ENST00000300027.12,ENST00000310775.12,ENST00000447611.6,ENST00000561894.1,ENST00000566895.5,ENST00000676003.1,ENST00000696717.1,ENST00000696718.1,ENST00000696719.1,ENST00000696721.1 | no | — | — |
| FGFR1 | SCZ | caudate_basal_ganglia | chr8:38429390-38429682(-) | 0.889 | abf | intermediate | ENST00000484370.5 | no | — | — |
| FGFR1 | SCZ | cerebellum | chr8:38428435-38429682(-) | 0.841 | abf | primary_prior | ENST00000335922.9,ENST00000341462.9,ENST00000397091.9,ENST00000397108.8,ENST00000397113.6,ENST00000425967.8,ENST00000447712.7,ENST00000474970.1,ENST00000525001.5,ENST00000532386.5,ENST00000532791.5,ENST00000649678.1,ENST00000674189.1,ENST00000674380.1,ENST00000674474.1,ENST00000683276.1,ENST00000683765.1,ENST00000683795.1,ENST00000683815.1,ENST00000703405.1 | yes | — | no annotated structural change |
| FOXN2 | SCZ | cerebellar_hemisphere | chr2:48346751-48359047(+) | 0.997 | abf | robust | ENST00000340553.8,ENST00000413569.5 | no | — | — |
| FOXN2 | SCZ | cerebellum | chr2:48346751-48359047(+) | 0.996 | abf | robust | ENST00000340553.8,ENST00000413569.5 | no | — | — |
| GABBR2 | SCZ | cerebellum | chr9:98299353-98303241(-) | 0.975 | susie | intermediate | ENST00000259455.4,ENST00000637410.1 | no | — | — |
| GALNT15 | SCZ | cerebellar_hemisphere | chr3:16222758-16227354(+) | 0.850 | abf | primary_prior | ENST00000339732.10 | no | — | — |
| GLYCTK | SCZ | anterior_cingulate_cortex_ba24 | chr3:52293192-52293307(+) | 0.812 | abf | primary_prior | ENST00000471180.5,ENST00000473032.5 | no | — | — |
| GLYCTK | SCZ | anterior_cingulate_cortex_ba24 | chr3:52293192-52293307(+) | 0.812 | abf | primary_prior | ENST00000471180.5,ENST00000473032.5 | no | — | — |
| GPM6A | SCZ | caudate_basal_ganglia | chr4:175701767-176002309(-) | 0.993 | susie | robust | ENST00000506894.5 | no | yes | — |
| GPM6A | SCZ | anterior_cingulate_cortex_ba24 | chr4:175701767-176002309(-) | 0.992 | susie | robust | ENST00000506894.5 | no | — | — |
| GPM6A | SCZ | putamen_basal_ganglia | chr4:175701767-175812191(-) | 0.908 | abf | intermediate | ENST00000280187.11,ENST00000393658.7,ENST00000513365.1 | yes | — | no annotated structural change |
| GPM6A | SCZ | amygdala | chr4:175701767-175812191(-) | 0.841 | abf | primary_prior | ENST00000280187.11,ENST00000393658.7,ENST00000513365.1 | no | — | — |
| GPR135 | SCZ | cerebellum | chr14:59459004-59460972(-) | 0.933 | susie | intermediate | ENST00000481661.1 | yes | — | no annotated structural change |
| GPR135 | SCZ | cerebellar_hemisphere | chr14:59473090-59474439(-) | 0.933 | susie | intermediate | — | — | — | — |
| HMOX2 | SCZ | frontal_cortex_ba9 | chr16:4483754-4505484(+) | 0.829 | susie | primary_prior | ENST00000458134.7,ENST00000570445.5 | yes | — | no annotated structural change |
| IDH3B | SCZ | cerebellum | chr20:2658837-2659109(-) | 0.860 | abf | primary_prior | — | — | — | — |
| IDH3B | SCZ | cerebellum | chr20:2658837-2659109(-) | 0.860 | abf | primary_prior | — | — | — | — |
| IDH3B | SCZ | cortex | chr20:2658837-2659109(-) | 0.859 | abf | primary_prior | — | — | — | — |
| IDH3B | SCZ | cortex | chr20:2658837-2659109(-) | 0.859 | abf | primary_prior | — | — | — | — |
| IDH3B | SCZ | cerebellar_hemisphere | chr20:2658837-2659109(-) | 0.840 | abf | primary_prior | — | — | — | — |
| IDH3B | SCZ | cerebellar_hemisphere | chr20:2658837-2659109(-) | 0.840 | abf | primary_prior | — | — | — | — |
| IDH3B | SCZ | hippocampus | chr20:2659312-2659525(-) | 0.826 | abf | primary_prior | ENST00000613370.1 | yes | — | no annotated structural change |
| IDH3B | SCZ | hippocampus | chr20:2659312-2659525(-) | 0.826 | abf | primary_prior | ENST00000613370.1 | yes | — | no annotated structural change |
| IKBIP | SCZ | frontal_cortex_ba9 | chr12:98634413-98644523(-) | 0.853 | abf | primary_prior | ENST00000299157.5,ENST00000342502.6 | no | — | — |
| INO80E | SCZ | cerebellar_hemisphere | chr16:30001040-30001212(+) | 0.958 | susie | intermediate | ENST00000540562.1,ENST00000562441.5,ENST00000567065.5,ENST00000569957.5,ENST00000620599.4 | yes | — | no annotated structural change |
| INO80E | SCZ | cerebellum | chr16:30004657-30005221(+) | 0.953 | susie | intermediate | ENST00000540562.1,ENST00000562441.5,ENST00000567987.5,ENST00000569957.5 | yes | — | no annotated structural change |
| INO80E | SCZ | putamen_basal_ganglia | chr16:30004657-30005221(+) | 0.943 | abf | intermediate | ENST00000540562.1,ENST00000562441.5,ENST00000567987.5,ENST00000569957.5 | no | — | — |
| INO80E | SCZ | nucleus_accumbens_basal_ganglia | chr16:30004657-30005221(+) | 0.936 | abf | intermediate | ENST00000540562.1,ENST00000562441.5,ENST00000567987.5,ENST00000569957.5 | no | — | — |
| INO80E | SCZ | hippocampus | chr16:30004657-30005221(+) | 0.921 | abf | intermediate | ENST00000540562.1,ENST00000562441.5,ENST00000567987.5,ENST00000569957.5 | no | — | — |
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
| KLC1 | SCZ | cortex | chr14:103679518-103684986(+) | 0.855 | abf | primary_prior | ENST00000380038.7,ENST00000445352.8,ENST00000553325.5,ENST00000555856.1 | yes | — | no annotated structural change |
| L3HYPDH | SCZ | cerebellar_hemisphere | chr14:59479351-59505189(-) | 0.944 | susie | intermediate | — | — | — | — |
| L3HYPDH | SCZ | cerebellum | chr14:59504669-59505189(-) | 0.937 | susie | intermediate | — | — | — | — |
| LPCAT4 | SCZ | putamen_basal_ganglia | chr15:34360209-34362196(-) | 0.884 | abf | primary_prior | ENST00000567507.1 | no | — | — |
| LPCAT4 | SCZ | putamen_basal_ganglia | chr15:34360209-34362196(-) | 0.884 | abf | primary_prior | ENST00000567507.1 | no | — | — |
| MAD1L1 | SCZ | anterior_cingulate_cortex_ba24 | chr7:1936897-1957629(-) | 0.992 | abf | robust | ENST00000265854.12,ENST00000399654.6,ENST00000402746.5,ENST00000406869.5,ENST00000450235.5 | no | — | — |
| MAD1L1 | SCZ | frontal_cortex_ba9 | chr7:1936897-1940088(-) | 0.966 | abf | intermediate | ENST00000437877.1 | yes | — | no annotated structural change |
| MAD1L1 | SCZ | cortex | chr7:1936897-1957629(-) | 0.933 | abf | intermediate | ENST00000265854.12,ENST00000399654.6,ENST00000402746.5,ENST00000406869.5,ENST00000450235.5 | no | — | — |
| MAP2K5 | SCZ | cortex | chr15:67563350-67585890(+) | 0.956 | abf | intermediate | — | — | — | — |
| MAP2K5 | SCZ | cerebellum | chr15:67563350-67585890(+) | 0.852 | abf | primary_prior | — | — | — | — |
| MAP2K5 | SCZ | anterior_cingulate_cortex_ba24 | chr15:67563350-67585890(+) | 0.822 | abf | primary_prior | — | — | — | — |
| MAP7D1 | SCZ | nucleus_accumbens_basal_ganglia | chr1:36157290-36170971(+) | 0.961 | abf | intermediate | ENST00000530729.1 | no | — | — |
| MAP7D1 | SCZ | substantia_nigra | chr1:36157290-36170971(+) | 0.911 | abf | intermediate | ENST00000530729.1 | no | — | — |
| MAP7D1 | SCZ | hippocampus | chr1:36157290-36170971(+) | 0.860 | abf | primary_prior | ENST00000530729.1 | no | no | — |
| MED19 | SCZ | hippocampus | chr11:57704396-57704719(-) | 0.810 | abf | primary_prior | ENST00000431606.5 | yes | — | no annotated structural change |
| MRPS33 | SCZ | hypothalamus | chr7:141006535-141014551(-) | 0.994 | abf | robust | ENST00000484502.1 | no | — | — |
| MRPS33 | SCZ | cerebellum | chr7:141006535-141014551(-) | 0.859 | abf | primary_prior | ENST00000484502.1 | no | — | — |
| MRPS33 | SCZ | nucleus_accumbens_basal_ganglia | chr7:141006535-141014551(-) | 0.848 | abf | primary_prior | ENST00000484502.1 | no | — | — |
| NDUFAF7 | SCZ | putamen_basal_ganglia | chr2:37248426-37253086(+) | 0.925 | abf | intermediate | ENST00000441905.1 | yes | — | no annotated structural change |
| NDUFAF7 | SCZ | amygdala | chr2:37248426-37253086(+) | 0.859 | abf | primary_prior | ENST00000441905.1 | no | — | — |
| NEK4 | SCZ | frontal_cortex_ba9 | chr3:52711869-52737586(-) | 0.820 | abf | primary_prior | ENST00000233027.10,ENST00000535191.5 | yes | yes | no annotated structural change |
| NMRAL1 | SCZ | cortex | chr16:4474166-4474554(-) | 0.932 | susie | intermediate | ENST00000283429.11,ENST00000571291.5,ENST00000573520.5,ENST00000575002.5 | yes | — | no annotated structural change |
| NMRAL1 | SCZ | putamen_basal_ganglia | chr16:4474166-4474554(-) | 0.931 | susie | intermediate | ENST00000283429.11,ENST00000571291.5,ENST00000573520.5,ENST00000575002.5 | no | — | — |
| NMRAL1 | SCZ | spinal_cord_cervical_c_1 | chr16:4474166-4474423(-) | 0.925 | susie | intermediate | ENST00000574425.5,ENST00000575995.1 | no | — | — |
| NMRAL1 | SCZ | amygdala | chr16:4474166-4474554(-) | 0.915 | susie | intermediate | ENST00000283429.11,ENST00000571291.5,ENST00000573520.5,ENST00000575002.5 | no | — | — |
| NMRAL1 | SCZ | frontal_cortex_ba9 | chr16:4474166-4474423(-) | 0.914 | susie | intermediate | ENST00000574425.5,ENST00000575995.1 | no | — | — |
| NMRAL1 | SCZ | anterior_cingulate_cortex_ba24 | chr16:4474166-4474423(-) | 0.912 | susie | intermediate | ENST00000574425.5,ENST00000575995.1 | no | — | — |
| NMRAL1 | SCZ | hypothalamus | chr16:4474166-4474423(-) | 0.909 | susie | intermediate | ENST00000574425.5,ENST00000575995.1 | no | — | — |
| NMRAL1 | SCZ | hippocampus | chr16:4474166-4474554(-) | 0.906 | susie | intermediate | ENST00000283429.11,ENST00000571291.5,ENST00000573520.5,ENST00000575002.5 | no | — | — |
| NMRAL1 | SCZ | nucleus_accumbens_basal_ganglia | chr16:4474166-4474423(-) | 0.904 | susie | intermediate | ENST00000574425.5,ENST00000575995.1 | no | — | — |
| NMRAL1 | SCZ | cerebellum | chr16:4474166-4474554(-) | 0.902 | susie | intermediate | ENST00000283429.11,ENST00000571291.5,ENST00000573520.5,ENST00000575002.5 | no | — | — |
| NMRAL1 | SCZ | cerebellar_hemisphere | chr16:4474166-4474554(-) | 0.899 | susie | intermediate | ENST00000283429.11,ENST00000571291.5,ENST00000573520.5,ENST00000575002.5 | no | — | — |
| NMRAL1 | SCZ | substantia_nigra | chr16:4474166-4474554(-) | 0.859 | susie | primary_prior | ENST00000283429.11,ENST00000571291.5,ENST00000573520.5,ENST00000575002.5 | yes | — | no annotated structural change |
| NT5C2 | SCZ | cerebellum | chr10:103174982-103181185(-) | 0.959 | susie | intermediate | ENST00000404739.8,ENST00000467380.1,ENST00000674860.1,ENST00000675020.1,ENST00000675326.1,ENST00000675436.1 | yes | — | no annotated structural change |
| NT5C2 | SCZ | cerebellar_hemisphere | chr10:103174982-103181185(-) | 0.944 | susie | intermediate | ENST00000404739.8,ENST00000467380.1,ENST00000674860.1,ENST00000675020.1,ENST00000675326.1,ENST00000675436.1 | yes | — | no annotated structural change |
| NT5C2 | SCZ | cortex | chr10:103174982-103181185(-) | 0.907 | susie | intermediate | ENST00000404739.8,ENST00000467380.1,ENST00000674860.1,ENST00000675020.1,ENST00000675326.1,ENST00000675436.1 | yes | yes | no annotated structural change |
| NT5C2 | SCZ | putamen_basal_ganglia | chr10:103174982-103181185(-) | 0.907 | abf | intermediate | ENST00000404739.8,ENST00000467380.1,ENST00000674860.1,ENST00000675020.1,ENST00000675326.1,ENST00000675436.1 | yes | — | no annotated structural change |
| NT5C2 | SCZ | caudate_basal_ganglia | chr10:103174982-103181185(-) | 0.901 | susie | intermediate | ENST00000404739.8,ENST00000467380.1,ENST00000674860.1,ENST00000675020.1,ENST00000675326.1,ENST00000675436.1 | no | — | — |
| NT5C2 | SCZ | amygdala | chr10:103174982-103181185(-) | 0.823 | abf | primary_prior | ENST00000404739.8,ENST00000467380.1,ENST00000674860.1,ENST00000675020.1,ENST00000675326.1,ENST00000675436.1 | no | — | — |
| NT5C2 | SCZ | frontal_cortex_ba9 | chr10:103174982-103181185(-) | 0.822 | abf | primary_prior | ENST00000404739.8,ENST00000467380.1,ENST00000674860.1,ENST00000675020.1,ENST00000675326.1,ENST00000675436.1 | no | yes | — |
| NUCB2 | SCZ | hypothalamus | chr11:17312120-17315386(+) | 0.808 | abf | primary_prior | ENST00000323688.10,ENST00000527735.1,ENST00000529010.6,ENST00000531242.1,ENST00000533773.5,ENST00000646648.1 | yes | — | no annotated structural change |
| NUP50 | SCZ | caudate_basal_ganglia | chr22:45183520-45184453(+) | 0.960 | abf | intermediate | ENST00000347635.9,ENST00000396096.6,ENST00000407019.6,ENST00000493456.5 | no | — | — |
| NUP50 | SCZ | nucleus_accumbens_basal_ganglia | chr22:45183520-45184453(+) | 0.956 | abf | intermediate | ENST00000347635.9,ENST00000396096.6,ENST00000407019.6,ENST00000493456.5 | no | — | — |
| NUP50 | SCZ | cerebellum | chr22:45183520-45184453(+) | 0.953 | abf | intermediate | ENST00000347635.9,ENST00000396096.6,ENST00000407019.6,ENST00000493456.5 | no | — | — |
| NUP50 | SCZ | cerebellar_hemisphere | chr22:45183520-45184453(+) | 0.945 | abf | intermediate | ENST00000347635.9,ENST00000396096.6,ENST00000407019.6,ENST00000493456.5 | no | — | — |
| NUP50 | SCZ | cortex | chr22:45164296-45168168(+) | 0.935 | abf | intermediate | ENST00000347635.9,ENST00000407019.6,ENST00000417702.5,ENST00000430547.5,ENST00000434760.5,ENST00000484186.5,ENST00000497960.5 | no | — | — |
| NUP50 | SCZ | frontal_cortex_ba9 | chr22:45164296-45168168(+) | 0.927 | abf | intermediate | ENST00000347635.9,ENST00000407019.6,ENST00000417702.5,ENST00000430547.5,ENST00000434760.5,ENST00000484186.5,ENST00000497960.5 | yes | — | no annotated structural change |
| NUP50 | SCZ | hypothalamus | chr22:45164296-45168168(+) | 0.918 | abf | intermediate | ENST00000347635.9,ENST00000407019.6,ENST00000417702.5,ENST00000430547.5,ENST00000434760.5,ENST00000484186.5,ENST00000497960.5 | no | — | — |
| NUP50 | SCZ | anterior_cingulate_cortex_ba24 | chr22:45164296-45168168(+) | 0.918 | abf | intermediate | ENST00000347635.9,ENST00000407019.6,ENST00000417702.5,ENST00000430547.5,ENST00000434760.5,ENST00000484186.5,ENST00000497960.5 | yes | — | no annotated structural change |
| NUP50 | SCZ | hippocampus | chr22:45164296-45168168(+) | 0.910 | abf | intermediate | ENST00000347635.9,ENST00000407019.6,ENST00000417702.5,ENST00000430547.5,ENST00000434760.5,ENST00000484186.5,ENST00000497960.5 | no | — | — |
| NUP50 | SCZ | amygdala | chr22:45183520-45184453(+) | 0.835 | abf | primary_prior | ENST00000347635.9,ENST00000396096.6,ENST00000407019.6,ENST00000493456.5 | yes | — | no annotated structural change |
| NUP50 | SCZ | putamen_basal_ganglia | chr22:45164296-45168168(+) | 0.805 | abf | primary_prior | ENST00000347635.9,ENST00000407019.6,ENST00000417702.5,ENST00000430547.5,ENST00000434760.5,ENST00000484186.5,ENST00000497960.5 | no | — | — |
| PAK6 | SCZ | cerebellum | chr15:40266495-40272224(+) | 0.866 | abf | primary_prior | ENST00000260404.8,ENST00000441369.6,ENST00000453867.7,ENST00000455577.6,ENST00000542403.3,ENST00000558658.6,ENST00000560346.6 | yes | — | no annotated structural change |
| PAK6 | SCZ | cerebellum | chr15:40266495-40272224(+) | 0.866 | abf | primary_prior | ENST00000260404.8,ENST00000441369.6,ENST00000453867.7,ENST00000455577.6,ENST00000542403.3,ENST00000558658.6,ENST00000560346.6 | yes | — | no annotated structural change |
| PAM16 | SCZ | amygdala | chr16:4355163-4355286(-) | 0.908 | susie | intermediate | — | — | — | — |
| PBRM1 | SCZ | cerebellum | chr3:52644789-52648344(-) | 0.985 | abf | robust | ENST00000296302.11,ENST00000337303.8,ENST00000356770.8,ENST00000394830.7,ENST00000409057.5,ENST00000409114.7,ENST00000409767.5,ENST00000410007.5,ENST00000412587.5,ENST00000423351.5,ENST00000446103.5,ENST00000707071.1 | no | — | — |
| PCBP3 | SCZ | cerebellar_hemisphere | chr21:45735429-45849961(+) | 0.819 | abf | primary_prior | ENST00000400314.5 | yes | — | no annotated structural change |
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
| PPIP5K1 | SCZ | cerebellar_hemisphere | chr15:43558932-43560413(-) | 0.926 | abf | intermediate | ENST00000420765.6,ENST00000439195.5,ENST00000644537.1 | yes | — | no annotated structural change |
| PPIP5K1 | SCZ | cerebellum | chr15:43558932-43564103(-) | 0.869 | abf | primary_prior | ENST00000381879.8,ENST00000381885.5,ENST00000396923.7 | yes | — | no annotated structural change |
| PRMT7 | SCZ | amygdala | chr16:68356215-68356701(+) | 0.934 | abf | intermediate | ENST00000675132.1,ENST00000685141.1,ENST00000692283.1 | no | — | — |
| PSMD6 | SCZ | hippocampus | chr3:64010955-64018599(-) | 0.976 | abf | robust | ENST00000480205.5 | no | no | — |
| RAI1 | SCZ | cerebellum | chr17:17793288-17803756(+) | 0.971 | abf | intermediate | — | — | — | — |
| RAI1 | SCZ | cerebellar_hemisphere | chr17:17798513-17803756(+) | 0.946 | susie | intermediate | ENST00000353383.6,ENST00000583166.1 | no | — | — |
| RBM6 | SCZ | nucleus_accumbens_basal_ganglia | chr3:49962685-49967470(+) | 0.869 | abf | primary_prior | ENST00000266022.9,ENST00000425608.5,ENST00000433811.1 | no | — | — |
| RCBTB1 | SCZ | hippocampus | chr13:49534262-49540876(-) | 0.886 | abf | primary_prior | ENST00000258646.3,ENST00000378302.7 | no | yes | — |
| REEP2 | SCZ | cortex | chr5:138441461-138443587(+) | 0.913 | abf | intermediate | ENST00000464751.6 | no | — | — |
| RERE | SCZ | nucleus_accumbens_basal_ganglia | chr1:8438300-8465925(-) | 0.856 | abf | primary_prior | ENST00000659924.1 | no | — | — |
| RERE | SCZ | cerebellar_hemisphere | chr1:8438300-8465925(-) | 0.818 | abf | primary_prior | ENST00000659924.1 | no | — | — |
| SETD6 | SCZ | anterior_cingulate_cortex_ba24 | chr16:58516343-58518401(+) | 0.926 | abf | intermediate | ENST00000447443.1 | no | — | — |
| SNAP91 | SCZ | cerebellum | chr6:83607808-83610650(-) | 0.974 | susie | intermediate | ENST00000195649.10,ENST00000369694.7,ENST00000439399.6,ENST00000518312.5,ENST00000520213.5,ENST00000520302.5,ENST00000521485.5,ENST00000521616.5,ENST00000521743.5,ENST00000521931.5 | yes | — | no annotated structural change |
| SNAP91 | SCZ | nucleus_accumbens_basal_ganglia | chr6:83607808-83610650(-) | 0.964 | abf | intermediate | ENST00000195649.10,ENST00000369694.7,ENST00000439399.6,ENST00000518312.5,ENST00000520213.5,ENST00000520302.5,ENST00000521485.5,ENST00000521616.5,ENST00000521743.5,ENST00000521931.5 | no | — | — |
| SNAP91 | SCZ | spinal_cord_cervical_c_1 | chr6:83607808-83610650(-) | 0.962 | abf | intermediate | ENST00000195649.10,ENST00000369694.7,ENST00000439399.6,ENST00000518312.5,ENST00000520213.5,ENST00000520302.5,ENST00000521485.5,ENST00000521616.5,ENST00000521743.5,ENST00000521931.5 | no | — | — |
| SNAP91 | SCZ | hypothalamus | chr6:83607808-83610650(-) | 0.961 | abf | intermediate | ENST00000195649.10,ENST00000369694.7,ENST00000439399.6,ENST00000518312.5,ENST00000520213.5,ENST00000520302.5,ENST00000521485.5,ENST00000521616.5,ENST00000521743.5,ENST00000521931.5 | yes | — | no annotated structural change |
| SNAP91 | SCZ | cortex | chr6:83607808-83610650(-) | 0.958 | abf | intermediate | ENST00000195649.10,ENST00000369694.7,ENST00000439399.6,ENST00000518312.5,ENST00000520213.5,ENST00000520302.5,ENST00000521485.5,ENST00000521616.5,ENST00000521743.5,ENST00000521931.5 | yes | — | no annotated structural change |
| SNAP91 | SCZ | cerebellar_hemisphere | chr6:83593017-83593478(-) | 0.955 | susie | intermediate | ENST00000518312.5 | no | — | — |
| SNAP91 | SCZ | caudate_basal_ganglia | chr6:83607808-83610650(-) | 0.954 | abf | intermediate | ENST00000195649.10,ENST00000369694.7,ENST00000439399.6,ENST00000518312.5,ENST00000520213.5,ENST00000520302.5,ENST00000521485.5,ENST00000521616.5,ENST00000521743.5,ENST00000521931.5 | no | yes | — |
| SYT5 | SCZ | frontal_cortex_ba9 | chr19:55179086-55179964(-) | 0.904 | abf | intermediate | ENST00000589172.5 | no | — | — |
| TAOK2 | SCZ | cerebellum | chr16:29986504-29986844(+) | 0.829 | abf | primary_prior | ENST00000543033.5 | yes | — | no annotated structural change |
| TEAD4 | SCZ | putamen_basal_ganglia | chr12:2959481-2959948(+) | 0.908 | abf | intermediate | ENST00000358409.7,ENST00000359864.8,ENST00000536826.2,ENST00000540314.2 | no | — | — |
| THAP3 | SCZ | cerebellar_hemisphere | chr1:6628691-6629643(+) | 0.964 | susie | intermediate | ENST00000487819.5 | yes | — | no annotated structural change |
| THAP3 | SCZ | cerebellum | chr1:6628691-6629643(+) | 0.860 | abf | primary_prior | ENST00000487819.5 | no | — | — |
| TMED4 | SCZ | hypothalamus | chr7:44579628-44581449(-) | 0.956 | abf | intermediate | ENST00000289577.10 | no | — | — |
| TMED4 | SCZ | cerebellum | chr7:44579628-44581093(-) | 0.943 | abf | intermediate | ENST00000457408.7 | yes | — | no annotated structural change |
| TMED4 | SCZ | nucleus_accumbens_basal_ganglia | chr7:44579628-44581449(-) | 0.941 | abf | intermediate | ENST00000289577.10 | no | — | — |
| TMED4 | SCZ | putamen_basal_ganglia | chr7:44579628-44581449(-) | 0.940 | abf | intermediate | ENST00000289577.10 | no | — | — |
| TMED4 | SCZ | caudate_basal_ganglia | chr7:44579628-44581449(-) | 0.939 | abf | intermediate | ENST00000289577.10 | no | yes | — |
| TMED4 | SCZ | hippocampus | chr7:44579628-44581093(-) | 0.939 | abf | intermediate | ENST00000457408.7 | no | — | — |
| TMED4 | SCZ | cortex | chr7:44579628-44581449(-) | 0.938 | abf | intermediate | ENST00000289577.10 | yes | yes | no annotated structural change |
| TMED4 | SCZ | cerebellar_hemisphere | chr7:44579628-44581449(-) | 0.937 | abf | intermediate | ENST00000289577.10 | no | — | — |
| TMED4 | SCZ | frontal_cortex_ba9 | chr7:44579628-44581449(-) | 0.937 | abf | intermediate | ENST00000289577.10 | yes | yes | no annotated structural change |
| TMED4 | SCZ | amygdala | chr7:44579628-44581449(-) | 0.936 | abf | intermediate | ENST00000289577.10 | no | — | — |
| TMED4 | SCZ | spinal_cord_cervical_c_1 | chr7:44579628-44581093(-) | 0.922 | abf | intermediate | ENST00000457408.7 | yes | — | no annotated structural change |
| TMED4 | SCZ | substantia_nigra | chr7:44579628-44581449(-) | 0.914 | abf | intermediate | ENST00000289577.10 | no | — | — |
| TMED4 | SCZ | anterior_cingulate_cortex_ba24 | chr7:44579628-44581449(-) | 0.850 | abf | primary_prior | ENST00000289577.10 | yes | — | no annotated structural change |
| TSPAN31 | SCZ | frontal_cortex_ba9 | chr12:57738201-57745745(+) | 0.813 | abf | primary_prior | — | — | — | — |
| TUBGCP4 | SCZ | cerebellum | chr15:43398179-43400044(+) | 0.882 | abf | primary_prior | ENST00000260383.11,ENST00000563147.5,ENST00000563963.1,ENST00000564079.6 | no | — | — |
| TUBGCP4 | SCZ | cerebellum | chr15:43398179-43400044(+) | 0.882 | abf | primary_prior | ENST00000260383.11,ENST00000563147.5,ENST00000563963.1,ENST00000564079.6 | no | — | — |
| YPEL1 | SCZ | cerebellar_hemisphere | chr22:21701218-21703370(-) | 0.985 | abf | robust | ENST00000339468.8,ENST00000477675.1,ENST00000672036.2 | no | — | — |
| YPEL1 | SCZ | cerebellum | chr22:21701218-21703370(-) | 0.985 | abf | robust | ENST00000339468.8,ENST00000477675.1,ENST00000672036.2 | no | — | — |
| YPEL1 | SCZ | caudate_basal_ganglia | chr22:21701218-21703370(-) | 0.984 | abf | robust | ENST00000339468.8,ENST00000477675.1,ENST00000672036.2 | no | — | — |
| YPEL1 | SCZ | putamen_basal_ganglia | chr22:21701218-21703370(-) | 0.982 | abf | robust | ENST00000339468.8,ENST00000477675.1,ENST00000672036.2 | no | — | — |
| YPEL1 | SCZ | frontal_cortex_ba9 | chr22:21701218-21703370(-) | 0.982 | abf | robust | ENST00000339468.8,ENST00000477675.1,ENST00000672036.2 | yes | — | no annotated structural change |
| YPEL1 | SCZ | nucleus_accumbens_basal_ganglia | chr22:21701218-21703370(-) | 0.904 | abf | intermediate | ENST00000339468.8,ENST00000477675.1,ENST00000672036.2 | no | — | — |
| YPEL1 | SCZ | anterior_cingulate_cortex_ba24 | chr22:21701218-21703370(-) | 0.899 | abf | intermediate | ENST00000339468.8,ENST00000477675.1,ENST00000672036.2 | yes | — | no annotated structural change |
| YPEL1 | SCZ | cortex | chr22:21701218-21703370(-) | 0.866 | abf | primary_prior | ENST00000339468.8,ENST00000477675.1,ENST00000672036.2 | no | — | — |
| YPEL3 | SCZ | anterior_cingulate_cortex_ba24 | chr16:30094897-30095252(-) | 0.845 | abf | primary_prior | — | — | — | — |
| YWHAB | SCZ | frontal_cortex_ba9 | chr20:44885886-44901531(+) | 0.891 | abf | intermediate | ENST00000353703.9,ENST00000479421.5 | yes | — | no annotated structural change |
| YWHAB | SCZ | nucleus_accumbens_basal_ganglia | chr20:44885886-44901531(+) | 0.881 | abf | primary_prior | ENST00000353703.9,ENST00000479421.5 | no | — | — |
| YWHAB | SCZ | spinal_cord_cervical_c_1 | chr20:44885886-44901531(+) | 0.861 | abf | primary_prior | ENST00000353703.9,ENST00000479421.5 | no | — | — |
| YWHAB | SCZ | cerebellar_hemisphere | chr20:44885886-44901531(+) | 0.860 | abf | primary_prior | ENST00000353703.9,ENST00000479421.5 | no | — | — |
| YWHAB | SCZ | cortex | chr20:44885886-44901531(+) | 0.856 | abf | primary_prior | ENST00000353703.9,ENST00000479421.5 | yes | — | no annotated structural change |
| YWHAB | SCZ | cerebellum | chr20:44885886-44901531(+) | 0.834 | abf | primary_prior | ENST00000353703.9,ENST00000479421.5 | no | — | — |
| YWHAB | SCZ | caudate_basal_ganglia | chr20:44885886-44901531(+) | 0.832 | abf | primary_prior | ENST00000353703.9,ENST00000479421.5 | no | — | — |
| YWHAB | SCZ | hippocampus | chr20:44885886-44901531(+) | 0.803 | abf | primary_prior | ENST00000353703.9,ENST00000479421.5 | yes | yes | no annotated structural change |
| YWHAB | SCZ | hypothalamus | chr20:44885886-44901531(+) | 0.802 | abf | primary_prior | ENST00000353703.9,ENST00000479421.5 | no | — | — |
| ZDHHC12 | SCZ | spinal_cord_cervical_c_1 | chr9:128721502-128721651(-) | 0.998 | abf | robust | ENST00000372663.9,ENST00000372667.9 | no | — | — |
| ZDHHC12 | SCZ | substantia_nigra | chr9:128721502-128721651(-) | 0.957 | abf | intermediate | ENST00000372663.9,ENST00000372667.9 | no | — | — |
| ZDHHC12 | SCZ | anterior_cingulate_cortex_ba24 | chr9:128721480-128721651(-) | 0.932 | abf | intermediate | ENST00000467312.1 | no | — | — |
| ZDHHC12 | SCZ | cerebellar_hemisphere | chr9:128721502-128721651(-) | 0.931 | abf | intermediate | ENST00000372663.9,ENST00000372667.9 | yes | — | no annotated structural change |
| ZDHHC12 | SCZ | caudate_basal_ganglia | chr9:128721502-128721651(-) | 0.931 | abf | intermediate | ENST00000372663.9,ENST00000372667.9 | no | — | — |
| ZDHHC12 | SCZ | nucleus_accumbens_basal_ganglia | chr9:128721480-128721651(-) | 0.931 | abf | intermediate | ENST00000467312.1 | no | — | — |
| ZDHHC12 | SCZ | hippocampus | chr9:128721502-128721651(-) | 0.930 | abf | intermediate | ENST00000372663.9,ENST00000372667.9 | no | — | — |
| ZDHHC12 | SCZ | amygdala | chr9:128721502-128721651(-) | 0.928 | abf | intermediate | ENST00000372663.9,ENST00000372667.9 | no | — | — |
| ZDHHC12 | SCZ | cortex | chr9:128721480-128721651(-) | 0.927 | abf | intermediate | ENST00000467312.1 | no | — | — |
| ZDHHC12 | SCZ | putamen_basal_ganglia | chr9:128721480-128721651(-) | 0.923 | abf | intermediate | ENST00000467312.1 | no | — | — |
| ZDHHC12 | SCZ | frontal_cortex_ba9 | chr9:128721480-128721651(-) | 0.923 | abf | intermediate | ENST00000467312.1 | no | — | — |
| ZDHHC12 | SCZ | cerebellum | chr9:128721502-128721651(-) | 0.915 | abf | intermediate | ENST00000372663.9,ENST00000372667.9 | yes | — | no annotated structural change |
| ZDHHC12 | SCZ | hypothalamus | chr9:128721502-128721651(-) | 0.884 | abf | primary_prior | ENST00000372663.9,ENST00000372667.9 | no | — | — |
| ZFYVE21 | SCZ | caudate_basal_ganglia | chr14:103729182-103732620(+) | 0.959 | abf | intermediate | ENST00000311141.7 | no | yes | — |
| ZFYVE21 | SCZ | caudate_basal_ganglia | chr14:103729182-103732620(+) | 0.955 | susie | intermediate | ENST00000311141.7 | no | yes | — |
| ZFYVE21 | SCZ | nucleus_accumbens_basal_ganglia | chr14:103729182-103732620(+) | 0.955 | susie | intermediate | ENST00000311141.7 | no | — | — |
| ZFYVE21 | SCZ | putamen_basal_ganglia | chr14:103729182-103729793(+) | 0.924 | susie | intermediate | ENST00000216602.10,ENST00000555163.2,ENST00000555501.1,ENST00000556610.5 | no | — | — |
