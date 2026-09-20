# Short-read junction confirmation of the anchored switch pairs

**30** genes carry a junction in the tissue-matched IsoGraph switch pair. **21** sit only in GTEx tissues BrainSEQ does not sequence and cannot be tested here at all. Of the **9** that reach a BrainSEQ region, **6** have the junction measured by a PSI event and **5** validate on at least one event.

alpha = 0.05, minor-form usage threshold = 0.05. Rows below are one per candidate PSI event (173 rows), not one per gene.

Every gene whose junction lands in the tissue-matched IsoGraph switch pair is
tested, not only the ones a figure panel highlights.

| verdict | rows |
| --- | ---: |
| validated | 63 |
| no_matched_brainseq_region | 60 |
| not_validated | 44 |
| junction_not_measured | 6 |

## Per gene, in the regions BrainSEQ can measure

`events` counts candidate PSI events, `validated` how many of them pass the usage threshold. A gene with a mixed count validates on some contrasts and not others, which is why the display-item rule below names one.

| gene | trait | region | match | events | validated | call |
| --- | --- | --- | --- | ---: | ---: | --- |
| B3GAT1 | scz | caudate | secondary | 6 | 6 | validates |
| B3GAT1 | scz | dlpfc | exact | 6 | 6 | validates |
| GSTO2 | scz | caudate | secondary | 6 | 6 | validates |
| GSTO2 | scz | dlpfc | exact | 6 | 6 | validates |
| PIGQ | als | caudate | adjacent | 21 | 8 | mixed |
| PPP6R2 | scz | caudate | secondary | 0 | 0 | junction not measured |
| PPP6R2 | scz | dlpfc | exact | 0 | 0 | junction not measured |
| PRDM2 | als | caudate | secondary | 0 | 0 | junction not measured |
| PRDM2 | als | dlpfc | exact | 0 | 0 | junction not measured |
| RPAIN | scz | caudate | secondary | 23 | 13 | mixed |
| RPAIN | scz | dlpfc | exact | 23 | 13 | mixed |
| TBC1D15 | pd | caudate | secondary | 0 | 0 | junction not measured |
| TBC1D15 | pd | dlpfc | exact | 0 | 0 | junction not measured |
| TMEM175 | pd | caudate | secondary | 6 | 3 | mixed |
| TMEM175 | pd | dlpfc | exact | 6 | 2 | mixed |
| VAMP2 | als | caudate | secondary | 2 | 0 | does not validate |
| VAMP2 | als | dlpfc | exact | 2 | 0 | does not validate |

`direct` rows test a single PSI event whose two arms contrast the anchored
junction against the alternative form, so PSI is the switch ratio itself and there
is no compositional-closure confound; the statistic is minor-form usage,
`min(median PSI, 1 - median PSI)`, which does not depend on PSI orientation.
`paired` rows test PSI anti-correlation against the closure baseline.
`no_matched_brainseq_region` rows are concordant events in a GTEx tissue BrainSEQ
does not sequence (chiefly cerebellum); they are untestable here, not negative.
`match` says how close the BrainSEQ region is to the GTEx tissue: exact, adjacent
(same division, different structure) or secondary (caudate carried alongside DLPFC).

| gene | trait | tissue | region | match | mode | n | minor-form usage | median PSI | rho | q_within | q_matched | verdict |
| --- | --- | --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| PPP6R2 | scz | Frontal_Cortex_BA9 | caudate | secondary | none | -- | -- | -- | -- | -- | -- | junction_not_measured |
| PPP6R2 | scz | Frontal_Cortex_BA9 | dlpfc | exact | none | -- | -- | -- | -- | -- | -- | junction_not_measured |
| PRDM2 | als | Cortex | caudate | secondary | none | -- | -- | -- | -- | -- | -- | junction_not_measured |
| PRDM2 | als | Cortex | dlpfc | exact | none | -- | -- | -- | -- | -- | -- | junction_not_measured |
| TBC1D15 | pd | Frontal_Cortex_BA9 | caudate | secondary | none | -- | -- | -- | -- | -- | -- | junction_not_measured |
| TBC1D15 | pd | Frontal_Cortex_BA9 | dlpfc | exact | none | -- | -- | -- | -- | -- | -- | junction_not_measured |
| BAIAP3 | als | Hypothalamus | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| CRELD2 | scz | Cerebellum | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| CRELD2 | scz | Cerebellum | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| CRELD2 | scz | Cerebellum | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| CRELD2 | scz | Cerebellum | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| CRELD2 | scz | Cerebellum | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| CTC1 | als | Cerebellum | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| CTC1 | als | Cerebellum | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| CTSB | pd | Amygdala | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| DLG1 | scz | Cerebellar_Hemisphere | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| DLG1 | scz | Cerebellar_Hemisphere | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| DNAJA3 | scz | Cerebellar_Hemisphere | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| DNAJA3 | scz | Cerebellar_Hemisphere | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| DOC2A | scz | Cerebellar_Hemisphere | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| FLCN | ad | Cerebellum | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| FLCN | ad | Cerebellum | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| FLCN | ad | Cerebellum | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| FLCN | scz | Cerebellar_Hemisphere | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| FLCN | scz | Cerebellar_Hemisphere | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| FLCN | scz | Cerebellar_Hemisphere | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| GPR135 | scz | Cerebellum | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| IFNAR2 | ad | Cerebellum | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| IFNAR2 | ad | Cerebellum | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| INO80E | scz | Cerebellar_Hemisphere | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| INO80E | scz | Cerebellar_Hemisphere | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| INO80E | scz | Cerebellar_Hemisphere | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| NADSYN1 | scz | Cerebellar_Hemisphere | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| NADSYN1 | scz | Cerebellar_Hemisphere | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| NADSYN1 | scz | Cerebellar_Hemisphere | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| NADSYN1 | scz | Cerebellar_Hemisphere | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| NADSYN1 | scz | Cerebellar_Hemisphere | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| NADSYN1 | scz | Cerebellar_Hemisphere | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| NADSYN1 | scz | Cerebellar_Hemisphere | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| NADSYN1 | scz | Cerebellar_Hemisphere | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| NT5C2 | scz | Cerebellum | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| PCGF3 | pd | Cerebellum | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| PCGF3 | pd | Cerebellum | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| PTPRN | als | Cerebellum | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| RBFA | scz | Cerebellum | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| SNAP91 | scz | Cerebellum | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| SPG7 | scz | Cerebellum | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| SPG7 | scz | Cerebellum | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| SPG7 | scz | Cerebellum | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| SPG7 | scz | Cerebellum | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| SPG7 | scz | Cerebellum | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| SPG7 | scz | Cerebellum | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| SPG7 | scz | Cerebellum | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| SPG7 | scz | Cerebellum | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| SPG7 | scz | Cerebellum | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| SPG7 | scz | Cerebellum | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| TARBP1 | scz | Cerebellum | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| TARBP1 | scz | Cerebellum | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| TARBP1 | scz | Cerebellum | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| THAP3 | scz | Cerebellar_Hemisphere | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| TMEM175 | ad | Cerebellum | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| TMEM175 | ad | Cerebellum | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| TMEM175 | als | Cerebellar_Hemisphere | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| TMEM175 | als | Cerebellar_Hemisphere | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| TMEM175 | als | Cerebellar_Hemisphere | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| TPP1 | als | Cerebellum | -- | none | none | -- | -- | -- | -- | -- | -- | no_matched_brainseq_region |
| PIGQ | als | Putamen_basal_ganglia | caudate | adjacent | direct | 238 | 0.0000 | 1.000 | -- | -- | -- | not_validated |
| PIGQ | als | Putamen_basal_ganglia | caudate | adjacent | direct | 235 | 0.0000 | 1.000 | -- | -- | -- | not_validated |
| PIGQ | als | Putamen_basal_ganglia | caudate | adjacent | direct | 214 | 0.0000 | 1.000 | -- | -- | -- | not_validated |
| PIGQ | als | Putamen_basal_ganglia | caudate | adjacent | direct | 117 | 0.0000 | 0.000 | -- | -- | -- | not_validated |
| PIGQ | als | Putamen_basal_ganglia | caudate | adjacent | direct | 214 | 0.0000 | 1.000 | -- | -- | -- | not_validated |
| PIGQ | als | Putamen_basal_ganglia | caudate | adjacent | direct | 238 | 0.0000 | 1.000 | -- | -- | -- | not_validated |
| PIGQ | als | Putamen_basal_ganglia | caudate | adjacent | direct | 214 | 0.0000 | 1.000 | -- | -- | -- | not_validated |
| PIGQ | als | Putamen_basal_ganglia | caudate | adjacent | direct | 212 | 0.0000 | 1.000 | -- | -- | -- | not_validated |
| PIGQ | als | Putamen_basal_ganglia | caudate | adjacent | direct | 204 | 0.0000 | 1.000 | -- | -- | -- | not_validated |
| PIGQ | als | Putamen_basal_ganglia | caudate | adjacent | direct | 201 | 0.0000 | 1.000 | -- | -- | -- | not_validated |
| PIGQ | als | Putamen_basal_ganglia | caudate | adjacent | direct | 198 | 0.0000 | 1.000 | -- | -- | -- | not_validated |
| PIGQ | als | Putamen_basal_ganglia | caudate | adjacent | direct | 195 | 0.0000 | 1.000 | -- | -- | -- | not_validated |
| PIGQ | als | Putamen_basal_ganglia | caudate | adjacent | direct | 195 | 0.0000 | 1.000 | -- | -- | -- | not_validated |
| RPAIN | scz | Cortex | caudate | secondary | direct | 238 | 0.0312 | 0.031 | -- | -- | -- | not_validated |
| RPAIN | scz | Cortex | caudate | secondary | direct | 238 | 0.0000 | 1.000 | -- | -- | -- | not_validated |
| RPAIN | scz | Cortex | caudate | secondary | direct | 238 | 0.0434 | 0.957 | -- | -- | -- | not_validated |
| RPAIN | scz | Cortex | caudate | secondary | direct | 238 | 0.0081 | 0.992 | -- | -- | -- | not_validated |
| RPAIN | scz | Cortex | caudate | secondary | direct | 236 | 0.0000 | 0.000 | -- | -- | -- | not_validated |
| RPAIN | scz | Cortex | caudate | secondary | direct | 236 | 0.0000 | 0.000 | -- | -- | -- | not_validated |
| RPAIN | scz | Cortex | caudate | secondary | direct | 228 | 0.0000 | 1.000 | -- | -- | -- | not_validated |
| RPAIN | scz | Cortex | caudate | secondary | direct | 238 | 0.0000 | 0.000 | -- | -- | -- | not_validated |
| RPAIN | scz | Cortex | caudate | secondary | direct | 229 | 0.0000 | 1.000 | -- | -- | -- | not_validated |
| RPAIN | scz | Cortex | caudate | secondary | direct | 228 | 0.0000 | 1.000 | -- | -- | -- | not_validated |
| RPAIN | scz | Cortex | dlpfc | exact | direct | 222 | 0.0179 | 0.018 | -- | -- | -- | not_validated |
| RPAIN | scz | Cortex | dlpfc | exact | direct | 222 | 0.0000 | 1.000 | -- | -- | -- | not_validated |
| RPAIN | scz | Cortex | dlpfc | exact | direct | 222 | 0.0431 | 0.957 | -- | -- | -- | not_validated |
| RPAIN | scz | Cortex | dlpfc | exact | direct | 222 | 0.0058 | 0.994 | -- | -- | -- | not_validated |
| RPAIN | scz | Cortex | dlpfc | exact | direct | 219 | 0.0000 | 0.000 | -- | -- | -- | not_validated |
| RPAIN | scz | Cortex | dlpfc | exact | direct | 218 | 0.0000 | 0.000 | -- | -- | -- | not_validated |
| RPAIN | scz | Cortex | dlpfc | exact | direct | 205 | 0.0000 | 1.000 | -- | -- | -- | not_validated |
| RPAIN | scz | Cortex | dlpfc | exact | direct | 222 | 0.0000 | 0.000 | -- | -- | -- | not_validated |
| RPAIN | scz | Cortex | dlpfc | exact | direct | 212 | 0.0000 | 1.000 | -- | -- | -- | not_validated |
| RPAIN | scz | Cortex | dlpfc | exact | direct | 205 | 0.0000 | 1.000 | -- | -- | -- | not_validated |
| TMEM175 | pd | Cortex | caudate | secondary | direct | 238 | 0.0270 | 0.027 | -- | -- | -- | not_validated |
| TMEM175 | pd | Cortex | caudate | secondary | direct | 238 | 0.0000 | 0.000 | -- | -- | -- | not_validated |
| TMEM175 | pd | Cortex | caudate | secondary | direct | 156 | 0.0000 | 0.000 | -- | -- | -- | not_validated |
| TMEM175 | pd | Cortex | dlpfc | exact | direct | 222 | 0.0269 | 0.027 | -- | -- | -- | not_validated |
| TMEM175 | pd | Cortex | dlpfc | exact | direct | 222 | 0.0265 | 0.974 | -- | -- | -- | not_validated |
| TMEM175 | pd | Cortex | dlpfc | exact | direct | 220 | 0.0000 | 0.000 | -- | -- | -- | not_validated |
| TMEM175 | pd | Cortex | dlpfc | exact | direct | 129 | 0.0000 | 0.000 | -- | -- | -- | not_validated |
| VAMP2 | als | Cortex | caudate | secondary | direct | 238 | 0.0137 | 0.014 | -- | -- | -- | not_validated |
| VAMP2 | als | Cortex | caudate | secondary | direct | 238 | 0.0180 | 0.982 | -- | -- | -- | not_validated |
| VAMP2 | als | Cortex | dlpfc | exact | direct | 222 | 0.0114 | 0.011 | -- | -- | -- | not_validated |
| VAMP2 | als | Cortex | dlpfc | exact | direct | 222 | 0.0165 | 0.984 | -- | -- | -- | not_validated |
| B3GAT1 | scz | Cortex | caudate | secondary | direct | 238 | 0.4521 | 0.452 | -- | -- | -- | validated |
| B3GAT1 | scz | Cortex | caudate | secondary | direct | 238 | 0.0872 | 0.913 | -- | -- | -- | validated |
| B3GAT1 | scz | Cortex | caudate | secondary | direct | 238 | 0.3362 | 0.336 | -- | -- | -- | validated |
| B3GAT1 | scz | Cortex | caudate | secondary | direct | 238 | 0.4521 | 0.452 | -- | -- | -- | validated |
| B3GAT1 | scz | Cortex | caudate | secondary | direct | 238 | 0.0872 | 0.913 | -- | -- | -- | validated |
| B3GAT1 | scz | Cortex | caudate | secondary | direct | 238 | 0.3362 | 0.336 | -- | -- | -- | validated |
| B3GAT1 | scz | Cortex | dlpfc | exact | direct | 222 | 0.4313 | 0.431 | -- | -- | -- | validated |
| B3GAT1 | scz | Cortex | dlpfc | exact | direct | 222 | 0.1042 | 0.896 | -- | -- | -- | validated |
| B3GAT1 | scz | Cortex | dlpfc | exact | direct | 222 | 0.2089 | 0.209 | -- | -- | -- | validated |
| B3GAT1 | scz | Cortex | dlpfc | exact | direct | 222 | 0.4313 | 0.431 | -- | -- | -- | validated |
| B3GAT1 | scz | Cortex | dlpfc | exact | direct | 222 | 0.1042 | 0.896 | -- | -- | -- | validated |
| B3GAT1 | scz | Cortex | dlpfc | exact | direct | 222 | 0.2089 | 0.209 | -- | -- | -- | validated |
| GSTO2 | scz | Frontal_Cortex_BA9 | caudate | secondary | direct | 238 | 0.4022 | 0.598 | -- | -- | -- | validated |
| GSTO2 | scz | Frontal_Cortex_BA9 | caudate | secondary | direct | 238 | 0.2166 | 0.217 | -- | -- | -- | validated |
| GSTO2 | scz | Frontal_Cortex_BA9 | caudate | secondary | direct | 238 | 0.2832 | 0.283 | -- | -- | -- | validated |
| GSTO2 | scz | Frontal_Cortex_BA9 | caudate | secondary | direct | 238 | 0.3699 | 0.370 | -- | -- | -- | validated |
| GSTO2 | scz | Frontal_Cortex_BA9 | caudate | secondary | direct | 238 | 0.3699 | 0.370 | -- | -- | -- | validated |
| GSTO2 | scz | Frontal_Cortex_BA9 | caudate | secondary | direct | 238 | 0.3699 | 0.370 | -- | -- | -- | validated |
| GSTO2 | scz | Frontal_Cortex_BA9 | dlpfc | exact | direct | 222 | 0.4504 | 0.550 | -- | -- | -- | validated |
| GSTO2 | scz | Frontal_Cortex_BA9 | dlpfc | exact | direct | 222 | 0.1797 | 0.180 | -- | -- | -- | validated |
| GSTO2 | scz | Frontal_Cortex_BA9 | dlpfc | exact | direct | 222 | 0.2696 | 0.270 | -- | -- | -- | validated |
| GSTO2 | scz | Frontal_Cortex_BA9 | dlpfc | exact | direct | 222 | 0.4318 | 0.432 | -- | -- | -- | validated |
| GSTO2 | scz | Frontal_Cortex_BA9 | dlpfc | exact | direct | 222 | 0.4318 | 0.432 | -- | -- | -- | validated |
| GSTO2 | scz | Frontal_Cortex_BA9 | dlpfc | exact | direct | 222 | 0.4318 | 0.432 | -- | -- | -- | validated |
| PIGQ | als | Putamen_basal_ganglia | caudate | adjacent | direct | 238 | 0.1800 | 0.180 | -- | -- | -- | validated |
| PIGQ | als | Putamen_basal_ganglia | caudate | adjacent | direct | 238 | 0.1877 | 0.812 | -- | -- | -- | validated |
| PIGQ | als | Putamen_basal_ganglia | caudate | adjacent | direct | 238 | 0.1877 | 0.812 | -- | -- | -- | validated |
| PIGQ | als | Putamen_basal_ganglia | caudate | adjacent | direct | 238 | 0.4684 | 0.468 | -- | -- | -- | validated |
| PIGQ | als | Putamen_basal_ganglia | caudate | adjacent | direct | 238 | 0.1601 | 0.160 | -- | -- | -- | validated |
| PIGQ | als | Putamen_basal_ganglia | caudate | adjacent | direct | 238 | 0.1877 | 0.812 | -- | -- | -- | validated |
| PIGQ | als | Putamen_basal_ganglia | caudate | adjacent | direct | 232 | 0.2087 | 0.209 | -- | -- | -- | validated |
| PIGQ | als | Putamen_basal_ganglia | caudate | adjacent | direct | 222 | 0.2923 | 0.292 | -- | -- | -- | validated |
| RPAIN | scz | Cortex | caudate | secondary | direct | 238 | 0.4727 | 0.527 | -- | -- | -- | validated |
| RPAIN | scz | Cortex | caudate | secondary | direct | 238 | 0.1755 | 0.825 | -- | -- | -- | validated |
| RPAIN | scz | Cortex | caudate | secondary | direct | 219 | 0.1932 | 0.807 | -- | -- | -- | validated |
| RPAIN | scz | Cortex | caudate | secondary | direct | 190 | 0.3469 | 0.347 | -- | -- | -- | validated |
| RPAIN | scz | Cortex | caudate | secondary | direct | 189 | 0.3533 | 0.353 | -- | -- | -- | validated |
| RPAIN | scz | Cortex | caudate | secondary | direct | 164 | 0.4273 | 0.427 | -- | -- | -- | validated |
| RPAIN | scz | Cortex | caudate | secondary | direct | 238 | 0.3589 | 0.641 | -- | -- | -- | validated |
| RPAIN | scz | Cortex | caudate | secondary | direct | 238 | 0.4787 | 0.479 | -- | -- | -- | validated |
| RPAIN | scz | Cortex | caudate | secondary | direct | 238 | 0.1942 | 0.194 | -- | -- | -- | validated |
| RPAIN | scz | Cortex | caudate | secondary | direct | 238 | 0.1755 | 0.825 | -- | -- | -- | validated |
| RPAIN | scz | Cortex | caudate | secondary | direct | 229 | 0.0694 | 0.069 | -- | -- | -- | validated |
| RPAIN | scz | Cortex | caudate | secondary | direct | 219 | 0.1932 | 0.807 | -- | -- | -- | validated |
| RPAIN | scz | Cortex | caudate | secondary | direct | 164 | 0.4273 | 0.427 | -- | -- | -- | validated |
| RPAIN | scz | Cortex | dlpfc | exact | direct | 222 | 0.2549 | 0.745 | -- | -- | -- | validated |
| RPAIN | scz | Cortex | dlpfc | exact | direct | 222 | 0.1099 | 0.890 | -- | -- | -- | validated |
| RPAIN | scz | Cortex | dlpfc | exact | direct | 185 | 0.2558 | 0.744 | -- | -- | -- | validated |
| RPAIN | scz | Cortex | dlpfc | exact | direct | 185 | 0.3078 | 0.308 | -- | -- | -- | validated |
| RPAIN | scz | Cortex | dlpfc | exact | direct | 178 | 0.4728 | 0.527 | -- | -- | -- | validated |
| RPAIN | scz | Cortex | dlpfc | exact | direct | 138 | 0.4476 | 0.448 | -- | -- | -- | validated |
| RPAIN | scz | Cortex | dlpfc | exact | direct | 222 | 0.1839 | 0.816 | -- | -- | -- | validated |
| RPAIN | scz | Cortex | dlpfc | exact | direct | 222 | 0.2730 | 0.727 | -- | -- | -- | validated |
| RPAIN | scz | Cortex | dlpfc | exact | direct | 222 | 0.3301 | 0.330 | -- | -- | -- | validated |
| RPAIN | scz | Cortex | dlpfc | exact | direct | 222 | 0.1099 | 0.890 | -- | -- | -- | validated |
| RPAIN | scz | Cortex | dlpfc | exact | direct | 214 | 0.0721 | 0.072 | -- | -- | -- | validated |
| RPAIN | scz | Cortex | dlpfc | exact | direct | 185 | 0.2558 | 0.744 | -- | -- | -- | validated |
| RPAIN | scz | Cortex | dlpfc | exact | direct | 138 | 0.4476 | 0.448 | -- | -- | -- | validated |
| TMEM175 | pd | Cortex | caudate | secondary | direct | 238 | 0.0609 | 0.939 | -- | -- | -- | validated |
| TMEM175 | pd | Cortex | caudate | secondary | direct | 238 | 0.2374 | 0.237 | -- | -- | -- | validated |
| TMEM175 | pd | Cortex | caudate | secondary | direct | 236 | 0.2211 | 0.779 | -- | -- | -- | validated |
| TMEM175 | pd | Cortex | dlpfc | exact | direct | 222 | 0.0863 | 0.914 | -- | -- | -- | validated |
| TMEM175 | pd | Cortex | dlpfc | exact | direct | 221 | 0.2185 | 0.218 | -- | -- | -- | validated |

## Events tested

- **B3GAT1** / dlpfc: `AL:chr11:134378504:134380747-134381924:134381645:134381817-134381924:-` (AL), n = 222
- **B3GAT1** / dlpfc: `AL:chr11:134378506:134380747-134381924:134381645:134381817-134381924:-` (AL), n = 222
- **B3GAT1** / dlpfc: `AL:chr11:134378508:134380747-134381924:134381645:134381817-134381924:-` (AL), n = 222
- **B3GAT1** / dlpfc: `AL:chr11:134378504:134380747-134381924:134381645:134381817-134381924:-` (AL), n = 222
- **B3GAT1** / dlpfc: `AL:chr11:134378506:134380747-134381924:134381645:134381817-134381924:-` (AL), n = 222
- **B3GAT1** / dlpfc: `AL:chr11:134378508:134380747-134381924:134381645:134381817-134381924:-` (AL), n = 222
- **GSTO2** / dlpfc: `AF:chr10:104274898:104274949-104277894:104275162:104275334-104277894:+` (AF), n = 222
- **GSTO2** / dlpfc: `AF:chr10:104274920:104274949-104277894:104275162:104275334-104277894:+` (AF), n = 222
- **GSTO2** / dlpfc: `SE:chr10:104274949-104275226:104275334-104277894:+` (SE), n = 222
- **GSTO2** / dlpfc: `SE:chr10:104278116-104279370:104279471-104297578:+` (SE), n = 222
- **GSTO2** / dlpfc: `SE:chr10:104278116-104279370:104279471-104297578:+` (SE), n = 222
- **GSTO2** / dlpfc: `SE:chr10:104278116-104279370:104279471-104297578:+` (SE), n = 222
- **RPAIN** / dlpfc: `AL:chr17:5426299-5428071:5431484:5426299-5432542:5432876:+` (AL), n = 222
- **RPAIN** / dlpfc: `MX:chr17:5422829-5425971:5426082-5428071:5422829-5426236:5426299-5428071:+` (MX), n = 222
- **RPAIN** / dlpfc: `RI:chr17:5425971:5426299-5428071:5428211:+` (RI), n = 222
- **RPAIN** / dlpfc: `SE:chr17:5422829-5426236:5426299-5428071:+` (SE), n = 222
- **RPAIN** / dlpfc: `SE:chr17:5426082-5426236:5426299-5428071:+` (SE), n = 222
- **RPAIN** / dlpfc: `SE:chr17:5426299-5428071:5428211-5432542:+` (SE), n = 222
- **RPAIN** / dlpfc: `AL:chr17:5426299-5428071:5428355:5426299-5432542:5432876:+` (AL), n = 219
- **RPAIN** / dlpfc: `AL:chr17:5426299-5428071:5428653:5426299-5432542:5432876:+` (AL), n = 218
- **RPAIN** / dlpfc: `MX:chr17:5422829-5426236:5426299-5432542:5422829-5428071:5428211-5432542:+` (MX), n = 205
- **RPAIN** / dlpfc: `MX:chr17:5426082-5426236:5426299-5432542:5426082-5428071:5428211-5432542:+` (MX), n = 185
- **RPAIN** / dlpfc: `SE:chr17:5422829-5425971:5426299-5428071:+` (SE), n = 185
- **RPAIN** / dlpfc: `A5:chr17:5426299-5428071:5426082-5428071:+` (A5), n = 178
- **RPAIN** / dlpfc: `MX:chr17:5422829-5425971:5426299-5432542:5422829-5428071:5428211-5432542:+` (MX), n = 138
- **RPAIN** / dlpfc: `A5:chr17:5428211-5432542:5426082-5432542:+` (A5), n = 222
- **RPAIN** / dlpfc: `A5:chr17:5428211-5432542:5426299-5432542:+` (A5), n = 222
- **RPAIN** / dlpfc: `SE:chr17:5422829-5425971:5428211-5432542:+` (SE), n = 222
- **RPAIN** / dlpfc: `SE:chr17:5422829-5428071:5428211-5432542:+` (SE), n = 222
- **RPAIN** / dlpfc: `SE:chr17:5426299-5428071:5428211-5432542:+` (SE), n = 222
- **RPAIN** / dlpfc: `SE:chr17:5426082-5428071:5428211-5432542:+` (SE), n = 214
- **RPAIN** / dlpfc: `MX:chr17:5422829-5425971:5426082-5432542:5422829-5428071:5428211-5432542:+` (MX), n = 212
- **RPAIN** / dlpfc: `MX:chr17:5422829-5426236:5426299-5432542:5422829-5428071:5428211-5432542:+` (MX), n = 205
- **RPAIN** / dlpfc: `MX:chr17:5426082-5426236:5426299-5432542:5426082-5428071:5428211-5432542:+` (MX), n = 185
- **RPAIN** / dlpfc: `MX:chr17:5422829-5425971:5426299-5432542:5422829-5428071:5428211-5432542:+` (MX), n = 138
- **TMEM175** / dlpfc: `A3:chr4:932540-947709:932540-948116:+` (A3), n = 222
- **TMEM175** / dlpfc: `SE:chr4:932540-945999:946126-947709:+` (SE), n = 222
- **TMEM175** / dlpfc: `SE:chr4:932540-947709:947892-948116:+` (SE), n = 222
- **TMEM175** / dlpfc: `SE:chr4:932540-947709:948154-950421:+` (SE), n = 221
- **TMEM175** / dlpfc: `SE:chr4:932540-947709:947892-950421:+` (SE), n = 220
- **TMEM175** / dlpfc: `MX:chr4:932540-947709:947892-950421:932540-948116:948154-950421:+` (MX), n = 129
- **VAMP2** / dlpfc: `A5:chr17:8162369-8162444:8162369-8162878:-` (A5), n = 222
- **VAMP2** / dlpfc: `AF:chr17:8162369-8162878:8162948:8162369-8163467:8163546:-` (AF), n = 222
- **B3GAT1** / caudate: `AL:chr11:134378504:134380747-134381924:134381645:134381817-134381924:-` (AL), n = 238
- **B3GAT1** / caudate: `AL:chr11:134378506:134380747-134381924:134381645:134381817-134381924:-` (AL), n = 238
- **B3GAT1** / caudate: `AL:chr11:134378508:134380747-134381924:134381645:134381817-134381924:-` (AL), n = 238
- **B3GAT1** / caudate: `AL:chr11:134378504:134380747-134381924:134381645:134381817-134381924:-` (AL), n = 238
- **B3GAT1** / caudate: `AL:chr11:134378506:134380747-134381924:134381645:134381817-134381924:-` (AL), n = 238
- **B3GAT1** / caudate: `AL:chr11:134378508:134380747-134381924:134381645:134381817-134381924:-` (AL), n = 238
- **GSTO2** / caudate: `AF:chr10:104274898:104274949-104277894:104275162:104275334-104277894:+` (AF), n = 238
- **GSTO2** / caudate: `AF:chr10:104274920:104274949-104277894:104275162:104275334-104277894:+` (AF), n = 238
- **GSTO2** / caudate: `SE:chr10:104274949-104275226:104275334-104277894:+` (SE), n = 238
- **GSTO2** / caudate: `SE:chr10:104278116-104279370:104279471-104297578:+` (SE), n = 238
- **GSTO2** / caudate: `SE:chr10:104278116-104279370:104279471-104297578:+` (SE), n = 238
- **GSTO2** / caudate: `SE:chr10:104278116-104279370:104279471-104297578:+` (SE), n = 238
- **PIGQ** / caudate: `A5:chr16:580972-582248:580263-582248:+` (A5), n = 238
- **PIGQ** / caudate: `SE:chr16:580972-581199:581341-582248:+` (SE), n = 238
- **PIGQ** / caudate: `SE:chr16:580972-582248:582309-582883:+` (SE), n = 238
- **PIGQ** / caudate: `SE:chr16:580263-580858:580972-582248:+` (SE), n = 235
- **PIGQ** / caudate: `MX:chr16:580263-580858:580972-582883:580263-582248:582309-582883:+` (MX), n = 214
- **PIGQ** / caudate: `MX:chr16:579180-580183:580263-582248:579180-580858:580972-582248:+` (MX), n = 117
- **PIGQ** / caudate: `SE:chr16:580972-582248:582309-582883:+` (SE), n = 238
- **PIGQ** / caudate: `MX:chr16:580263-580858:580972-582883:580263-582248:582309-582883:+` (MX), n = 214
- **PIGQ** / caudate: `AF:chr16:581761:582309-582883:582490:582611-582883:+` (AF), n = 238
- **PIGQ** / caudate: `AL:chr16:582309-582484:582647:582309-582883:583095:+` (AL), n = 238
- **PIGQ** / caudate: `AL:chr16:582309-582484:582647:582309-582883:584108:+` (AL), n = 238
- **PIGQ** / caudate: `SE:chr16:580972-582248:582309-582883:+` (SE), n = 238
- **PIGQ** / caudate: `AL:chr16:582309-582484:582647:582309-582883:584109:+` (AL), n = 232
- **PIGQ** / caudate: `AL:chr16:582309-582484:582647:582309-582883:584083:+` (AL), n = 222
- **PIGQ** / caudate: `MX:chr16:580263-580858:580972-582883:580263-582248:582309-582883:+` (MX), n = 214
- **PIGQ** / caudate: `AL:chr16:582309-582484:582647:582309-582883:583223:+` (AL), n = 212
- **PIGQ** / caudate: `AL:chr16:582309-582484:582647:582309-582883:583334:+` (AL), n = 204
- **PIGQ** / caudate: `AL:chr16:582309-582484:582647:582309-582883:582893:+` (AL), n = 201
- **PIGQ** / caudate: `AL:chr16:582309-582484:582647:582309-582883:584094:+` (AL), n = 198
- **PIGQ** / caudate: `AL:chr16:582309-582484:582647:582309-582883:583035:+` (AL), n = 195
- **PIGQ** / caudate: `AL:chr16:582309-582484:582647:582309-582883:583455:+` (AL), n = 195
- **RPAIN** / caudate: `AL:chr17:5426299-5428071:5431484:5426299-5432542:5432876:+` (AL), n = 238
- **RPAIN** / caudate: `MX:chr17:5422829-5425971:5426082-5428071:5422829-5426236:5426299-5428071:+` (MX), n = 238
- **RPAIN** / caudate: `RI:chr17:5425971:5426299-5428071:5428211:+` (RI), n = 238
- **RPAIN** / caudate: `SE:chr17:5422829-5426236:5426299-5428071:+` (SE), n = 238
- **RPAIN** / caudate: `SE:chr17:5426082-5426236:5426299-5428071:+` (SE), n = 238
- **RPAIN** / caudate: `SE:chr17:5426299-5428071:5428211-5432542:+` (SE), n = 238
- **RPAIN** / caudate: `AL:chr17:5426299-5428071:5428355:5426299-5432542:5432876:+` (AL), n = 236
- **RPAIN** / caudate: `AL:chr17:5426299-5428071:5428653:5426299-5432542:5432876:+` (AL), n = 236
- **RPAIN** / caudate: `MX:chr17:5422829-5426236:5426299-5432542:5422829-5428071:5428211-5432542:+` (MX), n = 228
- **RPAIN** / caudate: `MX:chr17:5426082-5426236:5426299-5432542:5426082-5428071:5428211-5432542:+` (MX), n = 219
- **RPAIN** / caudate: `A5:chr17:5426299-5428071:5426082-5428071:+` (A5), n = 190
- **RPAIN** / caudate: `SE:chr17:5422829-5425971:5426299-5428071:+` (SE), n = 189
- **RPAIN** / caudate: `MX:chr17:5422829-5425971:5426299-5432542:5422829-5428071:5428211-5432542:+` (MX), n = 164
- **RPAIN** / caudate: `A5:chr17:5428211-5432542:5426082-5432542:+` (A5), n = 238
- **RPAIN** / caudate: `A5:chr17:5428211-5432542:5426299-5432542:+` (A5), n = 238
- **RPAIN** / caudate: `SE:chr17:5422829-5425971:5428211-5432542:+` (SE), n = 238
- **RPAIN** / caudate: `SE:chr17:5422829-5428071:5428211-5432542:+` (SE), n = 238
- **RPAIN** / caudate: `SE:chr17:5426299-5428071:5428211-5432542:+` (SE), n = 238
- **RPAIN** / caudate: `MX:chr17:5422829-5425971:5426082-5432542:5422829-5428071:5428211-5432542:+` (MX), n = 229
- **RPAIN** / caudate: `SE:chr17:5426082-5428071:5428211-5432542:+` (SE), n = 229
- **RPAIN** / caudate: `MX:chr17:5422829-5426236:5426299-5432542:5422829-5428071:5428211-5432542:+` (MX), n = 228
- **RPAIN** / caudate: `MX:chr17:5426082-5426236:5426299-5432542:5426082-5428071:5428211-5432542:+` (MX), n = 219
- **RPAIN** / caudate: `MX:chr17:5422829-5425971:5426299-5432542:5422829-5428071:5428211-5432542:+` (MX), n = 164
- **TMEM175** / caudate: `SE:chr4:932540-945999:946126-947709:+` (SE), n = 238
- **TMEM175** / caudate: `SE:chr4:932540-947709:947892-948116:+` (SE), n = 238
- **TMEM175** / caudate: `SE:chr4:932540-947709:947892-950421:+` (SE), n = 238
- **TMEM175** / caudate: `SE:chr4:932540-947709:948154-950421:+` (SE), n = 238
- **TMEM175** / caudate: `A3:chr4:932540-947709:932540-948116:+` (A3), n = 236
- **TMEM175** / caudate: `MX:chr4:932540-947709:947892-950421:932540-948116:948154-950421:+` (MX), n = 156
- **VAMP2** / caudate: `A5:chr17:8162369-8162444:8162369-8162878:-` (A5), n = 238
- **VAMP2** / caudate: `AF:chr17:8162369-8162878:8162948:8162369-8163467:8163546:-` (AF), n = 238

## Display-item contrast (the row the decision rests on)

No event contrasts the anchored junction against the arm the figure claims; the decision rule cannot be applied.

## Pre-registered decision rule

Applies only to genes with a display-item contrast in `_REFERENCE_ARM`; every other gene is reported above, not decided.

Decision rule not applicable -- no display-item contrast measured.

Per-target reasons:

- B3GAT1 / dlpfc: minor-form usage 0.4313 >= threshold 0.05
- B3GAT1 / dlpfc: minor-form usage 0.1042 >= threshold 0.05
- B3GAT1 / dlpfc: minor-form usage 0.2089 >= threshold 0.05
- B3GAT1 / dlpfc: minor-form usage 0.4313 >= threshold 0.05
- B3GAT1 / dlpfc: minor-form usage 0.1042 >= threshold 0.05
- B3GAT1 / dlpfc: minor-form usage 0.2089 >= threshold 0.05
- GSTO2 / dlpfc: minor-form usage 0.4504 >= threshold 0.05
- GSTO2 / dlpfc: minor-form usage 0.1797 >= threshold 0.05
- GSTO2 / dlpfc: minor-form usage 0.2696 >= threshold 0.05
- GSTO2 / dlpfc: minor-form usage 0.4318 >= threshold 0.05
- GSTO2 / dlpfc: minor-form usage 0.4318 >= threshold 0.05
- GSTO2 / dlpfc: minor-form usage 0.4318 >= threshold 0.05
- PPP6R2 / dlpfc: gene is quantified (5 PSI events) but none carries the anchored junction; a coverage gap in the LIBD event catalogue, not evidence against it
- PRDM2 / dlpfc: gene is quantified (19 PSI events) but none carries the anchored junction; a coverage gap in the LIBD event catalogue, not evidence against it
- RPAIN / dlpfc: minor-form usage 0.2549 >= threshold 0.05
- RPAIN / dlpfc: minor-form usage 0.0179 < threshold 0.05
- RPAIN / dlpfc: minor-form usage 0.0000 < threshold 0.05
- RPAIN / dlpfc: minor-form usage 0.0431 < threshold 0.05
- RPAIN / dlpfc: minor-form usage 0.0058 < threshold 0.05
- RPAIN / dlpfc: minor-form usage 0.1099 >= threshold 0.05
- RPAIN / dlpfc: minor-form usage 0.0000 < threshold 0.05
- RPAIN / dlpfc: minor-form usage 0.0000 < threshold 0.05
- RPAIN / dlpfc: minor-form usage 0.0000 < threshold 0.05
- RPAIN / dlpfc: minor-form usage 0.2558 >= threshold 0.05
- RPAIN / dlpfc: minor-form usage 0.3078 >= threshold 0.05
- RPAIN / dlpfc: minor-form usage 0.4728 >= threshold 0.05
- RPAIN / dlpfc: minor-form usage 0.4476 >= threshold 0.05
- RPAIN / dlpfc: minor-form usage 0.1839 >= threshold 0.05
- RPAIN / dlpfc: minor-form usage 0.2730 >= threshold 0.05
- RPAIN / dlpfc: minor-form usage 0.3301 >= threshold 0.05
- RPAIN / dlpfc: minor-form usage 0.0000 < threshold 0.05
- RPAIN / dlpfc: minor-form usage 0.1099 >= threshold 0.05
- RPAIN / dlpfc: minor-form usage 0.0721 >= threshold 0.05
- RPAIN / dlpfc: minor-form usage 0.0000 < threshold 0.05
- RPAIN / dlpfc: minor-form usage 0.0000 < threshold 0.05
- RPAIN / dlpfc: minor-form usage 0.2558 >= threshold 0.05
- RPAIN / dlpfc: minor-form usage 0.4476 >= threshold 0.05
- TBC1D15 / dlpfc: gene is quantified (24 PSI events) but none carries the anchored junction; a coverage gap in the LIBD event catalogue, not evidence against it
- TMEM175 / dlpfc: minor-form usage 0.0863 >= threshold 0.05
- TMEM175 / dlpfc: minor-form usage 0.0269 < threshold 0.05
- TMEM175 / dlpfc: minor-form usage 0.0265 < threshold 0.05
- TMEM175 / dlpfc: minor-form usage 0.2185 >= threshold 0.05
- TMEM175 / dlpfc: minor-form usage 0.0000 < threshold 0.05
- TMEM175 / dlpfc: minor-form usage 0.0000 < threshold 0.05
- VAMP2 / dlpfc: minor-form usage 0.0114 < threshold 0.05
- VAMP2 / dlpfc: minor-form usage 0.0165 < threshold 0.05
- B3GAT1 / caudate: minor-form usage 0.4521 >= threshold 0.05
- B3GAT1 / caudate: minor-form usage 0.0872 >= threshold 0.05
- B3GAT1 / caudate: minor-form usage 0.3362 >= threshold 0.05
- B3GAT1 / caudate: minor-form usage 0.4521 >= threshold 0.05
- B3GAT1 / caudate: minor-form usage 0.0872 >= threshold 0.05
- B3GAT1 / caudate: minor-form usage 0.3362 >= threshold 0.05
- GSTO2 / caudate: minor-form usage 0.4022 >= threshold 0.05
- GSTO2 / caudate: minor-form usage 0.2166 >= threshold 0.05
- GSTO2 / caudate: minor-form usage 0.2832 >= threshold 0.05
- GSTO2 / caudate: minor-form usage 0.3699 >= threshold 0.05
- GSTO2 / caudate: minor-form usage 0.3699 >= threshold 0.05
- GSTO2 / caudate: minor-form usage 0.3699 >= threshold 0.05
- PIGQ / caudate: minor-form usage 0.0000 < threshold 0.05
- PIGQ / caudate: minor-form usage 0.1800 >= threshold 0.05
- PIGQ / caudate: minor-form usage 0.1877 >= threshold 0.05
- PIGQ / caudate: minor-form usage 0.0000 < threshold 0.05
- PIGQ / caudate: minor-form usage 0.0000 < threshold 0.05
- PIGQ / caudate: minor-form usage 0.0000 < threshold 0.05
- PIGQ / caudate: minor-form usage 0.1877 >= threshold 0.05
- PIGQ / caudate: minor-form usage 0.0000 < threshold 0.05
- PIGQ / caudate: minor-form usage 0.0000 < threshold 0.05
- PIGQ / caudate: minor-form usage 0.4684 >= threshold 0.05
- PIGQ / caudate: minor-form usage 0.1601 >= threshold 0.05
- PIGQ / caudate: minor-form usage 0.1877 >= threshold 0.05
- PIGQ / caudate: minor-form usage 0.2087 >= threshold 0.05
- PIGQ / caudate: minor-form usage 0.2923 >= threshold 0.05
- PIGQ / caudate: minor-form usage 0.0000 < threshold 0.05
- PIGQ / caudate: minor-form usage 0.0000 < threshold 0.05
- PIGQ / caudate: minor-form usage 0.0000 < threshold 0.05
- PIGQ / caudate: minor-form usage 0.0000 < threshold 0.05
- PIGQ / caudate: minor-form usage 0.0000 < threshold 0.05
- PIGQ / caudate: minor-form usage 0.0000 < threshold 0.05
- PIGQ / caudate: minor-form usage 0.0000 < threshold 0.05
- PPP6R2 / caudate: gene is quantified (5 PSI events) but none carries the anchored junction; a coverage gap in the LIBD event catalogue, not evidence against it
- PRDM2 / caudate: gene is quantified (19 PSI events) but none carries the anchored junction; a coverage gap in the LIBD event catalogue, not evidence against it
- RPAIN / caudate: minor-form usage 0.4727 >= threshold 0.05
- RPAIN / caudate: minor-form usage 0.0312 < threshold 0.05
- RPAIN / caudate: minor-form usage 0.0000 < threshold 0.05
- RPAIN / caudate: minor-form usage 0.0434 < threshold 0.05
- RPAIN / caudate: minor-form usage 0.0081 < threshold 0.05
- RPAIN / caudate: minor-form usage 0.1755 >= threshold 0.05
- RPAIN / caudate: minor-form usage 0.0000 < threshold 0.05
- RPAIN / caudate: minor-form usage 0.0000 < threshold 0.05
- RPAIN / caudate: minor-form usage 0.0000 < threshold 0.05
- RPAIN / caudate: minor-form usage 0.1932 >= threshold 0.05
- RPAIN / caudate: minor-form usage 0.3469 >= threshold 0.05
- RPAIN / caudate: minor-form usage 0.3533 >= threshold 0.05
- RPAIN / caudate: minor-form usage 0.4273 >= threshold 0.05
- RPAIN / caudate: minor-form usage 0.3589 >= threshold 0.05
- RPAIN / caudate: minor-form usage 0.4787 >= threshold 0.05
- RPAIN / caudate: minor-form usage 0.1942 >= threshold 0.05
- RPAIN / caudate: minor-form usage 0.0000 < threshold 0.05
- RPAIN / caudate: minor-form usage 0.1755 >= threshold 0.05
- RPAIN / caudate: minor-form usage 0.0000 < threshold 0.05
- RPAIN / caudate: minor-form usage 0.0694 >= threshold 0.05
- RPAIN / caudate: minor-form usage 0.0000 < threshold 0.05
- RPAIN / caudate: minor-form usage 0.1932 >= threshold 0.05
- RPAIN / caudate: minor-form usage 0.4273 >= threshold 0.05
- TBC1D15 / caudate: gene is quantified (24 PSI events) but none carries the anchored junction; a coverage gap in the LIBD event catalogue, not evidence against it
- TMEM175 / caudate: minor-form usage 0.0270 < threshold 0.05
- TMEM175 / caudate: minor-form usage 0.0609 >= threshold 0.05
- TMEM175 / caudate: minor-form usage 0.0000 < threshold 0.05
- TMEM175 / caudate: minor-form usage 0.2374 >= threshold 0.05
- TMEM175 / caudate: minor-form usage 0.2211 >= threshold 0.05
- TMEM175 / caudate: minor-form usage 0.0000 < threshold 0.05
- VAMP2 / caudate: minor-form usage 0.0137 < threshold 0.05
- VAMP2 / caudate: minor-form usage 0.0180 < threshold 0.05

- 60 targets over 22 genes sit in a tissue BrainSEQ does not sequence: Cerebellum (14 genes), Cerebellar_Hemisphere (8 genes), Amygdala (1 gene), Hypothalamus (1 gene).