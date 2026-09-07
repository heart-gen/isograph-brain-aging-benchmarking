# Short-read junction confirmation of the anchored switch pairs

Targets: **32**  |  validated: **18**  |  alpha = 0.05, minor-form usage threshold = 0.05

`direct` rows test a single PSI event whose two arms contrast the anchored
junction against the alternative form, so PSI is the switch ratio itself and there
is no compositional-closure confound; the statistic is minor-form usage,
`min(median PSI, 1 - median PSI)`, which does not depend on PSI orientation.
`paired` rows test PSI anti-correlation against the closure baseline.

| gene | trait | region | mode | n | minor-form usage | median PSI | rho | p_within | p_matched | verdict |
| --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| CTSH | ad | hippocampus | direct | 238 | 0.0159 | 0.016 | -- | -- | -- | not_validated |
| CTSH | ad | hippocampus | direct | 238 | 0.0000 | 1.000 | -- | -- | -- | not_validated |
| CTSH | ad | hippocampus | direct | 238 | 0.0000 | 0.000 | -- | -- | -- | not_validated |
| CTSH | ad | hippocampus | direct | 238 | 0.0203 | 0.020 | -- | -- | -- | not_validated |
| SNCA | lbd | dlpfc | direct | 222 | 0.1884 | 0.812 | -- | -- | -- | validated |
| SNCA | lbd | dlpfc | direct | 222 | 0.3177 | 0.682 | -- | -- | -- | validated |
| SNCA | lbd | dlpfc | direct | 222 | 0.0209 | 0.979 | -- | -- | -- | not_validated |
| SNCA | lbd | dlpfc | direct | 222 | 0.1988 | 0.801 | -- | -- | -- | validated |
| SNCA | lbd | dlpfc | direct | 222 | 0.1886 | 0.811 | -- | -- | -- | validated |
| SNCA | lbd | dlpfc | direct | 221 | 0.0054 | 0.995 | -- | -- | -- | not_validated |
| SNCA | lbd | dlpfc | direct | 221 | 0.0017 | 0.998 | -- | -- | -- | not_validated |
| SNCA | pd | dlpfc | direct | 222 | 0.1884 | 0.812 | -- | -- | -- | validated |
| SNCA | pd | dlpfc | direct | 222 | 0.3177 | 0.682 | -- | -- | -- | validated |
| SNCA | pd | dlpfc | direct | 222 | 0.0209 | 0.979 | -- | -- | -- | not_validated |
| SNCA | pd | dlpfc | direct | 222 | 0.1988 | 0.801 | -- | -- | -- | validated |
| SNCA | pd | dlpfc | direct | 222 | 0.1886 | 0.811 | -- | -- | -- | validated |
| SNCA | pd | dlpfc | direct | 221 | 0.0054 | 0.995 | -- | -- | -- | not_validated |
| SNCA | pd | dlpfc | direct | 221 | 0.0017 | 0.998 | -- | -- | -- | not_validated |
| SNCA | lbd | caudate | direct | 238 | 0.3050 | 0.695 | -- | -- | -- | validated |
| SNCA | lbd | caudate | direct | 238 | 0.3822 | 0.618 | -- | -- | -- | validated |
| SNCA | lbd | caudate | direct | 238 | 0.0524 | 0.948 | -- | -- | -- | validated |
| SNCA | lbd | caudate | direct | 238 | 0.3386 | 0.661 | -- | -- | -- | validated |
| SNCA | lbd | caudate | direct | 238 | 0.2343 | 0.766 | -- | -- | -- | validated |
| SNCA | lbd | caudate | direct | 236 | 0.0000 | 1.000 | -- | -- | -- | not_validated |
| SNCA | lbd | caudate | direct | 236 | 0.0000 | 1.000 | -- | -- | -- | not_validated |
| SNCA | pd | caudate | direct | 238 | 0.3050 | 0.695 | -- | -- | -- | validated |
| SNCA | pd | caudate | direct | 238 | 0.3822 | 0.618 | -- | -- | -- | validated |
| SNCA | pd | caudate | direct | 238 | 0.0524 | 0.948 | -- | -- | -- | validated |
| SNCA | pd | caudate | direct | 238 | 0.3386 | 0.661 | -- | -- | -- | validated |
| SNCA | pd | caudate | direct | 238 | 0.2343 | 0.766 | -- | -- | -- | validated |
| SNCA | pd | caudate | direct | 236 | 0.0000 | 1.000 | -- | -- | -- | not_validated |
| SNCA | pd | caudate | direct | 236 | 0.0000 | 1.000 | -- | -- | -- | not_validated |

## Events tested

- **CTSH** / hippocampus: `A3:chr15:78937496-78939140:78937423-78939140:-` (A3), n = 238
- **CTSH** / hippocampus: `SE:chr15:78935750-78937318:78937423-78939140:-` (SE), n = 238
- **CTSH** / hippocampus: `SE:chr15:78937423-78937686:78937792-78939140:-` (SE), n = 238
- **CTSH** / hippocampus: `SE:chr15:78937423-78937686:78937823-78939140:-` (SE), n = 238
- **SNCA** / dlpfc: `AF:chr4:89835692-89836127:89836213:89835692-89836743:89837139:-` (AF), n = 222
- **SNCA** / dlpfc: `AF:chr4:89835692-89836127:89836213:89835692-89836962:89837161:-` (AF), n = 222
- **SNCA** / dlpfc: `AF:chr4:89835692-89836127:89836213:89835692-89836962:89837310:-` (AF), n = 222
- **SNCA** / dlpfc: `AF:chr4:89835692-89836127:89836213:89835692-89837098:89837168:-` (AF), n = 222
- **SNCA** / dlpfc: `AF:chr4:89835692-89836127:89836213:89835692-89838252:89838315:-` (AF), n = 222
- **SNCA** / dlpfc: `AF:chr4:89835692-89836127:89836213:89835692-89836743:89836976:-` (AF), n = 221
- **SNCA** / dlpfc: `AF:chr4:89835692-89836127:89836213:89835692-89836962:89837199:-` (AF), n = 221
- **SNCA** / dlpfc: `AF:chr4:89835692-89836127:89836213:89835692-89836743:89837139:-` (AF), n = 222
- **SNCA** / dlpfc: `AF:chr4:89835692-89836127:89836213:89835692-89836962:89837161:-` (AF), n = 222
- **SNCA** / dlpfc: `AF:chr4:89835692-89836127:89836213:89835692-89836962:89837310:-` (AF), n = 222
- **SNCA** / dlpfc: `AF:chr4:89835692-89836127:89836213:89835692-89837098:89837168:-` (AF), n = 222
- **SNCA** / dlpfc: `AF:chr4:89835692-89836127:89836213:89835692-89838252:89838315:-` (AF), n = 222
- **SNCA** / dlpfc: `AF:chr4:89835692-89836127:89836213:89835692-89836743:89836976:-` (AF), n = 221
- **SNCA** / dlpfc: `AF:chr4:89835692-89836127:89836213:89835692-89836962:89837199:-` (AF), n = 221
- **SNCA** / caudate: `AF:chr4:89835692-89836127:89836213:89835692-89836743:89837139:-` (AF), n = 238
- **SNCA** / caudate: `AF:chr4:89835692-89836127:89836213:89835692-89836962:89837161:-` (AF), n = 238
- **SNCA** / caudate: `AF:chr4:89835692-89836127:89836213:89835692-89836962:89837310:-` (AF), n = 238
- **SNCA** / caudate: `AF:chr4:89835692-89836127:89836213:89835692-89837098:89837168:-` (AF), n = 238
- **SNCA** / caudate: `AF:chr4:89835692-89836127:89836213:89835692-89838252:89838315:-` (AF), n = 238
- **SNCA** / caudate: `AF:chr4:89835692-89836127:89836213:89835692-89836743:89836976:-` (AF), n = 236
- **SNCA** / caudate: `AF:chr4:89835692-89836127:89836213:89835692-89836962:89837199:-` (AF), n = 236
- **SNCA** / caudate: `AF:chr4:89835692-89836127:89836213:89835692-89836743:89837139:-` (AF), n = 238
- **SNCA** / caudate: `AF:chr4:89835692-89836127:89836213:89835692-89836962:89837161:-` (AF), n = 238
- **SNCA** / caudate: `AF:chr4:89835692-89836127:89836213:89835692-89836962:89837310:-` (AF), n = 238
- **SNCA** / caudate: `AF:chr4:89835692-89836127:89836213:89835692-89837098:89837168:-` (AF), n = 238
- **SNCA** / caudate: `AF:chr4:89835692-89836127:89836213:89835692-89838252:89838315:-` (AF), n = 238
- **SNCA** / caudate: `AF:chr4:89835692-89836127:89836213:89835692-89836743:89836976:-` (AF), n = 236
- **SNCA** / caudate: `AF:chr4:89835692-89836127:89836213:89835692-89836962:89837199:-` (AF), n = 236

## Display-item contrast (the row the decision rests on)

| gene | region | event | n | minor-form usage | verdict |
| --- | --- | --- | ---: | ---: | --- |
| CTSH | hippocampus | `SE:chr15:78937423-78937686:78937792-78939140:-` | 238 | 0.0000 | not_validated |
| CTSH | hippocampus | `SE:chr15:78937423-78937686:78937823-78939140:-` | 238 | 0.0203 | not_validated |
| SNCA | dlpfc | `AF:chr4:89835692-89836127:89836213:89835692-89838252:89838315:-` | 222 | 0.1886 | validated |
| SNCA | caudate | `AF:chr4:89835692-89836127:89836213:89835692-89838252:89838315:-` | 238 | 0.2343 | validated |

## Pre-registered decision rule

- **CTSH does not validate.** Per the pre-registered rule it stays off any main figure and the set-level result stands alone.
- **SNCA validates.** Report the short-read result as the orthogonal confirmation and cite the long-read failure as an assay limitation; SNCA may appear on a main figure.

Per-target reasons:

- CTSH / hippocampus: minor-form usage 0.0159 < threshold 0.05
- CTSH / hippocampus: minor-form usage 0.0000 < threshold 0.05
- CTSH / hippocampus: minor-form usage 0.0000 < threshold 0.05
- CTSH / hippocampus: minor-form usage 0.0203 < threshold 0.05
- SNCA / dlpfc: minor-form usage 0.1884 >= threshold 0.05
- SNCA / dlpfc: minor-form usage 0.3177 >= threshold 0.05
- SNCA / dlpfc: minor-form usage 0.0209 < threshold 0.05
- SNCA / dlpfc: minor-form usage 0.1988 >= threshold 0.05
- SNCA / dlpfc: minor-form usage 0.1886 >= threshold 0.05
- SNCA / dlpfc: minor-form usage 0.0054 < threshold 0.05
- SNCA / dlpfc: minor-form usage 0.0017 < threshold 0.05
- SNCA / dlpfc: minor-form usage 0.1884 >= threshold 0.05
- SNCA / dlpfc: minor-form usage 0.3177 >= threshold 0.05
- SNCA / dlpfc: minor-form usage 0.0209 < threshold 0.05
- SNCA / dlpfc: minor-form usage 0.1988 >= threshold 0.05
- SNCA / dlpfc: minor-form usage 0.1886 >= threshold 0.05
- SNCA / dlpfc: minor-form usage 0.0054 < threshold 0.05
- SNCA / dlpfc: minor-form usage 0.0017 < threshold 0.05
- SNCA / caudate: minor-form usage 0.3050 >= threshold 0.05
- SNCA / caudate: minor-form usage 0.3822 >= threshold 0.05
- SNCA / caudate: minor-form usage 0.0524 >= threshold 0.05
- SNCA / caudate: minor-form usage 0.3386 >= threshold 0.05
- SNCA / caudate: minor-form usage 0.2343 >= threshold 0.05
- SNCA / caudate: minor-form usage 0.0000 < threshold 0.05
- SNCA / caudate: minor-form usage 0.0000 < threshold 0.05
- SNCA / caudate: minor-form usage 0.3050 >= threshold 0.05
- SNCA / caudate: minor-form usage 0.3822 >= threshold 0.05
- SNCA / caudate: minor-form usage 0.0524 >= threshold 0.05
- SNCA / caudate: minor-form usage 0.3386 >= threshold 0.05
- SNCA / caudate: minor-form usage 0.2343 >= threshold 0.05
- SNCA / caudate: minor-form usage 0.0000 < threshold 0.05
- SNCA / caudate: minor-form usage 0.0000 < threshold 0.05