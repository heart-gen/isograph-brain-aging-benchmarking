# GPU VAE reproducibility probe

`gpu_reproducibility_probe.py` over 24 runs stratified by scenario. Each run was re-executed from the current grid config into an isolated root — GPU twice, its CPU twin once — and compared with the stored telemetry. Stored runs were not touched and no dataset was regenerated.

IsoGraph versions seen: `{"stored_gpu__isograph_version": ["0.1.2", "0.1.4"], "stored_cpu__isograph_version": ["0.1.2", "0.1.4"], "gpu_rep1__isograph_version": ["0.1.5"], "gpu_rep2__isograph_version": ["0.1.5"], "cpu_rep1__isograph_version": ["0.1.5"]}`.

## How to read it

- **GPU run-to-run** separates nondeterminism from everything else: if rep1 and rep2 differ, re-running cannot make the arm reproducible.
- **Current vs stored, per device** separates code drift from device effects: if the CPU twin also fails to reproduce its stored value, the drift is grid-wide.
- **GPU vs CPU** tests the documented claim that the backends are numerically matched.

| contrast | metric | n | n_identical | max_abs_diff | median_abs_diff |
|---|---|---|---|---|---|
| GPU run-to-run (rep1 vs rep2, current code) | module_recovery | 24 | 24 | 0 | 0 |
| GPU run-to-run (rep1 vs rep2, current code) | ari_planted | 21 | 21 | 0 | 0 |
| GPU run-to-run (rep1 vs rep2, current code) | n_predicted_modules | 24 | 24 | 0 | 0 |
| GPU run-to-run (rep1 vs rep2, current code) | n_edges | 24 | 24 | 0 | 0 |
| GPU run-to-run (rep1 vs rep2, current code) | switch_gene_detection_rate | 24 | 24 | 0 | 0 |
| GPU current code vs stored GPU | module_recovery | 24 | 3 | 0.6466 | 0.032 |
| GPU current code vs stored GPU | ari_planted | 21 | 15 | 0.7668 | 0 |
| GPU current code vs stored GPU | n_predicted_modules | 24 | 0 | 44 | 5 |
| GPU current code vs stored GPU | n_edges | 24 | 0 | 1.034e+04 | 30.5 |
| GPU current code vs stored GPU | switch_gene_detection_rate | 24 | 21 | 0.5 | 0 |
| CPU current code vs stored CPU | module_recovery | 24 | 2 | 0.3812 | 0.02471 |
| CPU current code vs stored CPU | ari_planted | 21 | 15 | 0.3382 | 0 |
| CPU current code vs stored CPU | n_predicted_modules | 24 | 5 | 32 | 2 |
| CPU current code vs stored CPU | n_edges | 24 | 0 | 2.934e+04 | 26 |
| CPU current code vs stored CPU | switch_gene_detection_rate | 24 | 22 | 0.07 | 0 |
| GPU vs CPU, current code | module_recovery | 24 | 17 | 0.2101 | 0 |
| GPU vs CPU, current code | ari_planted | 21 | 18 | 0.2069 | 0 |
| GPU vs CPU, current code | n_predicted_modules | 24 | 7 | 39 | 2.5 |
| GPU vs CPU, current code | n_edges | 24 | 1 | 1.899e+04 | 31.5 |
| GPU vs CPU, current code | switch_gene_detection_rate | 24 | 20 | 1 | 0 |
| GPU vs CPU, stored | module_recovery | 24 | 4 | 0.4755 | 0.00799 |
| GPU vs CPU, stored | ari_planted | 21 | 17 | 0.7436 | 0 |
| GPU vs CPU, stored | n_predicted_modules | 24 | 12 | 10 | 0.5 |
| GPU vs CPU, stored | n_edges | 24 | 2 | 5638 | 14 |
| GPU vs CPU, stored | switch_gene_detection_rate | 24 | 23 | 1 | 0 |
