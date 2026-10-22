# Recorded runtimes

The original simulation runs used 30 worker processes. Approximate runtimes in hours are summarized from the original logs. Method and replication coverage are recorded in `runtime_batches.csv`.

## Latest 100-replication nonparametric batches

These batches compute QR1, QR2, TR, and O-TR. Times are hours for one configuration, using 30 workers. Full estimation also computes the parametric methods.

| Setting | N = 600 | N = 1000 | N = 1400 | N = 1800 |
| --- | --- | --- | --- | --- |
| (i) | 1.6 | 3.1 | 6.6 | 12.3 |
| (ii) | 1.1 | 2.1 | 3.2 | 7.4 |
| (iii) | 1.5 | 2.5 | 4.1 | 6.5 |
| (iv) | 0.6 | 1.1 | 1.8 | 2.8 |

| Setting (v), signal | N | Hours |
| --- | --- | --- |
| 0.25 | 1800 | 7.8 |
| 0.5 | 1800 | 7.7 |
| 2 | 1800 | 7.7 |

| Setting (vi), sparsity | N | Hours |
| --- | --- | --- |
| (4,12) | 1800 | 6.4 |
| (7,18) | 1800 | 7.0 |

## Configurations run in separate batches

| Configuration | Methods | Saved replications | Hours |
| --- | --- | --- | --- |
| Setting (v), signal = 1 | QR1, QR2, TR, O-TR | 1-60 | 6.8 |
| Setting (v), signal = 1 | QR1, QR2, TR, O-TR | 60-89 | 2.3 |
| Setting (v), signal = 1 | QR1, QR2, TR, O-TR | 90-100 | 1.3 |
| Setting (vi), (13,30) | MQR1, MQR2, O-TR', TR, O-TR | 1-30 | 1.3 |
| Setting (vi), (13,30) | MQR1, MQR2, O-TR', TR, O-TR | 31-100 | 3.8 |
| Setting (vi), (13,30) | QR1, QR2 | 31-100 | 3.9 |
| Setting (vi), (17,39) | MQR1, MQR2, O-TR', TR, O-TR | 1-30 | 1.5 |
| Setting (vi), (17,39) | MQR1, MQR2, O-TR', TR, O-TR | 31-100 | 4.3 |
| Setting (vi), (17,39) | QR1, QR2 | 1-30 | 1.8 |
| Setting (vi), (17,39) | QR1, QR2 | 31-100 | 4.7 |

For signal = 1, replication 60 is taken once from the first batch; the next contributes 61-89. The runs overlap in time, so their durations are not added as a single wall-clock runtime. For (13,30), the separate QR log for replications 1-30 was not retained.

## ACTG175 analysis

The QR2, TR, and O-TR batch took about 2 hours for 100 sample splits. The real-data implementation specifies 31 worker processes.

## Progress during a new run

Method-by-method progress is written to the log file rather than continuously printed to the R console or terminal. The run directory is printed at startup. Open `internal_raw/logs/` within it (`internal_raw/ld1/` for Setting (iv), `internal_raw/logs_cd/` for real data). When execution finishes, the run directory contains a `Final_*.RData` and a summary CSV.

Regenerating all seven numerical tables from the saved results takes less than one minute. Full estimation requires more than eight hours; the recorded batch timings above provide a guide for planning individual runs.
