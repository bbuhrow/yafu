# Typical-quality GNFS polys for test sieving, C150-C185

For each 5-digit size this is one NFS@Home job polynomial (degree 5) of exactly that size, with msieve Murphy E closest to the median of all such jobs within ±1 digit in `../jobs.csv`. These are typical normal-effort polys, not record ones. Each poly comes from the job's "QUEUED AS" post in the mersenneforum queue thread (source URL in each file). All eight pass a check that f(m) ≡ 0 (mod n), with m = −Y0/Y1 mod n.

The files contain only the polynomial. NFS@Home's sieving parameters for each job are kept as `# nfs@home ...` comment lines, so the parameters for a test sieve have to be chosen explicitly:
- **Item 4 in `../local_collection.md` (yafu as it is today):** use the yafu table row (last column below).
- **Item 5 (optional comparison):** uncomment the NFS@Home parameters.

| size | file | Murphy E | NFS@Home lims / lpb / mfb | NFS@Home unique rels, matrix | yafu table row (lims / lpb / mfb / siever / start-q) |
|---|---|---|---|---|---|
| 150 | [c150_C150_829332_3736.poly](c150_C150_829332_3736.poly) | 5.159e-12 | 67M/67M, 29/29, 58/58 | 59M, 2.5M | 25M/25M 29/29 58/58 14e q0 12.5M |
| 155 | [c155_28887_226.poly](c155_28887_226.poly) | 2.807e-12 | 33M/33M, 29/29, 58/84 | 35M, 4.6M | 32M/32M 29/29 58/58 14e q0 16M |
| 160 | [c160_C160_79_97_30.poly](c160_C160_79_97_30.poly) | 1.506e-12 | 67M/67M, 30/30, 60/60 | 118M, 4.0M | 40M/40M 30/30 60/60 14e q0 20M |
| 165 | [c165_C165_159978_12955.poly](c165_C165_159978_12955.poly) | 7.822e-13 | 67M/50M, 31/32, 60/61 | 200M, 6.1M | 49M/49M 30/30 60/60 14e q0 24.5M |
| 170 | [c170_C170_2381080_10275.poly](c170_C170_2381080_10275.poly) | 3.958e-13 | 67M/67M, 30/30, 60/60 | n/a, 5.3M | 59M/59M 31/31 62/62 14e q0 29.5M |
| 175 | [c175_C175_933436_12624.poly](c175_C175_933436_12624.poly) | 1.823e-13 | 134M/134M, 31/31, 62/62 | n/a, 9.2M | 70M/70M 31/31 62/62 15e q0 35M |
| 180 | [c180_C180_933436_12589.poly](c180_C180_933436_12589.poly) | 8.710e-14 | 134M/134M, 31/31, 62/62 | n/a, 12.3M | 82M/82M 31/31 62/62 15e q0 41M |
| 185 | [c185_C185_842592_8066.poly](c185_C185_842592_8066.poly) | 4.150e-14 | 180M/134M, 32/32, 64/94 | n/a, 29.3M | 100M/100M 32/32 64/64 16e q0 50M |

Notes
- The NFS@Home unique relation and matrix columns are what that particular job used. NFS@Home usually oversieves, and target density varies from job to job, so treat these as context, not as the relations needed (see `../summary.md` §1). "n/a" means the job's log or forum report had no unique count.
- The C170 pick ran at lpb 30/30, below both yafu's 31/31 and the typical NFS@Home 31/32 at that size. The poly is still typical. Only that job's parameters were unusual.
- The yafu table rows are `ggnfs_table_Gimarel` 150-185 in `factor/nfs/nfs.c`; check that `start-q` and the siever match what yafu actually picks, since yafu may interpolate between rows.
- Median E per size (±1 digit): 150: 5.07e-12, 155: 2.81e-12, 160: 1.51e-12, 165: 7.82e-13, 170: 3.88e-13, 175: 1.86e-13, 180: 8.71e-14, 185: 4.14e-14. These are below the linear fit in `../summary.md` §3 at 175+, where the fit is pulled up by smaller jobs.
