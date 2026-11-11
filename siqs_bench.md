## SIQS Benchmark Results

Times are wall-clock seconds: **qs\_time** (sieve phase only),
with **total\_time** (including pre/post-processing) in parentheses.
Digit counts refer to the number of decimal digits in the input.

| System | Version | 45d | 50d | 55d | 60d | 65d | 70d | 75d | 80d | 85d | 90d | 95d | 100d |
|:-------|:-------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|
| **11th Gen Intel(R) Core(TM) i5-11600K @ 3.90GHz** <br>Linux (WSL/Ubuntu-22.04)<br>GCC 11.4<br>AVX512IFMA<br>4 threads | 3.1.9 | 0.04<br>(0.09) | 0.05<br>(0.10) | 0.08<br>(0.17) | 0.20<br>(0.38) | 0.52<br>(0.80) | 0.97<br>(1.48) | 3.48<br>(4.34) | 6.86<br>(9.06) | 21.80<br>(25.04) | 96.23<br>(101.46) | | |
| **11th Gen Intel(R) Core(TM) i5-11600K @ 3.90GHz** <br>Windows<br>MSVC 1944 + Clang 19.1.5 <br>AVX512BW<br>4 threads | 3.1.9 | 0.03<br>(0.05) | 0.03<br>(0.06) | 0.06<br>(0.11) | 0.19<br>(0.28) | 0.51<br>(0.65) | 0.95<br>(1.22) | 3.62<br>(4.08) | 6.49<br>(7.57) | 25.42<br>(27.12) | 97.39<br>(100.32) | | |
| **11th Gen Intel(R) Core(TM) i5-11600K @ 3.90GHz** <br>Windows<br>GCC 15.2<br>AVX512BW<br>4 threads | 3.1.9 | 0.03<br>(0.05) | 0.04<br>(0.07) | 0.06<br>(0.12) | 0.17<br>(0.27) | 0.46<br>(0.62) | 0.89<br>(1.19) | 3.25<br>(3.77) | 6.07<br>(7.24) | 19.52<br>(21.04) | 88.43<br>(91.77) | | |
