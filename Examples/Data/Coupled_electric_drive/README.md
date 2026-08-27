# Coupled Electric Drives Dataset

This is a brief explanation of the **Coupled Electric Drives (CED)** benchmark data used by the NL-FRA examples.

The benchmark system consists of two electric motors driving a pulley through a flexible belt. The pulley is restrained by a spring, which introduces a lightly damped dynamic mode. The two electric drives can be controlled independently, allowing the belt tension and belt speed to be controlled simultaneously. The benchmark considered here focuses on the speed-control dynamics. The pulley angular speed is measured using a pulse counter; the sensor is insensitive to the sign of the velocity. The relatively short available data records make the benchmark useful for testing nonlinear system-identification methods.

The data set and reference models were published by Torbjörn Wigren and Maarten Schoukens:

> T. Wigren and M. Schoukens, *Coupled Electric Drives Data Set and Reference Models*, Technical Report 2017-024, Department of Information Technology, Uppsala University, 2017.

## Original Sources

- **Nonlinear Benchmark page:**  
  https://www.nonlinearbenchmark.org/benchmarks/coupled-electric-drives

- **Uppsala University / DiVA record:**  
  https://uu.diva-portal.org/smash/record.jsf?aq2=%5B%5B%5D%5D&c=216&af=%5B%5D&searchType=SIMPLE&sortOrder2=title_sort_asc&query=Wigren+Torbj%C3%B6rn&language=sv&pid=diva2%3A1165531&aq=%5B%5B%5D%5D&sf=all&aqe=%5B%5D&sortOrder=author_sort_asc&onlyFullText=false&noOfRows=50&dswid=8047

- **Original data download:**  
  https://uu.diva-portal.org/smash/get/diva2:1165531/FULLTEXT01.zip

- `DATAUNIF.MAT` / `DATAUNIF.csv` — uniformly distributed-input data records.
- `DATAPRBS.MAT` / `DATAPRBS.csv` — pseudo-random binary sequence (PRBS) input data records.
- `Coupled Electric Drives Data Set and Reference Models.pdf` — the accompanying technical report.

The NL-FRA example [`NARX_band_limited_NOFRF.m`](../../NARX_band_limited_NOFRF.m) uses the uniformly distributed-input data to identify and independently validate a polynomial NARX model before evaluating its NOFRFs using a general band-limited probing input.
