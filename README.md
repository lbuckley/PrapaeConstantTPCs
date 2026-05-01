# P. rapae constant temperature TPC experiment

Repeating a study of the thermal sensitivity of P. rapae caterpillar growth

# GENERAL INFORMATION

This README.txt file was updated on May 1, 2026 by Lauren Buckley

## A. Paper associated with this archive

Citation: Taylor M. Hatcher, J. Gwen Shlichta, Lauren B. Buckley, and Joel G. Kingsolver. 2026. Contemporary evolution of insect thermal sensitivity across timescales and ontogeny

Synopsis: The interplay of ectotherm growth and developmental responses to temperature across timescales and ontogeny shapes fitness responses to climate change and variability. Fitness responses depend on whether temperature sensitivity evolves. We repeat a study of thermal sensitivity of larval growth rates in a population of the Cabbage White Butterfly (Pieris rapae) to assess evolutionary responses to 25 years of climate change. Modern caterpillars grow faster, particularly at warm temperatures that have become more prevalent in recent decades. Growth rate increases are most pronounced at shorter timescales and earlier developmental stages. Our study points to the potential for contemporary evolution consistent with a response to recent climate change even in a widespread and dispersive species.

## B. Originators

Taylor Hatcher and Lauren B. Buckley, Department of Biology, University of Washington, Seattle, WA 98195-1800, USA

## C. Contact information

Lauren Buckley. Department of Biology, University of Washington, Seattle, WA 98195-1800, USA. [lbuckley\@uw.edu](mailto:lbuckley@uw.edu){.email}

## D. Dates of data collection

1999 and 2024.

## E. Geographic Location(s) of data collection

Seattle, WA

## F. Funding Sources

This work was supported by the National Science Foundation (IOS-2222089 to L.B.B., IOS-2222090 to J.G.K.).

# ACCESS INFORMATION

## 1. Licenses/restrictions placed on the data or code

CC0 1.0 Universal (CC0 1.0). Public Domain Dedication

## 2. Data derived from other sources

Data from the 1999 experiment be were provided by Joel Kingsolver and are reported in bioRxiv:

Kingsolver JG, Shlichta JG & Moore ME. 2023. Heat stress and the temporal dynamics of insect growth. bioRxiv: 2023--09.

The data are also published here: <https://datadryad.org/dataset/doi:10.5061/dryad.nzs7h44zz>.

## 3. Recommended citation for this data/code archive

Data: source publications above

Code: Hatcher TM et al. 2026. Contemporary evolution of insect thermal sensitivity across timescales and ontogeny. <https://github.com/lbuckley/PrapaeConstantTPCst>.

# DATA & CODE FILE OVERVIEW

This data repository consist of 3 data files, 2 code scripts, this README document, and an archive code folder containing two scripts. The repository includes the following data and code filenames and variables.

Raw data is in the data folder and figures are written to the figures folder.

## Data files and variables

data/:

PrapaeGardenTemps_WARP.csv: Caterpillar model temperatures

| Column | Definition                       |
|--------|----------------------------------|
| dt     | Decimal day of year              |
| Date   | Day of year                      |
| Time   | Time in HHMM format              |
| hr     | Hour of day                      |
| T      | Temperature channel              |
| value  | Recorded measurement (degrees C) |
| Year   | Calendar year                    |

PastPresentFilteredConstantTpc2024.csv: Short term growth rate data from past and recent TPC constant experiment

| Column   | Definition                                                           |
|--------------|----------------------------------------------------------|
| UniID    | Unique identifier for each individual                                |
| mom      | Maternal identifier                                                  |
| ID       | individual identifier within maternal family                         |
| temp     | Measurement temperature (C)                                          |
| active   | Activity status during the interval (e.g., y = active, n = inactive) |
| instar   | Larval instar at observation                                         |
| time     | Measurement time duration (h)                                        |
| duration | Length of the observation interval (same units as `time`)            |
| Mo       | Initial mass at start of experiment (mg)                             |
| fw       | Final mass at end of interval (mg)                                   |
| time.per | Label for time period (e.g., "current", "previous")                  |
| mgain    | Mass gain during interval (fw − Mo, mg)                              |
| Jdate    | Julian date of observation (day of year)                             |
| rgrlog   | Log10-transformed relative growth rate over the interval             |
| gr       | Untransformed growth rate                                            |

GHCNdata/USW00094290_2025.csv: GHCN weather station data, see GHCN format

## Code scripts and workflow

1.  Figs1-2_PastFiguresAnalysis.R. Analyzes past and environmental data and produces figures 1-2 and supplementary figures.

2.  Figs3-4_CompFiguresAnalysis.R: Analyzes experimental data and produces figures 3-4 and supplementary figues.

Archive code:

DataAssembly.R: Assembles and prepares data for analysis.

TimeDependenceAnalysis.R: Assesses sensitivity to past variation in measurement time.

# SOFTWARE VERSIONS sessionInfo()

R version 4.3.1 (2023-06-16)

Packages:

library(ggplot2) #ggplot2_3.5.0

library(data.table) #data.table_1.14.8

library(patchwork) #patchwork_1.2.0.9000

library(reshape2) #reshape2_1.4.4

library(viridis) #viridis_0.6.4

library(nlme) #nlme_3.1-162

library(lme4) #lme4_1.1-34

library(ggridges) #ggridges_0.5.4

library(lubridate) #lubridate_1.9.2

library(dplyr) #dplyr_1.1.2

library(car) #car_3.1-2

library(tidyverse) #tidyverse_2.0.0

library(sjPlot) #sjPlot_2.9.0

library(zoo) #zoo_1.8-12

library(TrenchR) #TrenchR_1.1.1

library(mgcv) #mgcv_1.8-42

library(coin) #coin_1.4-3

# REFERENCES

Kingsolver JG, Shlichta JG & Moore ME. 2023. Heat stress and the temporal dynamics of insect growth. bioRxiv: 2023--09.
