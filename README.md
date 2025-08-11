# Estimated Abundance of Northern Pike in the Chatanika River Overwintering Area, 2025

Minto Flats, a large wetland complex located approximately 50 km west of 
Fairbanks, supports the largest sport fishery for northern pike within the 
Arctic-Yukon-Kuskokwim Management Area, as well as a substantial subsistence 
fishery on an overwintering aggregation of fish near the confluence of Goldstream 
Creek and the Chatanika River. Information on the abundance of the population 
present in the Chatanika River Overwintering Area (CROA) during winter is needed 
to ensure sport and subsistence fishery harvests are sustainable. This project 
will estimate the abundance of northern pike ≥600 and ≥720 mm FL for this population.

## Operational Plan

The operational plan is archived here:
https://www.adfg.alaska.gov/FedAidPDFs/ROP.SF.3F.2024.07.pdf

## Folder structure

### /OP_2024

This folder provides materials relevant to the operational plan, which was published 
in 2024.  **/OP_2024/R** and **/OP_2024/R_output** contain fairly minimal R code 
and R-generated tables and figures, respectively; all code and output pertain to
sample size determination.

### /FDS_2025

This folder provides materials relevant to the FDS report written in 2025 after 
project completion.

#### /FDS_2025/R

**1_MintoPike_MR.R**: This file contains most of the analysis code and exploration.
Most of the relevant pieces were copied over into **MR_summary.Rmd**, which is
more current.

**recapr_prep.R**: This file contains many helper functions for Mark-Recapture 
projects, and was copied from another project folder.  The intent is ultimately
to add them to the `recapr` package for R, but they are still incomplete at present.

#### /FDS_2025/flat_data

Raw sampling data are included here for both events, as .csv files.  Note: these 
were created from the .xlsx file provided by the project biologist.

### MR_summary.Rmd 

This was a summary document reporting results of the Mark-Recapture experiment.
The embedded code chunks in this document should be considered the most recent
version of analysis code.
