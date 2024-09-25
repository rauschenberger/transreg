
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/rauschenberger/transreg?svg=true)](https://ci.appveyor.com/project/rauschenberger/transreg)
[![R-CMD-check](https://github.com/rauschenberger/transreg/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/rauschenberger/transreg/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/rauschenberger/transreg/graph/badge.svg)](https://app.codecov.io/gh/rauschenberger/transreg)

# Penalised regression with multiple sets of prior effects

Improves the predictive performance of ridge and lasso regression exploiting one or more sources of prior information on the importance and direction of effects ("transfer learning").

## Installation

Install the current release from
[CRAN](https://CRAN.R-project.org/package=transreg):

``` r
install.packages("transreg")
```

or the latest development version from [GitHub](https://github.com/lcsb-bds/transreg) or [GitLab](https://gitlab.lcsb.uni.lu/bds/transreg):

``` r
#install.packages("remotes")
remotes::install_github("lcsb-bds/transreg") # upstream
remotes::install_github("rauschenberger/transreg") # fork
remotes::install_gitlab("bds/transreg",host="gitlab.lcsb.uni.lu") # mirror
```

## Reference

Armin Rauschenberger 
[![AR](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0000-0001-6498-4801),
Zied Landoulsi
[![ZL](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0000-0002-2327-3904),
Mark A. van de Wiel 
[![MvdW](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0000-0003-4780-8472),
and Enrico Glaab
[![EG](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0000-0003-3977-7469) (2023).
"Penalized regression with multiple sets of prior effects".
*Bioinformatics* 39(12):btad680. [doi: 10.1093/bioinformatics/btad680](https://doi.org/10.1093/bioinformatics/btad680).

## Reproducibility

The code for reproducing the simulations and applications shown in the manuscript is available in a vignette (<https://lcsb-bds.github.io/transreg/articles/analysis.html>). After installing the package with `remotes::install_github("lcsb-bds/transreg",build_vignettes=TRUE)` and restarting R, the vignette can also be loaded with `vignette(topic="analysis",package="transreg")`.

[![CRAN version](https://www.r-pkg.org/badges/version/transreg)](https://CRAN.R-project.org/package=transreg)
[![CRAN RStudio mirror downloads](https://cranlogs.r-pkg.org/badges/transreg)](https://CRAN.R-project.org/package=transreg)
[![Total CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/transreg)](https://CRAN.R-project.org/package=transreg)

## Disclaimer

The R package `transreg` implements penalised regression with multiple sources of prior effects ([Rauschenberger et al., 2023](https://doi.org/10.1093/bioinformatics/btad680)).

Copyright &copy; 2021 Armin Rauschenberger, University of Luxembourg, Luxembourg Centre for Systems Biomedicine (LCSB), Biomedical Data Science (BDS)

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.
