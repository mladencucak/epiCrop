# epiCrop
Set of functions and data for crop disease epidemiology work

## What is in it
1. A function to import, check, organise and interpolate the weather data. (under development)
2. Models:  
-BlightR: PLB model  
-IrishRules: PLB model
3. Template weather data as an example:  
-weather: Example weather data needed to run the model.
4. Vignette "run_BlightR" explains how to import data, run the model and visualise the results. (under development)

## How to install
The package is not on CRAN so it can be installed directly from this repository.
1. Restart r sesion (Ctrl/Cmd+Shift+F10).
2. Install `remotes` package.
3. Install the `BligtR` package.

This could be accomplished by running the following code:
``` r
.rs.restartR()
Sys.sleep(2)
if (!"remotes"%in% installed.packages()) {
  install.packages("remotes", repos = "http://cran.rstudio.com/")
  library("remotes")
}
remotes::install_github("mladencucak/epiCrop",dependencies = TRUE)
```
#### If you run into problems
Installing packages that are not on CRAN can be a pain, so there are a  few general notes on what to do.  
If you are Windows user it is recommended that you install [Rtools](https://cran.r-project.org/bin/windows/Rtools/).
Generally there are always packages that are not If you get an error looking like this while updating your packages:
``` r
Error: (converted from warning) cannot remove prior installation of package ‘name_of_the_package’
```
You might need to update some packages manually to match the *source* version. So, just replace the `name_of_the_package` to match the name of the package you are missing/need to update.   
``` r
install.packages("name of the package", type = "source")
```
After try run the code snippet from the beginning of this section.
