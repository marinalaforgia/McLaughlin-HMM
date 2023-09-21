# McLaughlin-HMM

Code for estimating and validating estimates of spatial and temporal dispersal in McLaughlin 80 site data. Utilizes code from Pluntz et al. 2018.

* **20221212_mcl_HMM_quadrat.RDS**: Model estimates produced from applying HMM_Functions.R to 80 site data
* **Scripts**: all scripts and other output from models
  * **HMM_Functions.R**: Creation of HMM functions
  * **HMM-McLaughlin.R**: Code for running HMMs on McLaughlin Data
  * **20230718_McLaughlin_HMM.R**: Code for post-HMM analyses
  * **Archive**: old analysis scripts, no longer used
  * **Validation**: model validation code and output (used JetStream2 to run as scripts use a lot of CPUs)
    * **simulation.R**: Code for simulating time series data from HMM estimates
    * **McLaughlin_validation-correlation.R**: code for checking correlation of data vs correlation from estimation procedure (produces Fig S2, corresponds to Pluntz et al. Fig A6); produced
      * **McL-sim-correlation.RDS**: output of correlation using simulated data from 400 patches (max number of patches a species can occur in)
      * **McL-sim-154.RDS**: output of correlation using simulated data from 154 patches (average number of patches a species occurred in)
    * **McLaughlin_validation-species.R**: code for checking accuracy of model estimates per species (produces Fig S3, corresponds to Pluntz et al. Fig A4); produced
      * **McL-sim-species-validation.RDS**: output using 400 patches (max number of patches a species can occur in)
      * **McL-sim-species-validation-50.RDS**: output using 50 patches (min number of patches a species can occur in)
    * **McLaughlin_validation-boot.R**: code for running HMM model 100 times to look at variation in model estimated parameters; produced:
      * **McL-sim-species-boot.RDS**: average and standard deviation of g, c, and s parameters per species; fairly constant estimation except for a handful of species (APHOCC and LEPNIT)
