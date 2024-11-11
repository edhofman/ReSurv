[![R-hub](https://github.com/gpitt71/ReSurv/actions/workflows/rhub.yaml/badge.svg)](https://github.com/gpitt71/ReSurv/actions/workflows/rhub.yaml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10418823.svg)](https://doi.org/10.5281/zenodo.10418823)

# ReSurv

`ReSurv` is an `R` software for predicting IBNR claims. The software includes tools for synthetic **data generation**, **data pre-processing**, **hyperparameters tuning**, **model estimation** and **prediction**.

The package is based on the approach illustrated in *Hiabu M., Hofman E., and Pittarello G. (2023)* and estimates feature dependent development factors using individual reserving data. 


## Available Machine Learning (ML) models

There is a one-to-one relationship between development factors and hazard rates (*Hiabu et al. (2023)*). The package implements extends the following machine learning algorithms for proportional hazard models:

* Cox model with splines (COX, *Gray (1992)*).

* Neural Networks (NN, *Katzman et al. (2018)*).

* eXtreme Gradient Boosting (XGB, *Chen et al. (2016)*).

`ReSurv` extends COX, NN, and XGB to account for ties in left-truncated and right-censored observations.

## Installation

### Developer Version

The developers version of the package can be installed from GitHub.

```
devtools::install_github('https://github.com/edhofman/ReSurv')
```

### Python Dependencies

For using the NN models we suggest to install a virtual environment using

```
install_pyresurv()
```


The default name of the virtual environment is `"pyresurv"`.

We then suggest to refresh the R session and to import the `ReSurv` package in `R` using 

```
library(ReSurv)
reticulate::use_virtualenv('pyresurv')
```

#### Managing Multiple Package Dependencies

This section is taken from the guidelines of the R package [reticulate](https://rstudio.github.io/reticulate/articles/python_dependencies.html) for handling the case of multiple packages in your session that used isolated-package-environments.
The most straightforward solution would be installing a dedicated environment for both.

```
envname <- "./venv"
ReSurv::install_pyresurv(envname = envname)
pysparklyr::install_pyspark(envname = envname)

```



## References 

- *Chen, T., & Guestrin, C. (2016, August). Xgboost: A scalable tree boosting system. In Proceedings of the 22nd acm sigkdd international conference on knowledge discovery and data mining (pp. 785-794).*

- *Gray, R. J. (1992). Flexible methods for analyzing survival data using splines, with applications to breast cancer prognosis. Journal of the American Statistical Association, 87(420), 942-951.*

- *Hiabu, M., Hofman, E., & Pittarello, G. (2023). A machine learning approach based on survival analysis for IBNR frequencies in non-life reserving. arXiv preprint arXiv:2312.14549.* 

- *Snoek, J., Larochelle, H., & Adams, R. P. (2012). Practical bayesian optimization of machine learning algorithms. Advances in neural information processing systems, 25.*

- *Katzman, J. L., Shaham, U., Cloninger, A., Bates, J., Jiang, T., & Kluger, Y. (2018). DeepSurv: personalized treatment recommender system using a Cox proportional hazards deep neural network. BMC medical research methodology, 18, 1-12.*



