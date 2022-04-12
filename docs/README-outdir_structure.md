
# Output directory structure

Below is the structure of the results directory. Pathways will be populated with parameters set in the [config_analysis](../config_analysis.json). An example of a parameters file is found in the [demo params](../demo/config_analysis.json).

```bash
flexible_variance_decomposition
├── covariates.tsv.gz
├── moltraits.tsv.gz
├── covariance_matrices
│   ├── mtx_1.tsv.gz
│   ├── mtx_2.tsv.gz
│   ... etc. ...
├── factor_0
│   ├── covariates.tsv.gz
│   ├── moltraits.tsv.gz
│   ├── norm=norm_method__resid=resid_distr_1
│   │   ├── mtc=effects_to_test_1
│   │   │   ├── var_decomp.tsv.gz
│   │   │   ├── plots
│   │   │   ├── [images: vardec results]
│   │   │   ... etc. ...
│   │   ├── mtc=effects_to_test_n
│   │   ... etc. ...
│   ├── norm=norm_method__resid=resid_distr_2
│   ... etc. ...
├── factor_5
... etc. ...
```
