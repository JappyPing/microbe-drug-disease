# Description

Benchmarking the methods of DRHGCN, BNNR, DR2DI, MBiRW and MSBMF for the Drug-Disease Associations (DgDsAs) prediction.

# DRHGCN
## Dependencies
* python == 3.6.12
* pytorch == 1.6.0
* pytorch-lightning==1.0.8
* torch_geometric
* scikit-learn 
* seaborn
* matplotlib
* openpyxl
* xlrd
## Usage
```
cd Durg_Disease_prediciton
python run_DRHGCN.py --dataset_name Fdataset --epochs 400
```
## Reference
Cai, L., Lu, C., Xu, J., Meng, Y., Wang, P., Fu, X., ... & Su, Y. (2021). Drug repositioning based on the heterogeneous information fusion graph convolutional network. Briefings in Bioinformatics, 22(6), bbab319.

# BNNR, DR2DI, MBiRW and MSBMF
## Dependencies
matlab 2016 or later
## Usage
```
cd Drug_Disease_prediciton/Models
```
Obtain and save the prediction results
```
matlab -nodesktop -nosplash -r ten_times_ten_fold_CV # using matlab in linux command
# or run 'ten_times_ten_fold_CV.m' in GUI
```
AUROC and AUPR calculation
```
cd ..
python run.py
```
## References
[1] Luo, H., Wang, J., Li, M., Luo, J., Peng, X., Wu, F. X., & Pan, Y. (2016). Drug repositioning based on comprehensive similarity measures and Bi-Random walk algorithm. Bioinformatics, 32(17), 2664-2671.

[2] Yang, M., Wu, G., Zhao, Q., Li, Y., & Wang, J. (2021). Computational drug repositioning based on multi-similarities bilinear matrix factorization. Briefings in Bioinformatics, 22(4), bbaa267.

[3] Lu, L., & Yu, H. (2018). DR2DI: a powerful computational tool for predicting novel drug-disease associations. Journal of Computer-Aided Molecular Design, 32(5), 633-642.

[4] Yang, M., Luo, H., Li, Y., & Wang, J. (2019). Drug repositioning based on bounded nuclear norm regularization. Bioinformatics, 35(14), i455-i463.
