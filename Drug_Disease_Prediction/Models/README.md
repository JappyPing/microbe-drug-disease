# Description

"drug_disease.mat": dataset
"ten_fold.m": 10-fold cross validation
"ten_times_ten_fold.m": 10-fold cross validation for ten times
"AUC_AUPR_calculate.py": AUC and AUPR calculation

# Usage

1.run ten_times_ten_fold.m
2.run "python demo.py --dataset_name Fdataset --epochs 400"

# References

The folder "src" and the code 'main.py' are from "Cai L, Lu C, Xu J, et al. Drug repositioning based on the heterogeneous information fusion graph convolutional network. Briefings in Bioinformatics 2021"
The code 'normFun.m','setparFun.m','nManiCluester.m','MBiRW.m' are from "Luo H, Wang J, Li M, et al. Drug repositioning based on comprehensive similarity measures and Bi-Random walk algorithm. Bioinformatics 2016; 32:2664?2671"
The code 'MSBMF.m' is from "Yang M, Wu G, Zhao Q, et al. Computational drug repositioning based on multi-similarities bilinear matrix factorization. Brief Bioinform 2021; 22:"
The code 'rls_kron_for_predict.m' is from "Lu L, Yu H. DR2DI: a powerful computational tool for predicting novel drug-disease associations. Journal of Computer-Aided Molecular Design 2018; 32:633?642"
The code 'BNNR.m','svt.m' are from "Yang M, Luo H, Li Y, Wang J. Drug repositioning based on bounded nuclear norm regularization[J]. Bioinformatics, 2019, 35(14): i455-i463."
