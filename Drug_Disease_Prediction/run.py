import Models.tools as tools
import os

dir = './ten_fold_predict_result/'
if not os.path.exists(dir):
    os.mkdir(dir)
subdir_list = ['BNNR', 'DR2DI', 'MBiRW', 'MSBMF']
for subdir in subdir_list:
    print(subdir)
    if not os.path.exists(dir + subdir):
        os.mkdir(os.path.join(dir, str(subdir)))

# tools.auc_aupr_ten('BNNR')
# tools.auc_aupr_ten('DR2DI')
# tools.auc_aupr_ten('MBiRW')
# tools.auc_aupr_ten('MSBMF')