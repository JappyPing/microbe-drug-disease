from DRHGCN import parse, train, report
from Models.DRHGCN.DRHGCN.model import DRHGCN
import os

if __name__=="__main__":
    # os.environ["CUDA_VISIBLE_DEVICES"] = "-1"
    dir = './ten_fold_predict_result/'
    if not os.path.exists(dir):
        os.mkdir(dir)
    subdir = './ten_fold_predict_result/DRHGCN'
    if not os.path.exists(subdir):
        os.mkdir(subdir)
    args = parse(print_help=True)
    train(args, DRHGCN)
    # report("runs")