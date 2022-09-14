import scanpy as sc
from sctriangulate import *
from sctriangulate.preprocessing import *
from sctriangulate.colors import *
import argparse


def main(args):
    # ungroup the args
    adata_path = args.adata_path
    dir_path = args.dir_path
    query = args.query
    predict_doublet = args.predict_doublet
    # main program
    sctriangulate_setting(backend='Agg')
    adata = sc.read(adata_path)
    sctri = ScTriangulate(dir=dir_path,adata=adata,query=query,predict_doublet=predict_doublet)
    sctri.lazy_run()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run ScTriangulate main functionalities in docker container')
    parser.add_argument('--adata_path',type=str,default=None,help='path to the h5ad file')
    parser.add_argument('--dir_path',type=str,default=None,help='path to the result folder, does not need to be present')
    parser.add_argument('--query',nargs='+',default=None,help='specifiy annotations to triangulate, seperated by whitespace')
    parser.add_argument('--predict_doublet',action='store_true')
    args = parser.parse_args()
    main(args)
