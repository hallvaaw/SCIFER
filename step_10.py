import glob
import pandas as pd

reads=glob.glob('*_plus.csv')

for f in reads:
    millmap=f.replace("_plus.csv","_millmap.csv")
    df_1=pd.read_csv(f,delimiter=',',header=None)
    df_2=pd.read_csv(millmap,delimiter=',',header=None)
    df_expand=pd.concat([df_2]*305, ignore_index=True)
    mappedreads6=df_expand*6
    FPKM=df_1/mappedreads6
    j=f+"_FPKM.csv"
    FPKM.to_csv(j, header=None)
    reads=glob.glob('*_minus.csv')

for f in reads:
    millmap=f.replace("_minus.csv","_millmap.csv")
    df_1=pd.read_csv(f,delimiter=',',header=None)
    df_2=pd.read_csv(millmap,delimiter=',',header=None)
    df_expand=pd.concat([df_2]*305, ignore_index=True)
    mappedreads6=df_expand*6
    FPKM=df_1/mappedreads6
    j=f+"_FPKM.csv"
    FPKM.to_csv(j, header=None)
