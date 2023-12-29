import os
import pandas as pd

# Read the files into two dataframes.
def join_standard_csv(foldername,number):
    df1 = pd.read_csv(foldername+'/3drna'+number+'.csv',sep="     ", header=None)
    df1.columns=['3drna']
    extracted_col = df1["3drna"]
    print(extracted_col)
    df2 = pd.read_csv(foldername+'/cgrnasdpc_R1'+number+'.csv',sep="     ", header=None)
    df2.columns=['ID','cgrnaspc']
    df3 = pd.read_csv(foldername+'/cgrnasdppc_R1'+number+'.csv',sep="     ", header=None)
    df3.columns=['ID','cgrnasppc']
    df4 = pd.read_csv(foldername+'/cgrnasp_cn_R1'+number+'.csv',sep="     ", header=None)
    df4.columns=['ID','cgrnasp_cn']
    df5 = pd.read_csv(foldername+'/cgrnasp_R1'+number+'.csv',sep="     ", header=None)
    df5.columns=['ID','cgrnasp']
    df6 = pd.read_csv(foldername+'/dfireR1'+number+'.csv',sep=";", header=None)
    df6.columns=['ID','dfire']
    df7 = pd.read_csv(foldername+'/R1'+number+'_rna3dcnn.csv',sep="	", header=None)
    df7.columns=['ID','rna3dcnn']
    df8 = pd.read_csv(foldername+'/rsrnasp'+number+'.csv',sep="     ", header=None)
    df8.columns=['ID','rsrnasp']
    # Merge the two dataframes, using _ID column as key
    dffinal = pd.merge(df2, df3, on = 'ID')
    dffinal = pd.merge(dffinal, df4, on = 'ID')
    dffinal = pd.merge(dffinal, df5, on = 'ID')
    dffinal = pd.merge(dffinal, df6, on = 'ID')
    dffinal = pd.merge(dffinal, df7, on = 'ID')
    dffinal = pd.merge(dffinal, df8, on = 'ID')
    # dffinal.set_index('ID', inplace = True)

    dffinal.insert(1, "3drna", extracted_col)
    # dffinal = pd.concat([dffinal, extracted_col.rename("3drna")], axis=1)
    # Write it to a new CSV file
    dffinal.to_csv(foldername+'_merged1.csv')
    print(dffinal)

# different approach to merging, helped find an issue
def join_standard_csv2(foldername,number):
    df1 = pd.read_csv(foldername+'/3drna'+number+'.csv',sep="     ", header=None)
    df1.columns=['3drna']
    extracted_col = df1["3drna"]
    print(extracted_col)
    df2 = pd.read_csv(foldername+'/cgrnasdpc_R1'+number+'.csv',sep="     ", header=None)
    df2.columns=['ID','cgrnaspc']
    df3 = pd.read_csv(foldername+'/cgrnasdppc_R1'+number+'.csv',sep="     ", header=None)
    df3.columns=['ID','cgrnasppc']
    extracted_col3 = df3["cgrnasppc"]
    df4 = pd.read_csv(foldername+'/cgrnasp_cn_R1'+number+'.csv',sep="     ", header=None)
    df4.columns=['ID','cgrnasp_cn']
    extracted_col4 = df4["cgrnasp_cn"]
    df5 = pd.read_csv(foldername+'/cgrnasp_R1'+number+'.csv',sep="     ", header=None)
    df5.columns=['ID','cgrnasp']
    extracted_col5 = df5["cgrnasp"]
    df6 = pd.read_csv(foldername+'/dfireR1'+number+'.csv',sep=";", header=None)
    df6.columns=['ID','dfire']
    extracted_col6 = df6["dfire"]
    df7 = pd.read_csv(foldername+'/R1'+number+'_rna3dcnn.csv',sep="	", header=None)
    df7.columns=['ID','rna3dcnn']
    extracted_col7 = df7["rna3dcnn"]
    df8 = pd.read_csv(foldername+'/rsrnasp'+number+'.csv',sep="     ", header=None)
    df8.columns=['ID','rsrnasp']
    extracted_col8 = df8["rsrnasp"]
    df2.insert(2, "3drna", extracted_col)
    df2.insert(2, "cgrnasppc", extracted_col3)
    df2.insert(2, "cgrnasp_cn", extracted_col4)
    df2.insert(2, "cgrnasp", extracted_col5)
    df2.insert(2, "dfire", extracted_col6)
    df2.insert(2, "rna3dcnn", extracted_col7)
    df2.insert(2, "rsrnasp", extracted_col8)

    # Write it to a new CSV file
    df2.to_csv(foldername+'_merged.csv')
    print(df2)

join_standard_csv("r1107","107")
