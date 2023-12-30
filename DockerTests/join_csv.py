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
    df7mc = pd.read_csv(foldername+'/R1'+number+'_rna3dcnnMC.csv',sep="	", header=None)
    df7mc.columns=['ID','rna3dcnn_mc']
    df8 = pd.read_csv(foldername+'/rsrnasp'+number+'.csv',sep="     ", header=None)
    df8.columns=['ID','rsrnasp']
    # Merge the two dataframes, using _ID column as key
    dffinal = pd.merge(df2, df3, on = 'ID')
    dffinal = pd.merge(dffinal, df4, on = 'ID')
    dffinal = pd.merge(dffinal, df5, on = 'ID')
    dffinal = pd.merge(dffinal, df6, on = 'ID')
    dffinal = pd.merge(dffinal, df7, on = 'ID')
    dffinal = pd.merge(dffinal, df7mc, on = 'ID')
    dffinal = pd.merge(dffinal, df8, on = 'ID')
    # dffinal.set_index('ID', inplace = True)

    dffinal.insert(1, "3drna", extracted_col)
    # dffinal = pd.concat([dffinal, extracted_col.rename("3drna")], axis=1)
    # Write it to a new CSV file
    dffinal=left_join_dataframe(dffinal,"lociparse.csv")
    dffinal.to_csv(foldername+'_merged.csv')
    print(dffinal)


def left_join_csv(destination,source):
    df1 = pd.read_csv(destination)
    df2 = pd.read_csv(source,sep="	", header=None)
    df2.columns=['ID','lociPARSE']
    dffinal = pd.merge(df1, df2, on = 'ID', how='left')
    dffinal.to_csv(destination)
    print(dffinal)
def left_join_dataframe(df,source):

    df2 = pd.read_csv(source,sep="	", header=None)
    df2.columns=['ID','lociPARSE']
    dffinal = pd.merge(df, df2, on = 'ID', how='left')
    # print(dffinal)
    return dffinal

# join_standard_csv("r1107","107")
join_standard_csv("r1108","108")
join_standard_csv("r1116","116")
join_standard_csv("r1117","117")
join_standard_csv("r1126","126")
join_standard_csv("r1128","128")
join_standard_csv("r1136","136")
join_standard_csv("r1138","138")
join_standard_csv("r1149","149")
join_standard_csv("r1156","156")
join_standard_csv("r1189","189")
join_standard_csv("r1190","190")

# left_join_csv("r1107_merged.csv","lociparse.csv")
# left_join_csv("r1108_merged.csv","lociparse.csv")
# left_join_csv("r1116_merged.csv","lociparse.csv")
# left_join_csv("r1117_merged.csv","lociparse.csv")
# left_join_csv("r1126_merged.csv","lociparse.csv")
# left_join_csv("r1128_merged.csv","lociparse.csv")
# left_join_csv("r1136_merged.csv","lociparse.csv")
# left_join_csv("r1138_merged.csv","lociparse.csv")
# left_join_csv("r1149_merged.csv","lociparse.csv")
# left_join_csv("r1156_merged.csv","lociparse.csv")
# left_join_csv("r1189_merged.csv","lociparse.csv")
# left_join_csv("r1190_merged.csv","lociparse.csv")