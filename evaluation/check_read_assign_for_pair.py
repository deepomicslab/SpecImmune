import pandas as pd

# load a csv file into dataframe
def load_csv(file_path):
    return pd.read_csv(file_path)

# given a csv file, and two column name, compare the two columns
def compare_columns(file_path, column1, column2):
    df = load_csv(file_path)
    # print column names
    # print (len(df.columns), df.columns)

    column1_num = 0
    column2_num = 0
    equal_num = 0
    x = []
    y = []
    for i in range(len(df)):
        ### each element is like 316/0.9783, extract the length value and the identity value
        field1 = df[column1][i].split("/")
        field2 = df[column2][i].split("/")
        x.append(float(field1[1]))
        y.append(float(field2[1]))

        if field1[1] > field2[1]:
            print (df['allele'][i], column1, field1[1], field2[1], field1[0], field2[0], sep="\t")
            column1_num += 1
        elif field1[1] < field2[1]:
            print (df['allele'][i], column2, field1[1], field2[1], field1[0], field2[0], sep="\t")
            column2_num += 1
        else:
            # print ("equal")
            equal_num += 1

    print (column1_num, column2_num, equal_num)

    # plot the data
    import matplotlib.pyplot as plt
    import numpy as np
    import seaborn as sns
    plt.figure()
    sns.scatterplot(x=x, y=y, hue=np.where(np.array(x) > np.array(y), 'x large', np.where(np.array(x) < np.array(y), 'y large', 'equal')))
    plt.xlabel(column1)
    plt.ylabel(column2)
    # set x, y lim to [0.9, 1.0]
    plt.xlim(0.85, 1.0)
    plt.ylim(0.85, 1.0)
    plt.plot([0, 1], [0, 1], transform=plt.gca().transAxes, color='red')
    plt.title(f"{column1} vs {column2}")
    plt.savefig(f"/mnt/d/HLAPro_backup/Nanopore_optimize/{column1}_{column2}.pdf")

    ## try to find if the data form one cluster or one cluster using sklearn. If there is two cluster, one cluster is above diagonal, the other is below diagonal
    from sklearn.cluster import KMeans
    import numpy as np
    X = np.array([x, y]).T
    kmeans = KMeans(n_clusters=2, random_state=0).fit(X)
    print (kmeans.labels_)
    print (kmeans.cluster_centers_)
    




# file_path = "/mnt/d/HLAPro_backup/Nanopore_optimize/output6/fredhutch-hla-GO85/fredhutch-hla-GO85.HLA-A.read.matrix.csv"
# column1 = 'HLA-A*01:01:01:01'
# column2 = 'HLA-A*01:01:01:94'

# file_path = "/mnt/d/HLAPro_backup/Nanopore_optimize/output6/fredhutch-hla-1347-4843/fredhutch-hla-1347-4843.HLA-B.read.matrix.csv"
# column1 = 'HLA-B*15:01:01:01'
# column2 = 'HLA-B*44:04'

# file_path = "/mnt/d/HLAPro_backup/Nanopore_optimize/output6/fredhutch-hla-KOSE/fredhutch-hla-KOSE.HLA-B.read.matrix.csv"
# column1 = 'HLA-B*35:01:01:60'
# column2 = 'HLA-B*35:03:01:03'

# file_path = "/mnt/d/HLAPro_backup/Nanopore_optimize/output6/fredhutch-hla-KOSE/fredhutch-hla-KOSE.HLA-C.read.matrix.csv"
# column1 = 'HLA-C*12:03:01:01'
# column2 = 'HLA-C*17:01:01:30'

file_path = "/mnt/d/HLAPro_backup/Nanopore_optimize/output6/fredhutch-hla-KOSE/fredhutch-hla-KOSE.HLA-A.read.matrix.csv"
column1 = 'HLA-A*02:01:01:01'
column2 = 'HLA-A*02:796N'

compare_columns(file_path, column1, column2)
