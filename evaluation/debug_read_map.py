import pandas as pd

# Read the CSV file
df = pd.read_csv('/mnt/d/HLAPro_backup/Nanopore_optimize/output5/fredhutch-hla-KT17/fredhutch-hla-KT17.HLA-C.read.matrix.csv', delimiter=',')

# Split the values and convert to integers
df['HLA-C*04:01:01:01_value'] = df['HLA-C*04:01:01:01'].apply(lambda x: int(x.split('/')[0]))
df['HLA-C*04:01:01:11_value'] = df['HLA-C*04:01:01:11'].apply(lambda x: int(x.split('/')[0]))

# Compare the columns and get the rows where 'HLA-C*04:01:01:01_value' is larger than 'HLA-C*04:01:01:11_value' over a cutoff
cutoff = 20
# result = df[df['HLA-C*04:01:01:01_value'] + cutoff < df['HLA-C*04:01:01:11_value'] ]

for index,row in df.iterrows():
    if row['HLA-C*04:01:01:01_value'] + cutoff < row['HLA-C*04:01:01:11_value']:
        print(row['allele'], row['HLA-C*04:01:01:01'], row['HLA-C*04:01:01:11'], row['HLA-C*04:01:01:01_value'],row['HLA-C*04:01:01:11_value'] )
