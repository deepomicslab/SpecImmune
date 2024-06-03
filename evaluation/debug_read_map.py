import pandas as pd


# # Split the values and convert to integers
# df['HLA-C*04:01:01:01_value'] = df['HLA-C*04:01:01:01'].apply(lambda x: int(x.split('/')[0]))
# df['HLA-C*04:01:01:11_value'] = df['HLA-C*04:01:01:11'].apply(lambda x: int(x.split('/')[0]))

# # Compare the columns and get the rows where 'HLA-C*04:01:01:01_value' is larger than 'HLA-C*04:01:01:11_value' over a cutoff
# cutoff = 20
# # result = df[df['HLA-C*04:01:01:01_value'] + cutoff < df['HLA-C*04:01:01:11_value'] ]

# for index,row in df.iterrows():
#     if row['HLA-C*04:01:01:01_value'] + cutoff < row['HLA-C*04:01:01:11_value']:
#         print(row['allele'], row['HLA-C*04:01:01:01'], row['HLA-C*04:01:01:11'], row['HLA-C*04:01:01:01_value'],row['HLA-C*04:01:01:11_value'] )


def compare_alleles(df, allele1, allele2, cutoff=0):
    df[f'{allele1}_value'] = df[allele1].apply(lambda x: int(x.split('/')[0]))
    df[f'{allele2}_value'] = df[allele2].apply(lambda x: int(x.split('/')[0]))

    df[f'{allele1}_ide'] = df[allele1].apply(lambda x: float(x.split('/')[1]))
    df[f'{allele2}_ide'] = df[allele2].apply(lambda x: float(x.split('/')[1]))

    for index, row in df.iterrows():
        # if row[f'{allele1}_value'] + cutoff < row[f'{allele2}_value']:
        if row[f'{allele1}_ide'] + cutoff < row[f'{allele2}_ide']:
            print(row['allele'], row[allele1], row[allele2], row[f'{allele1}_value'], row[f'{allele2}_value'])

# Read the CSV file
df = pd.read_csv('/mnt/d/HLAPro_backup/Nanopore_optimize/output5/fredhutch-hla-KOSE/fredhutch-hla-KOSE.HLA-A.read.matrix.csv', delimiter=',')


# Usage
compare_alleles(df, 'HLA-A*02:01:01:01', 'HLA-A*02:01:01:110')