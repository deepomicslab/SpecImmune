# geneEDGE_308449_length_2161_cov_5.311016-       NZ_CP040529.1   100.000 2161    33036   2161    0       0       1       2161    27509   29669   0.0     3991
# geneEDGE_308449_length_2161_cov_5.311016-       NC_004703.1     100.000 2161    33038   2161    0       0       1       2161    22157   24317   0.0     3991
# geneEDGE_308449_length_2161_cov_5.311016-       NZ_CP072217.1   99.893  2161    30813   935     1       0       1227    2161    13369   12435   0.0     1722
# geneEDGE_308449_length_2161_cov_5.311016-       NZ_CP072217.1   77.966  2161    30813   295     59      5       680     971     13892   13601   1.79e-46        18
# "6 qaccver saccver pident qlen slen length mismatch gapopen qstart qend sstart send evalue bitscore
import sys

if __name__ == "__main__":
    blast = open(sys.argv[1], 'r')
    res = open(sys.argv[2], "w")
    record = {}
    for idx, it in enumerate(blast.readlines()):
        print(it)
        items = it.strip().split()
        print(items)
        record[idx]={"pos":min(int(items[10]),int(items[11])), "nd_name": items[0],"start": items[10], "end":items[11]}
    sorted_record = sorted(record.items(), key=lambda x: x[1]["pos"], reverse=False)
    # print(sorted_record)
    # visited = []
    for r in sorted_record:
        res.write(r[1]['nd_name']+'_s_'+r[1]['start']+'_e_'+r[1]['end']+'\n')
    blast.close()
    res.close()
