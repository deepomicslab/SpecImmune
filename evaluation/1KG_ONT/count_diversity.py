import numpy as np

diversity = "hla/diversity_df.csv"


div_dict = {}

for line in open(diversity):
    line = line.strip().split(",")
    if line[0] == "Population":
        continue
    if line[1] not in div_dict:
        div_dict[line[1]] = []
    div_dict[line[1]].append(float(line[2]))

for key in div_dict:
    div_dict[key] = np.mean(div_dict[key])

## sort by value
sorted_div = sorted(div_dict.items(), key=lambda x: x[1])
for item in sorted_div:
    print(item[0], item[1])
    
