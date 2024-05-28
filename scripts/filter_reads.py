from collections import defaultdict
import numpy as np


def remove_items(dictionary, keys_to_remove):
    for key in keys_to_remove:
        dictionary.pop(key, None)

def get_quantile_number(data, quantile):
    sorted_data = np.sort(data)
    total_elements = len(sorted_data)
    index = (quantile / 100) * (total_elements - 1)
    floor_index = int(index)
    ceil_index = floor_index + 1
    quantile_number = (sorted_data[floor_index] + sorted_data[ceil_index]) / 2
    return quantile_number

def examine_reads(record_read_allele_dict, min_ave_len = 500):
    identity_dict = defaultdict()
    match_dict = defaultdict()

    all_read_ave_match = []

    for read_name in record_read_allele_dict:
        identity_list = []
        match_list = []
        for allele_name in record_read_allele_dict[read_name]:
            identity_list.append(record_read_allele_dict[read_name][allele_name].identity)
            match_list.append(record_read_allele_dict[read_name][allele_name].match_num)
    

        sort_match_list = sorted(match_list)
        # print (read_name, np.mean(match_list), np.median(match_list), sort_match_list[0], sort_match_list[-1])
        median_match = np.median(match_list)
        all_read_ave_match.append(median_match)
        match_dict[read_name] = median_match

    quantile_number = get_quantile_number(all_read_ave_match, 10)
    print ("quantile_number", quantile_number)

    remove_reads = []
    for read_name in record_read_allele_dict:
        if match_dict[read_name] < quantile_number:
            remove_reads.append(read_name)
    print (len(record_read_allele_dict))
    remove_items(record_read_allele_dict, remove_reads)
    print (len(record_read_allele_dict))

    return record_read_allele_dict


    