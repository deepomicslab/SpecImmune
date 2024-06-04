import scipy
from scipy.stats import binom


def get_binom_pmf(k1, k2):
    p = 0.5
    n = k1 + k2
    b = scipy.stats.binom(n, p)


    # print (binom.logpmf(k1, n, p, loc=0))
    return b.pmf(k1)
# print (binom.logpmf(k, n, p, loc=0))

def if_homo(record_allele_pair_sep_match, first_pair):
    allele_name_list = first_pair.split("&")
    # print (record_allele_pair_sep_match[first_pair])
    dp1 = round(record_allele_pair_sep_match[first_pair][allele_name_list[0]]["depth"])
    dp2 = round(record_allele_pair_sep_match[first_pair][allele_name_list[1]]["depth"])
    p_value = get_binom_pmf(dp1, dp2)
    print (dp1, dp2, p_value)
    # if p_value < 0.01:
    #     return True
    # else:
    #     return False
    return p_value


if __name__ == "__main__":

    k1 = 742
    k2 = 738



    print (get_binom_pmf(k1, k2))


