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

def if_homo2(record_allele_pair_sep_match, first_pair):
    allele_name_list = first_pair.split("&")
    # print (record_allele_pair_sep_match[first_pair])
    dp1 = round(record_allele_pair_sep_match[first_pair][allele_name_list[0]]["depth"])
    dp2 = round(record_allele_pair_sep_match[first_pair][allele_name_list[1]]["depth"])
    

    # Prior parameters
    alpha_prior = 1
    beta_prior = 1

    # Data
    n_trials = dp1 + dp2
    n_successes = min(dp1, dp2)

    # Posterior parameters
    alpha_post = alpha_prior + n_successes
    beta_post = beta_prior + n_trials - n_successes

    # Mean of the posterior distribution
    posterior_mean = alpha_post / (alpha_post + beta_post)
    print (dp1, dp2, posterior_mean)
    return posterior_mean


if __name__ == "__main__":

    k1 = 742
    k2 = 738



    print (get_binom_pmf(k1, k2))


