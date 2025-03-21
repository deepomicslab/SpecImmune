#!/usr/bin/env python3

import argparse
import pandas as pd
import os

def read_hla_file(filename):
    """Read the HLA typing file and extract key information."""
    data = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('Locus'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 7:  # Ensure there are enough columns
                locus = parts[0]
                chromosome = parts[1]
                one_guess = parts[6]
                data.append([locus, chromosome, one_guess])
    
    # Create a DataFrame
    df = pd.DataFrame(data, columns=['Locus', 'Chromosome', 'Allele'])
    
    # Pivot data to wide format, with two alleles per locus (for each chromosome)
    pivot_df = df.pivot(index='Locus', columns='Chromosome', values='Allele')
    pivot_df.columns = ['Allele_1', 'Allele_2']
    
    return pivot_df.reset_index()

def check_mendelian_inheritance(father_df, mother_df, child_df):
    """Evaluate if HLA loci follow Mendelian inheritance and phase alleles consistently.

    The logic is as follows:
      1. For each locus, determine which child's allele comes from the father and which from the mother by checking against the parent's alleles.
      2. Regardless of input order, always assign:
            - Child_Phased_1 as the allele transmitted from the father.
            - Child_Phased_2 as the allele transmitted from the mother.
      3. Record the parent's allele that was transmitted to the child.
    """
    # Initialize the result DataFrame by copying child data and adding extra columns.
    result_df = child_df.copy()
    result_df['Child_Allele_1'] = result_df['Allele_1']
    result_df['Child_Allele_2'] = result_df['Allele_2']
    
    # Keep only locus and allele information while adding more detailed columns.
    result_df = result_df[['Locus', 'Child_Allele_1', 'Child_Allele_2']]
    result_df['Father_Allele_1'] = None
    result_df['Father_Allele_2'] = None
    result_df['Mother_Allele_1'] = None
    result_df['Mother_Allele_2'] = None
    result_df['Mendelian'] = False
    result_df['Child_Phased_1'] = None  # Always stores the allele from the father.
    result_df['Child_Phased_2'] = None  # Always stores the allele from the mother.
    result_df['Father_Phased_1'] = None
    result_df['Father_Phased_2'] = None
    result_df['Mother_Phased_1'] = None
    result_df['Mother_Phased_2'] = None
    
    for i, row in result_df.iterrows():
        locus = row['Locus']
        
        # Retrieve the parent's alleles for the locus.
        father_row = father_df[father_df['Locus'] == locus]
        mother_row = mother_df[mother_df['Locus'] == locus]
        
        if len(father_row) > 0 and len(mother_row) > 0:
            father_alleles = [father_row['Allele_1'].values[0], father_row['Allele_2'].values[0]]
            mother_alleles = [mother_row['Allele_1'].values[0], mother_row['Allele_2'].values[0]]
            child_alleles = [row['Child_Allele_1'], row['Child_Allele_2']]
            
            # Save the parent's original allele information.
            result_df.at[i, 'Father_Allele_1'] = father_alleles[0]
            result_df.at[i, 'Father_Allele_2'] = father_alleles[1]
            result_df.at[i, 'Mother_Allele_1'] = mother_alleles[0]
            result_df.at[i, 'Mother_Allele_2'] = mother_alleles[1]
            
            # Check Mendelian inheritance:
            # Each child's allele should come from one of the parents.
            mendelian_ok = False
            child_from_father = None
            child_from_mother = None

            if (child_alleles[0] in father_alleles) and (child_alleles[1] in mother_alleles):
                mendelian_ok = True
                child_from_father = child_alleles[0]
                child_from_mother = child_alleles[1]
            elif (child_alleles[0] in mother_alleles) and (child_alleles[1] in father_alleles):
                mendelian_ok = True
                # Force assignment: first allele from father, second allele from mother.
                child_from_father = child_alleles[1]
                child_from_mother = child_alleles[0]
            
            result_df.at[i, 'Mendelian'] = mendelian_ok
            
            if mendelian_ok:
                # Always assign the child's phased allele:
                result_df.at[i, 'Child_Phased_1'] = child_from_father + " (from Father)"
                result_df.at[i, 'Child_Phased_2'] = child_from_mother + " (from Mother)"
                
                # For the father, determine the allele not transmitted.
                if father_alleles[0] == child_from_father:
                    father_nontransmitted = father_alleles[1]
                else:
                    father_nontransmitted = father_alleles[0]
                result_df.at[i, 'Father_Phased_1'] = child_from_father + " (to Child)"
                result_df.at[i, 'Father_Phased_2'] = father_nontransmitted
                
                # For the mother, determine the allele not transmitted.
                if mother_alleles[0] == child_from_mother:
                    mother_nontransmitted = mother_alleles[1]
                else:
                    mother_nontransmitted = mother_alleles[0]
                result_df.at[i, 'Mother_Phased_1'] = child_from_mother + " (to Child)"
                result_df.at[i, 'Mother_Phased_2'] = mother_nontransmitted
    
    return result_df

def consistent_phasing(results):
    """Check if phasing is consistent across all loci and perform genome-wide phasing.

    Since the previous function enforces a consistent order (Child_Phased_1 from Father, 
    Child_Phased_2 from Mother), this step mainly verifies the consistency. In the rare 
    event inconsistency is detected, the phased columns are swapped.
    """
    if not any(results['Mendelian']):
        return results
    
    mendelian_results = results[results['Mendelian']]
    if len(mendelian_results) <= 1:
        return results  # Only one Mendelian-locus; no need for consistency check.
    
    # Use the first Mendelian locus as reference and verify that Child_Phased_1 contains "from Father".
    first_mendelian = mendelian_results.iloc[0]
    ref_is_father = 'from Father' in first_mendelian['Child_Phased_1']
    
    consistent = True
    for _, row in mendelian_results.iloc[1:].iterrows():
        current_is_father = 'from Father' in row['Child_Phased_1']
        if current_is_father != ref_is_father:
            consistent = False
            break
    
    if not consistent:
        # If any discrepancies are found, swap the phasing.
        for i, row in results.iterrows():
            if row['Mendelian']:
                results.at[i, 'Child_Phased_1'], results.at[i, 'Child_Phased_2'] = \
                    results.at[i, 'Child_Phased_2'], results.at[i, 'Child_Phased_1']
                
                results.at[i, 'Father_Phased_1'], results.at[i, 'Father_Phased_2'] = \
                    results.at[i, 'Father_Phased_2'], results.at[i, 'Father_Phased_1']
                
                results.at[i, 'Mother_Phased_1'], results.at[i, 'Mother_Phased_2'] = \
                    results.at[i, 'Mother_Phased_2'], results.at[i, 'Mother_Phased_1']
    
    return results

def main():
    parser = argparse.ArgumentParser(description='Phase HLA genotypes for a trio family')
    parser.add_argument("-p1", '--father', required=True, help='Father HLA genotype file')
    parser.add_argument("-p2", '--mother', required=True, help='Mother HLA genotype file')
    parser.add_argument("-c", '--child', required=True, help='Child HLA genotype file')
    parser.add_argument("-o", '--output', required=True, help='Output file name')
    
    args = parser.parse_args()
    
    # Read the input files.
    father_df = read_hla_file(args.father)
    mother_df = read_hla_file(args.mother)
    child_df = read_hla_file(args.child)
    
    # Check for Mendelian inheritance and perform the initial phasing.
    result_df = check_mendelian_inheritance(father_df, mother_df, child_df)
    
    # Perform a global consistency check for the phased data.
    final_results = consistent_phasing(result_df)
    
    # Save the results.
    final_results.to_csv(args.output, index=False, sep='\t')
    print(f"Results saved to {args.output}")
    
    # Print a summary.
    mendelian_count = sum(final_results['Mendelian'])
    total_count = len(final_results)
    print(f"\nSummary:")
    print(f"Total HLA loci: {total_count}")
    print(f"Loci following Mendelian inheritance: {mendelian_count} ({mendelian_count/total_count*100:.1f}%)")
    print(f"Loci with successful phasing: {mendelian_count}")
    
    # Report the loci that do not follow Mendelian inheritance.
    if mendelian_count < total_count:
        non_mendelian = final_results[~final_results['Mendelian']]
        print("\nLoci not following Mendelian inheritance:")
        for _, row in non_mendelian.iterrows():
            print(f"  {row['Locus']}: Child ({row['Child_Allele_1']}, {row['Child_Allele_2']}), "
                  f"Father ({row['Father_Allele_1']}, {row['Father_Allele_2']}), "
                  f"Mother ({row['Mother_Allele_1']}, {row['Mother_Allele_2']})")

if __name__ == "__main__":
    main()