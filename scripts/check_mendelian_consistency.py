#!/usr/bin/env python3

import argparse
import pandas as pd
import os

def read_hla_file(filename):
    """Read HLA typing file and extract key information"""
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
    
    # Create DataFrame
    df = pd.DataFrame(data, columns=['Locus', 'Chromosome', 'Allele'])
    
    # Pivot data to wide format, with two alleles per locus (chromosome 1 and 2)
    pivot_df = df.pivot(index='Locus', columns='Chromosome', values='Allele')
    pivot_df.columns = ['Allele_1', 'Allele_2']
    
    return pivot_df.reset_index()

def check_mendelian_inheritance(father_df, mother_df, child_df):
    """Check if HLA loci follow Mendelian inheritance"""
    result_df = child_df.copy()
    result_df['Child_Allele_1'] = result_df['Allele_1']
    result_df['Child_Allele_2'] = result_df['Allele_2']
    
    # Remove the original columns after copying their values
    result_df = result_df[['Locus', 'Child_Allele_1', 'Child_Allele_2']]
    result_df['Father_Allele_1'] = None
    result_df['Father_Allele_2'] = None
    result_df['Mother_Allele_1'] = None
    result_df['Mother_Allele_2'] = None
    result_df['Mendelian'] = False
    result_df['Child_Phased_1'] = None
    result_df['Child_Phased_2'] = None
    result_df['Father_Phased_1'] = None
    result_df['Father_Phased_2'] = None
    result_df['Mother_Phased_1'] = None
    result_df['Mother_Phased_2'] = None
    
    for i, row in result_df.iterrows():
        locus = row['Locus']
        
        # Get parents' alleles at this locus
        father_row = father_df[father_df['Locus'] == locus]
        mother_row = mother_df[mother_df['Locus'] == locus]
        
        if len(father_row) > 0 and len(mother_row) > 0:
            father_alleles = [father_row['Allele_1'].values[0], father_row['Allele_2'].values[0]]
            mother_alleles = [mother_row['Allele_1'].values[0], mother_row['Allele_2'].values[0]]
            child_alleles = [row['Child_Allele_1'], row['Child_Allele_2']]
            
            result_df.at[i, 'Father_Allele_1'] = father_alleles[0]
            result_df.at[i, 'Father_Allele_2'] = father_alleles[1]
            result_df.at[i, 'Mother_Allele_1'] = mother_alleles[0]
            result_df.at[i, 'Mother_Allele_2'] = mother_alleles[1]
            
            # Check Mendelian inheritance
            # Each child allele must come from either father or mother
            child1_from_father = child_alleles[0] in father_alleles
            child1_from_mother = child_alleles[0] in mother_alleles
            child2_from_father = child_alleles[1] in father_alleles
            child2_from_mother = child_alleles[1] in mother_alleles
            
            is_mendelian = ((child1_from_father and child2_from_mother) or 
                            (child1_from_mother and child2_from_father))
            
            result_df.at[i, 'Mendelian'] = is_mendelian
            
            # Attempt phasing
            if is_mendelian:
                if child1_from_father and child2_from_mother:
                    # Find the exact parental alleles
                    father_match = father_alleles[0] if father_alleles[0] == child_alleles[0] else father_alleles[1]
                    mother_match = mother_alleles[0] if mother_alleles[0] == child_alleles[1] else mother_alleles[1]
                    
                    result_df.at[i, 'Child_Phased_1'] = child_alleles[0] + " (from Father)"
                    result_df.at[i, 'Child_Phased_2'] = child_alleles[1] + " (from Mother)"
                    
                    # Mark parents' phased haplotypes
                    result_df.at[i, 'Father_Phased_1'] = father_match + " (to Child)"
                    result_df.at[i, 'Father_Phased_2'] = father_alleles[0] if father_alleles[1] == father_match else father_alleles[1]
                    
                    result_df.at[i, 'Mother_Phased_1'] = mother_match + " (to Child)"
                    result_df.at[i, 'Mother_Phased_2'] = mother_alleles[0] if mother_alleles[1] == mother_match else mother_alleles[1]
                    
                elif child1_from_mother and child2_from_father:
                    # Find the exact parental alleles
                    father_match = father_alleles[0] if father_alleles[0] == child_alleles[1] else father_alleles[1]
                    mother_match = mother_alleles[0] if mother_alleles[0] == child_alleles[0] else mother_alleles[1]
                    
                    result_df.at[i, 'Child_Phased_1'] = child_alleles[0] + " (from Mother)"
                    result_df.at[i, 'Child_Phased_2'] = child_alleles[1] + " (from Father)"
                    
                    # Mark parents' phased haplotypes
                    result_df.at[i, 'Father_Phased_1'] = father_match + " (to Child)"
                    result_df.at[i, 'Father_Phased_2'] = father_alleles[0] if father_alleles[1] == father_match else father_alleles[1]
                    
                    result_df.at[i, 'Mother_Phased_1'] = mother_match + " (to Child)"
                    result_df.at[i, 'Mother_Phased_2'] = mother_alleles[0] if mother_alleles[1] == mother_match else mother_alleles[1]
    
    return result_df

def consistent_phasing(results):
    """Check if phasing is consistent across all loci and perform genome-wide phasing"""
    # If no loci follow Mendelian inheritance, consistent phasing is not possible
    if not any(results['Mendelian']):
        return results
    
    mendelian_results = results[results['Mendelian']]
    if len(mendelian_results) <= 1:
        return results  # Only one locus follows Mendelian inheritance, no need to check consistency
    
    # Check phasing consistency
    # Use the first Mendelian locus as reference
    first_mendelian = mendelian_results.iloc[0]
    ref_pattern = 'from Father' in first_mendelian['Child_Phased_1']
    
    # Check if other loci are consistent with the reference
    consistent = True
    for _, row in mendelian_results.iloc[1:].iterrows():
        current_pattern = 'from Father' in row['Child_Phased_1']
        if current_pattern != ref_pattern:
            consistent = False
            break
    
    # If inconsistent, adjust phasing
    if not consistent:
        for i, row in results.iterrows():
            if row['Mendelian']:
                # Swap child's phased haplotypes
                results.at[i, 'Child_Phased_1'], results.at[i, 'Child_Phased_2'] = \
                    results.at[i, 'Child_Phased_2'], results.at[i, 'Child_Phased_1']
                
                # Adjust parents' phased haplotypes
                if 'from Father' in results.at[i, 'Child_Phased_1']:
                    results.at[i, 'Child_Phased_1'] = results.at[i, 'Child_Phased_1'].replace('from Father', 'from Mother')
                    results.at[i, 'Child_Phased_2'] = results.at[i, 'Child_Phased_2'].replace('from Mother', 'from Father')
                else:
                    results.at[i, 'Child_Phased_1'] = results.at[i, 'Child_Phased_1'].replace('from Mother', 'from Father')
                    results.at[i, 'Child_Phased_2'] = results.at[i, 'Child_Phased_2'].replace('from Father', 'from Mother')
                
                # Swap parents' phased haplotypes
                results.at[i, 'Father_Phased_1'], results.at[i, 'Father_Phased_2'] = \
                    results.at[i, 'Father_Phased_2'], results.at[i, 'Father_Phased_1']
                
                results.at[i, 'Mother_Phased_1'], results.at[i, 'Mother_Phased_2'] = \
                    results.at[i, 'Mother_Phased_2'], results.at[i, 'Mother_Phased_1']
    
    return results

def main():
    parser = argparse.ArgumentParser(description='Phase HLA genotypes for a trio family')
    parser.add_argument("-p1", "--father", required=True, help="Father's HLA typing file")
    parser.add_argument("-p2", "--mother", required=True, help="Mother's HLA typing file")
    parser.add_argument("-c", "--child", required=True, help="Child's HLA typing file")
    parser.add_argument("-o", "--output", required=True, help="Output result file")
    
    args = parser.parse_args()
    
    # Read input files
    father_df = read_hla_file(args.father)
    mother_df = read_hla_file(args.mother)
    child_df = read_hla_file(args.child)
    
    # Check Mendelian inheritance and perform initial phasing
    result_df = check_mendelian_inheritance(father_df, mother_df, child_df)
    
    # Attempt consistent phasing
    final_results = consistent_phasing(result_df)
    
    # Save results
    final_results.to_csv(args.output, index=False, sep='\t')
    print(f"Results saved to {args.output}")
    
    # Output summary
    mendelian_count = sum(final_results['Mendelian'])
    total_count = len(final_results)
    print(f"\nSummary:")
    print(f"Total HLA loci: {total_count}")
    print(f"Loci following Mendelian inheritance: {mendelian_count} ({mendelian_count/total_count*100:.1f}%)")
    print(f"Loci with successful phasing: {mendelian_count}")
    
    # Output loci not following Mendelian inheritance
    if mendelian_count < total_count:
        non_mendelian = final_results[~final_results['Mendelian']]
        print("\nLoci not following Mendelian inheritance:")
        for _, row in non_mendelian.iterrows():
            print(f"  {row['Locus']}: Child ({row['Child_Allele_1']}, {row['Child_Allele_2']}), "
                  f"Father ({row['Father_Allele_1']}, {row['Father_Allele_2']}), "
                  f"Mother ({row['Mother_Allele_1']}, {row['Mother_Allele_2']})")

if __name__ == "__main__":
    main()