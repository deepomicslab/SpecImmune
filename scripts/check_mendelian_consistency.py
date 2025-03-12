import argparse

def parse_genotype(genotype_str):
    """Parse genotype string and return a set of alleles"""
    if not genotype_str or genotype_str == "NA":
        return set()
    return set(genotype_str.split(";"))

def read_hla_file(file_path):
    """Read HLA typing file and organize data by locus"""
    hla_data = {}  # Structure: {locus: {alleles for all chromosomes}}
    
    with open(file_path, 'r', encoding='utf-8') as file:
        for line in file:
            line = line.strip()
            if not line or line.startswith("#"):
                continue  # Skip empty lines and comments
                
            fields = line.split("\t")
            if len(fields) < 3:
                continue  # Ensure data integrity
                
            locus, chromosome, genotype = fields[0], fields[1], fields[2]
            
            # Initialize data structure for this locus
            if locus not in hla_data:
                hla_data[locus] = {"alleles": set(), "is_na": False}
                
            # Add alleles from this chromosome to the locus set
            if genotype == "NA":
                hla_data[locus]["is_na"] = True
            else:
                alleles = parse_genotype(genotype)
                hla_data[locus]["alleles"].update(alleles)
            
    return hla_data

def check_mendelian_consistency(parent1_file, parent2_file, child_file, output_file):
    """Check Mendelian consistency and output results"""
    # Read data files
    parent1_data = read_hla_file(parent1_file)
    parent2_data = read_hla_file(parent2_file)
    child_data = read_hla_file(child_file)
    
    results = []
    
    # Debug info
    print(f"Father loci count: {len(parent1_data)}")
    print(f"Mother loci count: {len(parent2_data)}")
    print(f"Child loci count: {len(child_data)}")
    
    # Check each locus present in the child's data
    for locus in child_data:
        # Skip loci with missing parental data
        if locus not in parent1_data or locus not in parent2_data:
            print(f"Skipping locus {locus}: incomplete parental data")
            continue
            
        child_alleles = child_data[locus]["alleles"]
        parent1_alleles = parent1_data[locus]["alleles"]
        parent2_alleles = parent2_data[locus]["alleles"]
        
        # If either parent has NA genotype, mark as unknown
        if parent1_data[locus]["is_na"] or parent2_data[locus]["is_na"]:
            results.append([locus, "unknown", 
                            "NA", 
                            "NA", 
                            ";".join(child_alleles) if child_alleles else "NA", 
                            "NA"])
            continue
            
        # Debug info
        print(f"Locus {locus}:")
        print(f"  Father alleles: {parent1_alleles}")
        print(f"  Mother alleles: {parent2_alleles}")
        print(f"  Child alleles: {child_alleles}")
        
        # If child has no alleles (NA), skip detailed checking
        if not child_alleles:
            results.append([locus, "unknown", 
                            ";".join(parent1_alleles) if parent1_alleles else "NA", 
                            ";".join(parent2_alleles) if parent2_alleles else "NA", 
                            "NA", 
                            "NA"])
            continue
        
        # All alleles from both parents
        all_parent_alleles = parent1_alleles.union(parent2_alleles)
        
        # Check if child's alleles are a subset of parents' alleles
        if child_alleles.issubset(all_parent_alleles):
            # Determine whether there's a possible combination consistent with Mendelian inheritance
            consistent = False
            allele_origin = "Cannot determine specific origin"
            
            # If child has exactly two alleles, try to determine origin
            if len(child_alleles) == 2:
                child_list = list(child_alleles)
                a1, a2 = child_list[0], child_list[1]
                
                # Try to determine the origin of each allele
                if (a1 in parent1_alleles and a2 in parent2_alleles):
                    consistent = True
                    allele_origin = f"{a1}(father),{a2}(mother)"
                elif (a1 in parent2_alleles and a2 in parent1_alleles):
                    consistent = True
                    allele_origin = f"{a1}(mother),{a2}(father)"
                # Parents may have shared alleles, origin cannot be precisely determined
                elif a1 in parent1_alleles and a1 in parent2_alleles:
                    consistent = True
                elif a2 in parent1_alleles and a2 in parent2_alleles:
                    consistent = True
            
            # More lenient check: if each allele comes from either parent, consider it consistent
            else:
                consistent = True
            
            if consistent:
                results.append([locus, "consistent", 
                                ";".join(parent1_alleles), 
                                ";".join(parent2_alleles), 
                                ";".join(child_alleles), 
                                allele_origin])
            else:
                results.append([locus, "inconsistent (abnormal combination)", 
                                ";".join(parent1_alleles), 
                                ";".join(parent2_alleles), 
                                ";".join(child_alleles), 
                                "Cannot form valid combination"])
        else:
            # Child has alleles not present in either parent
            new_alleles = child_alleles - all_parent_alleles
            results.append([locus, "inconsistent", 
                            ";".join(parent1_alleles), 
                            ";".join(parent2_alleles), 
                            ";".join(child_alleles), 
                            ";".join(new_alleles)])
    
    # Output results
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write("Locus\tMendelian_Consistency\tParent1_Alleles\tParent2_Alleles\tChild_Alleles\tAllele_Origin\n")
        for row in results:
            f.write("\t".join(row) + "\n")
    
    # Print warning if no results
    if not results:
        print("Warning: No results generated. Please check input file format and content.")
    else:
        print(f"Generated {len(results)} result records")
        
    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Check HLA Mendelian consistency")
    parser.add_argument("-p1", "--parent1", required=True, help="Father's HLA typing file")
    parser.add_argument("-p2", "--parent2", required=True, help="Mother's HLA typing file")
    parser.add_argument("-c", "--child", required=True, help="Child's HLA typing file")
    parser.add_argument("-o", "--output", required=True, help="Output result file")

    args = parser.parse_args()

    check_mendelian_consistency(args.parent1, args.parent2, args.child, args.output)