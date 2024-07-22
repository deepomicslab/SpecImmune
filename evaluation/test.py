def calculate_total_read_length(cigar_string):
    """
    Calculate the total length of a read based on a CIGAR string,
    skipping deletion lengths.

    Parameters:
    - cigar_string: A string representing the CIGAR operations.

    Returns:
    - Total length of the read.
    """
    total_length = 0
    current_number = ''
    count = 0  # To count valid operations

    for char in cigar_string:
        if char.isdigit():
            current_number += char  # Build the number
        else:
            if current_number:
                length = int(current_number)
                if char in ["S", "M", "I"]:  # Count for soft clip, match, and insertion
                    total_length += length
                    count += 1  # Increment count for valid operations
                elif char == "D":
                    # Skip deletion; do not increment total_length or count
                    pass
                else:
                    raise ValueError(f"Invalid CIGAR operation: {char}")
                
                current_number = ''  # Reset for the next number

    print("Count of valid operations:", count)
    return total_length

# # Example usage
# cigar_string = "70M10I5D15S"
# total_length = calculate_total_read_length(cigar_string)
# print("Total read length:", total_length)

# Example usage
cigar_string = "23481S340M1D100M3D60M1I25M2D106M1D62M3D8M1D13M2I2M1D148M1D2M1D225M1D85M2I55M2D204M1D26M2D3M2D201M2I35M1D121M1I10M1D9M1I39M1D26M1I4M1I38M1D8M6D2M1D31M1I407M2D13M1D65M1D260M"
total_length = calculate_total_read_length(cigar_string)
print("Total read length:", total_length)