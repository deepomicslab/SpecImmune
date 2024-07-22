def calculate_total_read_length(cigar_string):
    """
    Calculate the total length of a read based on a CIGAR string,
    skipping deletion and hard clipping lengths.

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
                elif char in ["D", "H"]:
                    # Skip deletion and hard clipping; do not increment total_length or count
                    pass
                else:
                    raise ValueError(f"Invalid CIGAR operation: {char}")
                
                current_number = ''  # Reset for the next number

    print("Count of valid operations:", count)
    return total_length

# Example usage
cigar_string = "9827H319M3I368M1D23M1D34M1I50M1I190M1I19M1I13M2I7M1I288M3D158M1D4M1D4M1D13M1I224M1D28M1I25M1D67M1D9M1D3M1D112M2D11M1D24M1I167M1I155M1I112M1I6M6D8M1D93M1D13M1D12M2I162M2D22M"
total_length = calculate_total_read_length(cigar_string)
print("Total read length:", total_length)