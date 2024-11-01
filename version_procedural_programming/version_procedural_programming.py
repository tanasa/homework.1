#!/usr/bin/env python
# coding: utf-8


# The steps to solve the problem are as follows:

# 1. First, we convert the CIGAR string into an expanded format, such as:
#    MMMMMMMM-------MMMMMM++MM-----------MMMMMMM
#    MMMMMMMMMMMMMMMMMMMM
#    This transformation allows us to track the position of each nucleotide within the expanded string representation.

# 2. Next, we align the positions in the expanded string representation to the genome sequence.

# 3. Finally, we print the required information: transcript name, nucleotide index, chromosome name,
#    and the genome position corresponding to each nucleotide index.


import os
import re  
import pandas as pd


def transform_cigar(cigar_string):

    '''
    Transforms a CIGAR string into an expanded representation of its alignment operations.

    The function converts the input CIGAR string, which specifies alignment operations
    (e.g., "8M7D6M2I2M11D7M", or "10M"), into a string where each operation is 
    represented by repeated characters. Specifically :
        - 'M' (match): Represented by 'M' repeated the specified number of times
        - 'D' (deletion): Represented by '-' repeated the specified number of times
        - 'I' (insertion): Represented by '+' repeated the specified number of times

    Parameters:
    -----------
    cigar_string : str
        The CIGAR string to be transformed. It consists of one or more 
        sections of an integer followed by a character (M, I, or D).

    Returns ::
    --------
    str
        A single string where each alignment operation in the CIGAR string 
        is expanded into its detailed character representation.

    Example:
    --------
    >>> transform_cigar("8M7D6M2I2M11D7M")
    'MMMMMMMM-------MMMMMM++MM-----------MMMMMMM'

    >>> transform_cigar("10M")
    'MMMMMMMMMM'

    '''
    
    # Initialize an empty list to store the transformed string
    transformed_string = []

    # Use a regular expression to find all matches in the CIGAR string
    matches = re.finditer(r'(\d+)([MID])', cigar_string)

    for match in matches:
        length = int(match.group(1))                 # Get the length
        operation = match.group(2)                   # Get the operation characters
        
        if operation == 'M':                         
            transformed_string.append('M' * length)  # Append 'M' length times
        elif operation == 'D':                       
            transformed_string.append('-' * length)  # Append '-' for deletions 
        elif operation == 'I':                       
            transformed_string.append('+' * length)  # Append '+' for insertions 

    # Join the list of characters into a single string and return
    return ''.join(transformed_string)


# print("Testing the function : transform_cigar():") 
# print("more information is available in the xls file.") 
# cigar_string1 = "8M7D6M2I2M11D7M"
# transformed_string1 = transform_cigar(cigar_string1)
# print(transformed_string1)  



def find_transcript_index(transformed_string, transcript_nucleotide_position):

    '''
    Finds the index of a specific NUCLEOTIDE NUMBER in the expanded CIGAR representation.

    The function takes a transformed CIGAR string (where 'M' indicates a match,
    '-' indicates a deletion, and '+' indicates an insertion) and a nucleotide number
    in the transcript, and returns the corresponding index in the transformed string.
    
    The function adjusts for 0-based indexing by subtracting 1 from the provided
    `transcript_nucleotide_position` if it is greater than 0. The function processes
    each character in the transformed string and keeps track of both the original transcript
    position and the current index in the transformed string.

    Parameters:
    -----------
    transformed_string : str
        The expanded representation of the CIGAR string, consisting of characters 'M',
        '-', and '+'.

    transcript_nucleotide_position : int
        The 1-based position of the nucleotide in the original transcript. 

    Returns:
    --------
    int
        The index of the nucleotide in the transformed string if found; otherwise,
        returns -1 if the nucleotide does not map to a position in the transformed string.

    Example:
    --------
    >>> transformed_string = transform_cigar("8M7D6M2I2M11D7M")
    >>> find_transcript_index(transformed_string, 19)
    36  # Returns the index in the transformed string for the given nucleotide position

    Notes:
    ------
    If the provided nucleotide position corresponds to an insertion, the function
    will indicate that mapping to an insertion is likely and will return -1.

    '''

    # Adjust for 0-based indexing
    if transcript_nucleotide_position > 0: 
                  transcript_nucleotide_position -= 1
    
    index_in_original_string = 0                      # Track the current transcript position
    index_in_transformed_string = 0                   # Track the index in the transformed string
    
    for char in transformed_string:
        if char == 'M' :                              # For matches or pluses (insertions), we advance both indexes
            if index_in_original_string == transcript_nucleotide_position:
               return index_in_transformed_string 
            index_in_original_string += 1
            index_in_transformed_string += 1 
        elif char == '-':                             # For deletions, only advance the index in the transformed string
            index_in_transformed_string += 1
        elif char == '+' :                            # For matches or pluses (insertions), we advance both indexes
            index_in_original_string += 1
            index_in_transformed_string += 1 
            # pass     
        
    print(f"Most likely, it maps to an insertion")
    return -1


# print("Testing the function find_transcript_index():")  
# print("more information is available in the xls file.") 
# transformed_string = transform_cigar("8M7D6M2I2M11D7M")
# print(transformed_string)
# nucleotide_number = 19  
# index = find_transcript_index(transformed_string, nucleotide_number)
# print(f'Nucleotide number {nucleotide_number} corresponds to index {index} in the expanded CIGAR string.')


def calculate_genome_position(transcript_modified_string, transcript_position_modified_string, position_start_transcript_on_genome):

    '''
    Calculates the genome position corresponding to a given nucleotide position in a modified transcript string.

    This function processes a modified transcript string, which consists of characters representing matches ('M'),
    deletions ('-'), and insertions ('+'). It maps a specified position in the modified transcript to its 
    corresponding position in the genome, taking into account the starting position of the transcript on the genome.

    Parameters:
    -----------
    transcript_modified_string : str
        The modified transcript string composed of characters 'M', '-', and '+' that represent matches, 
        deletions, and insertions, respectively.

    transcript_position_modified_string : int
        The 0-based position in the modified transcript for which the corresponding genome position is to be found.

    position_start_transcript_on_genome : int
        The starting position of the transcript in the genome, serving as a reference point for position calculations.

    Returns:
    --------
    int or None
        The genome position corresponding to the specified transcript position if found; otherwise, returns None
        to indicate that the position was not found in the modified transcript string.

    Example:
    --------
    >>> position_start_transcript_on_genome = 3
    >>> transcript_string_modified_string = "MMMMMMMM-------MMMMMM++MM-----------MMMMMMM"
    >>> transcript_position_modified_string = 42
    >>> calculate_genome_position(transcript_string_modified_string, transcript_position_modified_string, position_start_transcript_on_genome)
    45  # Returns the corresponding genome index

    Notes:
    ------
    - The function will return the genome position when the specified transcript position is encountered
      during the traversal of the modified transcript string.
    - If the specified position is not found, the function returns None.
    - The function assumes that the provided `transcript_position_modified_string` is within the bounds of the 
      `transcript_modified_string`.
    '''
    
    current_genome_position = position_start_transcript_on_genome
    current_position_on_alignment = 0                   # Track the current position in the alignment (expanded CIGAR string) 

    # Loop through the expanded transcript string
    for char in transcript_modified_string:
        if char == 'M':
            # For matches, advance both expnded string and genome indexes
            if current_position_on_alignment == transcript_position_modified_string:
                return current_genome_position 
            current_position_on_alignment += 1    # Advance the alignment position
            current_genome_position += 1          # Advance the genome position
        elif char == '-':
            # For deletions, skip the expanded string position and advance the genome position
            if current_position_on_alignment == transcript_position_modified_string:
                return current_genome_position
            current_position_on_alignment += 1    # Advance the alignment position
            current_genome_position += 1          # Advance the genome position
        elif char == '+':
            current_position_on_alignment += 1    # Advance the alignment position

    # If the transcript_position is not found, return None or 
    return None         # some appropriate value to indicate not found



# print("Testing the function find_transcript_index():") 
# print("more information is available in the xls file.")
# position_start_transcript_on_genome = 3 # for TR1
# transcript_string_modified_string = "MMMMMMMM-------MMMMMM++MM-----------MMMMMMM"
# transcript_position_modified_string = 19  
# genome_position = calculate_genome_position(transcript_string_modified_string, transcript_position_modified_string, position_start_transcript_on_genome)
# print(f'Transcript position in the expanded CIGAR string {transcript_position_modified_string} corresponds to genome index {genome_position}.')

# Verifying the entire pipeline on TR1 :

print("Verifying the entire pipeline on TR1 :")
position_start_transcript_on_genome = 3  # The index of the genome sequence that coincides with the transcript_start

transformed_string = transform_cigar("8M7D6M2I2M11D7M")
print(transformed_string)
nucleotide_number = 19                         # the NUCLEOTIDE NUMBER in the original transcript
transcript_position = nucleotide_number        # to make it equivalent with the INDEX

index = find_transcript_index(transformed_string, transcript_position)
print(f'Nucleotide number {nucleotide_number} corresponds to index {index} on the expanded CIGAR string.')
transcript_string_modified_string = transformed_string
transcript_position_modified_string = index    

genome_position = calculate_genome_position(transcript_string_modified_string, index, position_start_transcript_on_genome)
print(f'Nucleotide number {transcript_position} corresponds to the genome index {genome_position}.')


# Verifying the entire pipeline on TR2 :

print("Verifying the entire pipeline on TR2 :")
position_start_transcript_on_genome = 10 # The index of the genome sequence that coincides with the transcript_start
transformed_string = transform_cigar("20M")
print(transformed_string)

transcript_position = 11  # The nucleotide position in the original transcript that is mapped on expanded CIGAR string
index = find_transcript_index(transformed_string, transcript_position)
print(f'Nucleotide number {nucleotide_number} corresponds to index {index} on the expanded CIGAR string.')
transcript_string_modified_string = transformed_string
transcript_position_modified_string = index 

genome_position = calculate_genome_position(transcript_string_modified_string, index, position_start_transcript_on_genome)
print(f'Nucleotide number {transcript_position} corresponds to the genome index {genome_position}.')



# Read the file that contains the transcript alignments

def read_transcript_alignments(filename):

    df = pd.read_csv(filename, sep='\t', header=0, names=['transcript_name', 'chromosome', 'transcript_start', 'cigar_string'])
    df['transcript_start'] = df['transcript_start'].astype(int)
    df = df.rename(columns={'transcript_start': 'start_transcript_on_genome'})

    return df

filename1 = 'input_transcripts_alignments.txt'  
transcript_alignments = read_transcript_alignments(filename1)
# print(transcript_alignments)


# Read the file that contains the nucleotide positions 

def read_nucleotide_positions(filename):
    
    df = pd.read_csv(filename, sep='\t', header=0, names=['transcript_name','transcript_coordinate'])

    df['transcript_coordinate'] = df['transcript_coordinate'].astype(int)
    df = df.rename(columns={'transcript_coordinate': 'index_in_transcript'})

    return df


filename2 = 'input_transcripts_positions.txt'  
nucleotide_positions = read_nucleotide_positions(filename2)
# print(nucleotide_positions)


# Merge these DataFrames based on 'transcript_name'

combined_data = pd.merge(transcript_alignments, nucleotide_positions, on='transcript_name', how='outer')
print("\nThe combined files that contain both transcript alignments and nucleotide positions :")
print(combined_data)
# output_file = 'input_transcripts_alignments_and_nucleotide_positions.txt' 
# combined_data.to_csv(output_file, sep='\t', index=False) 


# Applying the 1st function : transform_cigar()


combined_data["transformed_string"] = combined_data["cigar_string"].apply(transform_cigar) 
# print(combined_data)


# Applying the 2nd function : find_transcript_index()
# Apply the second function with additional argument using a lambda function


combined_data["index_in_transformed_string"] = combined_data.apply(
    lambda row: find_transcript_index(row["transformed_string"], row["index_in_transcript"] + 1),
    axis=1
)
# print(combined_data)


# Applying the 3rd function : calculate_genome_position()


combined_data["genome_index"] = combined_data.apply(
    lambda row: calculate_genome_position(
        row["transformed_string"],
        row["index_in_transformed_string"]  ,
        row["start_transcript_on_genome"]
    ),
    axis=1
)
print(combined_data)


# Printing the columns of interest on the screen and in a file :

print("Results :")
combined_data = combined_data.rename(columns={'index_in_transcript': 'transcript_index'})
print(combined_data.loc[:, ["transcript_name", "transcript_index", "chromosome", "genome_index"]])

combined_data[["transcript_name", "transcript_index", "chromosome", "genome_index"]].to_csv("output_genome_positions.txt", sep='\t', index=False)