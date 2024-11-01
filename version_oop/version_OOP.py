#!/usr/bin/env python
# coding: utf-8

# To compare it with the functional version of the code,  we propose an object - oriented version of the same algorithm. 
# Those 3 functions that we have presented in the functional version are methods in a class entitled transcriptToGenome.

# The class transcriptToGenome() contains these 3 methods :

# transformCigar
# findTranscriptIndex
# calculateGenomePosition

# that operate on the input data.         

import re
import pandas as pd
from dataclasses import dataclass

@dataclass
class transcriptToGenome:

    def transformCigar(self, cigar_string):
        '''
        Transforms a CIGAR string into an expanded representation of its alignment operations.

        Parameters:
        -----------
        cigar_string : str
            The CIGAR string to be transformed.

        Returns:
        --------
        str
            A single string where each alignment operation in the CIGAR string 
            is expanded into its detailed character representation.
        '''
        transformed_string = []

        # Use a regular expression to find all matches in the CIGAR string
        matches = re.finditer(r'(\d+)([MID])', cigar_string)

        for match in matches:
            length = int(match.group(1))  # Get the length
            operation = match.group(2)     # Get the operation characters
            
            if operation == 'M':                         
                transformed_string.append('M' * length)  # Append 'M' length times
            elif operation == 'D':                       
                transformed_string.append('-' * length)  # Append '-' for deletions 
            elif operation == 'I':                       
                transformed_string.append('+' * length)  # Append '+' for insertions 

        return ''.join(transformed_string)

    def findTranscriptIndex(self, transformed_string, transcript_nucleotide_position):
        '''
        Finds the index of a specific NUCLEOTIDE NUMBER in the expanded CIGAR representation.

        Parameters:
        -----------
        transformed_string : str
            The expanded representation of the CIGAR string.

        transcript_nucleotide_position : int
            The 0-based position of the nucleotide in the original transcript. 

        Returns:
        --------
        int
            The index of the nucleotide in the transformed string if found; otherwise,
            returns -1 if the nucleotide does not map to a position in the transformed string.
        '''
        if transcript_nucleotide_position > 0: 
            transcript_nucleotide_position -= 1
        
        index_in_original_string = 0
        index_in_transformed_string = 0
        
        for char in transformed_string:
            if char == 'M':
                if index_in_original_string == transcript_nucleotide_position:
                    return index_in_transformed_string 
                index_in_original_string += 1
                index_in_transformed_string += 1 
            elif char == '-':
                index_in_transformed_string += 1
            elif char == '+':
                index_in_original_string += 1
                index_in_transformed_string += 1 
        
        print(f"Most likely, it maps to an insertion")
        return -1

    def calculateGenomePosition(self, transcript_modified_string, transcript_position_modified_string, position_start_transcript_on_genome):
        '''
        Calculates the genome position corresponding to a given nucleotide position in a modified transcript string.

        Parameters:
        -----------
        transcript_modified_string : str
            The modified transcript string composed of characters 'M', '-', and '+'.

        transcript_position_modified_string : int
            The 0-based position in the modified transcript.

        position_start_transcript_on_genome : int
            The starting position of the transcript in the genome.

        Returns:
        --------
        int or None
            The genome position corresponding to the specified transcript position if found; otherwise, returns None.
        '''
        
        current_genome_position = position_start_transcript_on_genome
        current_position_on_alignment = 0

        for char in transcript_modified_string:
            if char == 'M':
                if current_position_on_alignment == transcript_position_modified_string:
                    return current_genome_position 
                current_position_on_alignment += 1    
                current_genome_position += 1          
            elif char == '-':
                if current_position_on_alignment == transcript_position_modified_string:
                    return current_genome_position
                current_position_on_alignment += 1    
                current_genome_position += 1          
            elif char == '+':
                current_position_on_alignment += 1    

        return None


if __name__ == "__main__":


    # Read file 1 
    filename1 = 'input_transcripts_alignments.txt'
    df1 = pd.read_csv(filename1, sep='\t', header=0, names=['transcript_name', 'chromosome', 'transcript_start', 'cigar_string'])
    df1['transcript_start'] = df1['transcript_start'].astype(int)
    df1 = df1.rename(columns={'transcript_start': 'start_transcript_on_genome'})


    # Read file 2
    filename2 = 'input_transcripts_positions.txt'  
    df2 = pd.read_csv(filename2, sep='\t', header=0, names=['transcript_name','transcript_coordinate'])
    df2['transcript_coordinate'] = df2['transcript_coordinate'].astype(int)
    df2 = df2.rename(columns={'transcript_coordinate': 'transcript_index'})


    # Merge these two files 
    combined_data = pd.merge(df1, df2, on='transcript_name', how='outer')


     # Initialize new columns in the combined dataframes to store the computed results
    combined_data['transformed_cigar'] = None
    combined_data['transcript_index_in_cigar'] = None
    combined_data['genome_position'] = None


    # Iterate over each distinct transcript_name
    for transcript_name, transcript_data in combined_data.groupby('transcript_name'):
       
       # Create an instance of the transcriptToGenome class
       transcript_genome = transcriptToGenome()
    

       # Process each row for the current transcript using iterrows() to access row data correctly
       for a_position, a_row in transcript_data.iterrows():         # Unpack position and row data
        
           # Transform the CIGAR string
           transformed_cigar = transcript_genome.transformCigar(a_row['cigar_string'])
        
           # Find the transcript index in the transformed CIGAR string
           transcript_index = transcript_genome.findTranscriptIndex(transformed_cigar, a_row['transcript_index'] + 1)
        
           # Calculate the genome position corresponding to the transcript index
           genome_position = transcript_genome.calculateGenomePosition(
               transformed_cigar,
               transcript_index, 
               a_row['start_transcript_on_genome'])
        

           # Update the combined data frame with computed values
           combined_data.at[a_position, 'transformed_cigar'] = transformed_cigar
           combined_data.at[a_position, 'transcript_index_in_cigar'] = transcript_index
           combined_data.at[a_position, 'genome_position'] = genome_position

    # Print the updated combined data frame to verify the results
    print(combined_data)

    # Print the columns of interest on the screen and in a file 
    print("Results :")
    print(combined_data.loc[:, ["transcript_name", "transcript_index", "chromosome", "genome_position"]])
    combined_data[["transcript_name", "transcript_index", "chromosome", "genome_position"]].to_csv("output_genome_positions_voop.txt", sep='\t', index=False)