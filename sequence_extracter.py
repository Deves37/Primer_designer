def parse_blastn_alignment(alignment_file):
    """Parse blastn alignment file and extract query and subject sequences"""
    
    query_lines = []
    subject_lines = []
    
    with open(alignment_file, 'r') as f:
        lines = [line.strip() for line in f.readlines() if line.strip()]
    
    i = 0
    while i < len(lines):
        current_line = lines[i]
        
        # Handle Query line
        if current_line.startswith('Query'):
            parts = current_line.split()
            if len(parts) >= 3:
                query_seq = parts[2]
                query_lines.append(query_seq)
            
            # Look ahead for subject line (could be immediately after or after match line)
            if i + 1 < len(lines):
                next_line = lines[i + 1]
                # If next line is a match line, then subject should be after that
                if any(char in next_line for char in ['|', ' ', '+', '.']) and not next_line.startswith('Sbjct'):
                    if i + 2 < len(lines) and lines[i + 2].startswith('Sbjct'):
                        subject_parts = lines[i + 2].split()
                        if len(subject_parts) >= 3:
                            subject_seq = subject_parts[2]
                            subject_lines.append(subject_seq)
                        i += 2  # Skip match and subject lines
                # If next line is directly a subject line (no match line)
                elif i + 1 < len(lines) and lines[i + 1].startswith('Sbjct'):
                    subject_parts = lines[i + 1].split()
                    if len(subject_parts) >= 3:
                        subject_seq = subject_parts[2]
                        subject_lines.append(subject_seq)
                    i += 1  # Skip subject line
        
        # Handle standalone Subject line (in case Query line was missing)
        elif current_line.startswith('Sbjct'):
            parts = current_line.split()
            if len(parts) >= 3:
                subject_seq = parts[2]
                subject_lines.append(subject_seq)
        
        i += 1
    
    # Combine all sequence parts and convert to uppercase
    query_sequence = ''.join(query_lines).upper()
    subject_sequence = ''.join(subject_lines).upper()
    
    return query_sequence, subject_sequence

def write_fasta(filename, header, sequence):
    """Write sequence to FASTA file with proper formatting"""
    with open(filename, 'w') as f:
        f.write(header + '\n')
        
        # Write sequence in chunks of 60 characters
        for i in range(0, len(sequence), 60):
            f.write(sequence[i:i+60] + '\n')

def main():
    # Parse the alignment file
    query_seq, subject_seq = parse_blastn_alignment('blastn_aligned.txt')
    
    # Write query to FASTA
    write_fasta('Query_pig.fasta', '>Query', query_seq)
    
    # Write subject to FASTA
    write_fasta('Subject_human.fasta', '>Subject', subject_seq)
    
    print("FASTA files created successfully!")
    print(f"Query sequence length: {len(query_seq)}")
    print(f"Subject sequence length: {len(subject_seq)}")
    
    # Print the last few characters to verify the end was captured
    print(f"\nLast 20 chars of query: {query_seq[-20:]}")
    print(f"Last 20 chars of subject: {subject_seq[-20:]}")

if __name__ == "__main__":
    main()