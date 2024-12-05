#!/usr/bin/env python 
    
def detect_quality_encoding(quality_scores):
    ascii_values = [ord(q) for q in quality_scores]
    min_q = min(ascii_values)
    max_q = max(ascii_values)
    
    # Define possible encoding schemes
    encodings = {
        'Phred+33 (Sanger, Illumina 1.8+)': range(33, 74),
        'Phred+64 (Illumina 1.3+)': range(64, 104),
        'Solexa+64 (Solexa/Illumina 1.0)': range(59, 104),
        'Illumina 1.5+ (Phred+64)': range(66, 104),
    }
    
    possible_encodings = []
    for encoding, ascii_range in encodings.items():
        if min_q in ascii_range and max_q in ascii_range:
            possible_encodings.append(encoding)
    
    if not possible_encodings:
        return 'Unknown'
    elif len(possible_encodings) == 1:
        return possible_encodings[0]
    else:
        return ', '.join(possible_encodings)

def analyze_fastq(fastq_path, num_reads=1000):
    """
    Analyze the FASTQ file to determine quality encoding and read length distribution.
    """
    # Handle gzipped files
    if fastq_path.endswith('.gz'):
        handle = gzip.open(fastq_path, 'rt')
    else:
        handle = open(fastq_path, 'r')
    
    read_lengths = []
    quality_chars = []
    count = 0
    
    try:
        for record in SeqIO.parse(handle, 'fastq'):
            read_lengths.append(len(record.seq))
            quality_chars.extend([chr(q + 33) for q in record.letter_annotations["phred_quality"]])
            count += 1
            if count >= num_reads:
                break
    except Exception as e:
        print(f"Error reading FASTQ file: {e}")
        sys.exit(1)
    finally:
        handle.close()
    
    # Detect quality encoding
    encoding = detect_quality_encoding(quality_chars)
    
    # Analyze read lengths
    length_counter = Counter(read_lengths)
    most_common_lengths = length_counter.most_common(10)
    
    return encoding, read_lengths, most_common_lengths

def infer_sequencing_technology(encoding, read_lengths):
    """
    Infer possible sequencing technologies based on encoding and read length.
    """
    # Define typical characteristics of sequencing technologies
    technologies = [
        {
            'name': 'Illumina',
            'encodings': ['Phred+33 (Sanger, Illumina 1.8+)'],
            'read_length': (50, 300)
        },
        {
            'name': 'Illumina (Phred+64)',
            'encodings': ['Phred+64 (Illumina 1.3+)', 'Illumina 1.5+ (Phred+64)'],
            'read_length': (50, 300)
        },
        {
            'name': 'Ion Torrent',
            'encodings': ['Phred+33 (Sanger, Illumina 1.8+)', 'Phred+64 (Illumina 1.3+)'],
            'read_length': (200, 400)
        },
        {
            'name': 'PacBio',
            'encodings': ['Unknown'],  # PacBio uses different quality metrics
            'read_length': (1000, 30000)
        },
        {
            'name': 'Oxford Nanopore',
            'encodings': ['Unknown'],  # Nanopore uses different quality metrics
            'read_length': (500, 200000)
        },
    ]
    
    # Determine read length statistics
    avg_length = sum(read_lengths) / len(read_lengths)
    median_length = sorted(read_lengths)[len(read_lengths) // 2]
    
    possible_techs = []
    for tech in technologies:
        if encoding in tech['encodings']:
            min_len, max_len = tech['read_length']
            if min_len <= avg_length <= max_len:
                possible_techs.append(tech['name'])
    
    # Handle technologies with unknown encoding but typical long read lengths
    if not possible_techs:
        max_len = max(read_lengths)
        if max_len > 10000:
            possible_techs.append('PacBio or Oxford Nanopore')
        elif 500 <= max_len <= 10000:
            possible_techs.append('Oxford Nanopore')
    
    if not possible_techs:
        possible_techs.append('Unknown or other technology')
    
    return possible_techs, avg_length, median_length

def main():
    ap = argparse.ArgumentParser(
        description='Infer sequencing technology from FASTQ file based on Phred scores and read lengths.'
    )
    ap.add_argument(
        '--fastq',
        type = str,
        required = True,
        help='Path to FASTQ file (can be gzipped)'
    )
    ap.add_argument('--num_reads', 
        type=int, 
        required = False,
        default=1000,
        help='Number of reads to analyze (default: 1000)'
    )
    args = ap.parse_args()


    encoding, read_lengths, most_common_lengths = analyze_fastq(args.fastq, args.num_reads)
    
    print(f"Detected Quality Encoding: {encoding}")
    print("\nMost Common Read Lengths:")
    for length, count in most_common_lengths:
        print(f"  Length {length}: {count} reads")
    
    possible_techs, avg_len, median_len = infer_sequencing_technology(encoding, read_lengths)
    
    print(f"\nRead Length Statistics:")
    print(f"  Average Length: {avg_len:.2f}")
    print(f"  Median Length: {median_len}")
    
    print("\nPossible Sequencing Technologies:")
    for tech in possible_techs:
        print(f"  - {tech}")

if __name__ == '__main__':
    import argparse, gzip
    from collections import defaultdict, Counter
    from Bio import SeqIO
    main()
