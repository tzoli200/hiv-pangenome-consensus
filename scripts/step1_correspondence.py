#!/usr/bin/env python3

"""
Step 1: build a mapping between graph positions (GFA) and MSA columns.

This script takes an aligned FASTA (MSA) and a GFA graph, then tries to
link positions between them based on the paths. It also does a basic
sanity check to see if the sequences actually match.

Originally written for a specific dataset, so some assumptions (especially
around naming) might not generalize perfectly.
"""

import sys
import argparse
import re
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description='Build graph-to-MSA correspondence table')
    parser.add_argument('--msa', required=True, help='Multiple sequence alignment FASTA file')
    parser.add_argument('--gfa', required=True, help='Pangenome graph GFA file')
    parser.add_argument('--output', default='correspondence.tsv', help='Output correspondence table')
    parser.add_argument('--skip-validation', action='store_true', help='Skip sequence validation (not recommended)')
    return parser.parse_args()

def normalize_name(name):
    # normalize things like CON_A(8) -> CON_A_8 (seen in some FASTA headers)
    return name.replace('(', '_').replace(')', '')

def load_msa(msa_file):
    """
    Load MSA from FASTA.

    Bails if the sequences aren't all the same length (it's an alignment).
    Also normalizes names because they don't always match the GFA paths directly.
    """
    print(f"Loading MSA from {msa_file}...")
    
    msa = {}
    current_name = None
    current_seq = []
    
    with open(msa_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_name:
                    seq_name = normalize_name(current_name)
                    msa[seq_name] = ''.join(current_seq).upper()
                
                current_name = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_name:
            seq_name = normalize_name(current_name)
            msa[seq_name] = ''.join(current_seq).upper()
    
    print(f"Loaded {len(msa)} sequences")
    
    # all sequences must be same length, bail otherwise
    lengths = set(len(seq) for seq in msa.values())
    if len(lengths) != 1:
        print(f"ERROR: sequences in MSA have different lengths: {lengths}")
        sys.exit(1)
    
    msa_length = list(lengths)[0]
    print(f"MSA length = {msa_length}")
    
    return msa, msa_length

def load_gfa(gfa_file):
    """
    Parse GFA file.

    Only cares about:
      - S lines (nodes)
      - P/W lines (paths)

    Not handling every possible GFA edge case here.
    """
    print(f"\nLoading GFA from {gfa_file}...")
    
    nodes = {}  # node_id -> sequence
    paths = {}  # path_name -> [(node_id, orientation), ...]
    
    with open(gfa_file, 'r') as f:
        for line in f:
            if line.startswith('S\t'):
                # segment line: node id + sequence
                fields = line.strip().split('\t')
                node_id = fields[1]
                sequence = fields[2].upper()
                nodes[node_id] = sequence
                
            elif line.startswith('P\t'):
                # path line, something like: 1+,2+,3-
                fields = line.strip().split('\t')
                path_name = fields[1]
                path_str = fields[2]
                
                path_nodes = []
                for segment in path_str.split(','):
                    node_id = segment[:-1]
                    orient = segment[-1]  # '+' or '-'
                    path_nodes.append((node_id, orient))
                
                paths[path_name] = path_nodes
                
            elif line.startswith('W\t'):
                # alternative path format (walk)
                fields = line.strip().split('\t')
                path_name = f"{fields[1]}#{fields[2]}#{fields[3]}"
                walk_str = fields[6]
                
                # walk looks like >1>2<3 etc.
                path_nodes = []
                parts = re.split(r'([><])', walk_str)
                orient = '+'
                for part in parts:
                    if part == '>':
                        orient = '+'
                    elif part == '<':
                        orient = '-'
                    elif part:
                        path_nodes.append((part, orient))
                
                paths[path_name] = path_nodes
    
    print(f"Loaded {len(nodes)} nodes")
    print(f"Loaded {len(paths)} paths")
    
    return nodes, paths

def get_reverse_complement(seq):
    # simple reverse complement; non-ACGTN silently becomes N
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement.get(b, 'N') for b in reversed(seq))

def extract_path_sequence(path_nodes, nodes):
    """
    Build full sequence for a path by concatenating node sequences.

    If a node is missing, just skip it (this shouldn't really happen,
    but seen it once during testing).
    """
    seq_parts = []
    
    for node_id, orient in path_nodes:
        if node_id not in nodes:
            print(f"WARNING: node {node_id} not found")
            continue
            
        node_seq = nodes[node_id]
        
        if orient == '-':
            node_seq = get_reverse_complement(node_seq)
        
        seq_parts.append(node_seq)
    
    return ''.join(seq_parts)

def find_matching_msa_name(path_name, msa_names):
    """
    Try to match a GFA path name to an MSA sequence.

    This part is a bit messy because naming isn't consistent between files.
    Covers a few cases seen in the current dataset, but not guaranteed to work
    universally.
    """
    if path_name in msa_names:
        return path_name
    
    normalized = normalize_name(path_name)
    if normalized in msa_names:
        return normalized
    
    # try removing common suffixes
    for suffix in ['#0#genome', '_#0#genome', '#0', '_0']:
        stripped = path_name.replace(suffix, '')
        if stripped in msa_names:
            return stripped
        stripped_norm = normalize_name(stripped)
        if stripped_norm in msa_names:
            return stripped_norm
    
    # last attempt: partial matches (loose, can pick wrong name)
    for msa_name in msa_names:
        if msa_name in path_name or path_name in msa_name:
            return msa_name
        
        msa_norm = normalize_name(msa_name)
        path_norm = normalize_name(path_name)
        if msa_norm in path_norm or path_norm in msa_norm:
            return msa_name
    
    return None

def build_correspondence_with_validation(msa, nodes, paths, skip_validation=False):
    """
    Build mapping between path sequence positions and MSA columns.

    Reconstructs sequences from GFA paths and compares them to the ungapped
    MSA. The mapping itself is derived from the MSA (gaps removed).

    If they don't match, continues anyway (results might suffer from this tho).
    """
    print("\nBuilding correspondence table...")
    
    correspondence = {}
    msa_names = set(msa.keys())
    matched_paths = 0
    
    for path_name, path_nodes in paths.items():
        msa_name = find_matching_msa_name(path_name, msa_names)
        
        if msa_name is None:
            print(f"No match for path '{path_name}'")
            continue
        
        print(f"Processing {path_name} -> {msa_name}")
        
        path_seq = extract_path_sequence(path_nodes, nodes)
        
        msa_seq_gapped = msa[msa_name]
        msa_seq_ungapped = msa_seq_gapped.replace('-', '')
        
        if not skip_validation:
            if path_seq != msa_seq_ungapped:
                print("  sequence mismatch")
                print(f"  GFA length: {len(path_seq)}")
                print(f"  MSA length: {len(msa_seq_ungapped)}")
                
                for i, (a, b) in enumerate(zip(path_seq, msa_seq_ungapped)):
                    if a != b:
                        print(f"  first diff at {i}: {a} vs {b}")
                        break
                
                if len(path_seq) != len(msa_seq_ungapped):
                    print("  lengths differ")
                
                print("  using MSA positions anyway")
        
        graph_pos = 0
        positions_mapped = 0
        
        for msa_col in range(len(msa_seq_gapped)):
            base = msa_seq_gapped[msa_col]
            
            if base != '-':
                msa_position = msa_col + 1
                correspondence[(msa_name, graph_pos)] = msa_position
                graph_pos += 1
                positions_mapped += 1
        
        print(f"  mapped {positions_mapped} positions")
        matched_paths += 1
    
    print(f"\nMatched {matched_paths}/{len(paths)} paths")
    
    return correspondence

def write_correspondence(correspondence, output_file):
    print(f"\nWriting output to {output_file}...")
    
    with open(output_file, 'w') as f:
        f.write("sequence_name\tgraph_position\tmsa_position\n")
        
        for (seq_name, graph_pos), msa_pos in sorted(correspondence.items()):
            f.write(f"{seq_name}\t{graph_pos}\t{msa_pos}\n")
    
    print(f"Wrote {len(correspondence)} rows")

def main():
    args = parse_args()
    
    print("=" * 60)
    print("STEP 1: Graph-to-MSA mapping")
    print("=" * 60)
    
    msa, msa_length = load_msa(args.msa)
    nodes, paths = load_gfa(args.gfa)
    
    correspondence = build_correspondence_with_validation(
        msa, nodes, paths, 
        skip_validation=args.skip_validation
    )
    
    write_correspondence(correspondence, args.output)
    
    print("\nDone.")

if __name__ == "__main__":
    main()