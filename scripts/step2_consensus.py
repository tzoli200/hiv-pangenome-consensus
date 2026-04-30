#!/usr/bin/env python3
"""
Pangenome-based consensus sequence generation, v5b.

Builds a consensus from HIV-1 short reads mapped to a pangenome variation
graph (GFA). Mapping against all subtypes at once instead of a single
reference avoids the bias you'd get from picking the wrong one.

Three passes:
  PASS 1: Best alignment selection
    For each read, keep only the highest-scoring alignment (by AS tag,
    then MAPQ). Pangenomes with ~80-90% inter-subtype conservation tend
    to multi-map a lot, so this matters.

  PASS 2: Vote extraction (regular + insertion bases)
    Walk each alignment's CIGAR along the graph path. Match/mismatch
    positions cast one vote at the corresponding MSA column. Insertions
    (bases in the read but not in the graph) get recorded separately,
    anchored to the preceding MSA position.

  PASS 3: Insertion clustering and filtering
    Indel sliding makes the same biological insertion show up at slightly
    different positions across reads, so we group nearby anchors within a
    window (default 5 bp), keep the highest-coverage one per cluster, and
    drop anything below a coverage threshold relative to the surrounding
    regular coverage (default 50%).

  Consensus: majority vote at each covered MSA position, with any surviving
  insertion bases appended after their anchor.

Inputs:
  --gfa              Pangenome variation graph (GFA, P-lines for paths)
  --gaf              Reads mapped to the graph (GAF, with cg:Z: CIGAR)
  --fastq1/--fastq2  Original read files (for sequence retrieval)
  --correspondence   TSV from step 1: sequence_name, graph_pos, msa_pos

Outputs:
  - Consensus FASTA (with insertions)
  - Same again without insertions (for comparison)
  - Optional insertion report TSV (--output-insertions)
"""

import argparse
import re
import gzip
from collections import defaultdict
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass



#  Constants


COMPLEMENT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}


#  Data structures


@dataclass
class NodeInfo:
    """A node in the variation graph."""
    node_id: str
    sequence: str
    length: int


class CoverageMatrix:
    """
    Per-position nucleotide vote counts at MSA coordinates.
    
    Each MSA position tracks how many reads voted for A, T, G, C, or gap.
    The consensus base is the one with the most votes.
    """
    def __init__(self):
        self.data = defaultdict(lambda: {'A': 0, 'T': 0, 'G': 0, 'C': 0, '-': 0})
    
    def add(self, msa_pos: int, base: str) -> bool:
        """Record a vote for `base` at `msa_pos`."""
        if base in 'ATGC-':
            self.data[msa_pos][base] += 1
            return True
        return False
    
    def consensus(self, msa_pos: int, min_cov: int = 1) -> Tuple[Optional[str], int]:
        """Return (consensus_base, total_coverage) at a position."""
        if msa_pos not in self.data:
            return (None, 0)
        counts = self.data[msa_pos]
        total = sum(counts.values())
        if total < min_cov:
            return (None, total)
        best_base = max(counts.keys(), key=lambda b: counts[b])
        if counts[best_base] == 0:
            return (None, total)
        return (best_base, total)
    
    def total_coverage(self, msa_pos: int) -> int:
        """Return total read depth at a position."""
        if msa_pos not in self.data:
            return 0
        return sum(self.data[msa_pos].values())
    
    def get_covered_positions(self) -> set:
        """Return set of all MSA positions with any coverage."""
        return set(self.data.keys())


class InsertionTracker:
    """
    Tracks insertion events anchored to MSA positions.
    
    An insertion at anchor position P with sub-indices 0, 1, 2 means
    3 extra bases were found between MSA positions P and P+1.
    """
    def __init__(self):
        # {anchor_msa_pos: {sub_index: {'A':0, 'T':0, 'G':0, 'C':0}}}
        self.data = defaultdict(lambda: defaultdict(lambda: {'A': 0, 'T': 0, 'G': 0, 'C': 0}))
        self.event_count = 0
    
    def add(self, anchor_msa_pos: int, sub_index: int, base: str):
        """Record an inserted base at anchor position, sub-index."""
        if base in 'ATGC':
            self.data[anchor_msa_pos][sub_index][base] += 1
            self.event_count += 1
    
    def get_anchors(self) -> set:
        """Return all anchor MSA positions with recorded insertions."""
        return set(self.data.keys())
    
    def get_sub_positions(self, anchor: int) -> Dict[int, Dict[str, int]]:
        """Return {sub_index: {base: count}} for an anchor."""
        return dict(self.data[anchor])
    
    def total_coverage_at(self, anchor: int, sub_index: int) -> int:
        """Total insertion reads at a specific anchor + sub-index."""
        if anchor in self.data and sub_index in self.data[anchor]:
            return sum(self.data[anchor][sub_index].values())
        return 0
    
    def consensus_at(self, anchor: int, sub_index: int) -> Tuple[Optional[str], int]:
        """Return (consensus_base, coverage) for an insertion sub-position."""
        if anchor not in self.data or sub_index not in self.data[anchor]:
            return (None, 0)
        counts = self.data[anchor][sub_index]
        total = sum(counts.values())
        if total == 0:
            return (None, 0)
        best = max(counts.keys(), key=lambda b: counts[b])
        return (best, total)



#  Parsing functions


def parse_cigar(cigar_str: str) -> List[Tuple[int, str]]:
    """Parse a CIGAR string into a list of (length, operation) tuples."""
    if not cigar_str:
        return []
    return [(int(m.group(1)), m.group(2)) for m in re.finditer(r'(\d+)([=XIDMNS])', cigar_str)]


def build_cigar_mappings(cigar_ops: List[Tuple[int, str]]) -> Tuple[Dict[int, int], Dict[int, List[int]]]:
    """
    Build path_pos → read_pos mapping AND insertion tracking from CIGAR.
    
    Returns:
        mapping:    {path_pos: read_pos} for matched/mismatched bases
        insertions: {path_pos: [read_pos_0, read_pos_1, ...]} for inserted bases
                    (the insertion occurs BEFORE path_pos, i.e., after path_pos-1)
    """
    mapping = {}
    insertions = {}
    read_pos = path_pos = 0
    
    for length, op in cigar_ops:
        if op in ('=', 'X', 'M'):  # Match / mismatch
            for _ in range(length):
                mapping[path_pos] = read_pos
                read_pos += 1
                path_pos += 1
        elif op == 'I':  # Insertion (extra bases in read, not in graph)
            ins_positions = list(range(read_pos, read_pos + length))
            insertions[path_pos] = ins_positions
            read_pos += length
        elif op == 'D':  # Deletion (graph bases absent in read)
            path_pos += length
        elif op == 'S':  # Soft clip
            read_pos += length
    
    return mapping, insertions


def parse_fastq(fastq_files: List[str]) -> Dict[str, str]:
    """Load read sequences from FASTQ files (supports .gz)."""
    sequences = {}
    for fq in fastq_files:
        if not fq:
            continue
        print(f"  Loading: {fq}")
        opener = gzip.open if fq.endswith('.gz') else open
        mode = 'rt' if fq.endswith('.gz') else 'r'
        with opener(fq, mode) as f:
            while True:
                header = f.readline()
                if not header:
                    break
                seq = f.readline().strip()
                f.readline()  # + line
                f.readline()  # quality line
                name = header[1:].split()[0]
                sequences[name] = seq
    print(f"  Loaded {len(sequences)} reads")
    return sequences


def parse_gfa(gfa_file: str) -> Tuple[Dict[str, NodeInfo], Dict[str, List[Tuple[str, str, int]]]]:
    """
    Parse a GFA file for nodes (S-lines) and paths (P-lines).
    
    Returns:
        nodes: {node_id: NodeInfo}
        paths: {path_name: [(node_id, orientation, length), ...]}
    """
    nodes = {}
    paths = {}
    
    print(f"  Loading GFA: {gfa_file}")
    with open(gfa_file, 'r') as f:
        for line in f:
            if line.startswith('S\t'):
                fields = line.strip().split('\t')
                node_id = fields[1]
                seq = fields[2]
                nodes[node_id] = NodeInfo(node_id, seq, len(seq))
            elif line.startswith('P\t'):
                fields = line.strip().split('\t')
                path_name = fields[1]
                path_str = fields[2]
                
                path_nodes = []
                for segment in path_str.split(','):
                    node_id = segment[:-1]
                    orient = segment[-1]
                    if node_id in nodes:
                        path_nodes.append((node_id, orient, nodes[node_id].length))
                paths[path_name] = path_nodes
    
    print(f"  Loaded {len(nodes)} nodes, {len(paths)} paths")
    return nodes, paths


def load_correspondence(corr_file: str) -> Dict[Tuple[str, int], int]:
    """
    Load the correspondence table mapping (sequence_name, graph_position) → MSA_position.
    
    This table is the bridge between the variation graph coordinate system
    and the multiple sequence alignment coordinate system.
    """
    corr = {}
    print(f"  Loading correspondence: {corr_file}")
    with open(corr_file, 'r') as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split('\t')
            seq_name, graph_pos, msa_pos = parts[0], int(parts[1]), int(parts[2])
            corr[(seq_name, graph_pos)] = msa_pos
    print(f"  Loaded {len(corr)} entries, {len(set(k[0] for k in corr))} sequences")
    return corr


def build_node_to_msa_lookup(gfa_nodes, paths, corr):
    """
    Build a lookup: node_id → {position_in_node: MSA_position}.
    
    """
    print("  Building node → MSA position lookup...")
    
    node_to_msa = {}
    
    for path_name, path_nodes in paths.items():
        # Resolve correspondence table naming
        corr_name = path_name.split('#')[0]
        test_key = (corr_name, 0)
        if test_key not in corr:
            parts = corr_name.rsplit('_', 1)
            if len(parts) > 1:
                corr_name = parts[0]
        
        graph_pos = 0
        for node_id, orient, length in path_nodes:
            if node_id not in node_to_msa:
                node_to_msa[node_id] = {}
                
            for i in range(length):
                if orient == '+':
                    pos_in_node = i
                else:
                    pos_in_node = length - 1 - i
                
                msa_pos = corr.get((corr_name, graph_pos + i))
                if msa_pos is not None:
                    if pos_in_node not in node_to_msa[node_id]:
                        node_to_msa[node_id][pos_in_node] = msa_pos
            
            graph_pos += length
    
    total_positions = sum(len(v) for v in node_to_msa.values())
    print(f"  Mapped {len(node_to_msa)} nodes with {total_positions} positions")
    
    return node_to_msa



#  Alignment processing


def parse_alignment_score(fields: List[str]) -> Tuple[float, int]:
    """Extract alignment score (AS tag) and MAPQ from GAF fields."""
    as_score = 0.0
    mapq = int(fields[11]) if len(fields) > 11 else 0
    
    for field in fields[12:]:
        if field.startswith('AS:f:'):
            as_score = float(field[5:])
            break
        elif field.startswith('AS:i:'):
            as_score = float(field[5:])
            break
    
    return as_score, mapq


def extract_votes_from_alignment(aln_fields, read_seq, gfa_nodes, node_to_msa):
    """
    Extract regular base votes AND insertion votes from a single alignment.
    
    Walks the alignment path through graph nodes, using the CIGAR string
    to determine which read bases correspond to which graph positions,
    then translates graph positions to MSA coordinates via node_to_msa.
    
    Returns:
        votes:     {msa_pos: base} for regular (match/mismatch) positions
        ins_votes: [(anchor_msa_pos, sub_index, base)] for insertion bases
    """
    read_start = int(aln_fields[2])
    read_end = int(aln_fields[3])
    strand = aln_fields[4]
    path_str = aln_fields[5]
    path_start = int(aln_fields[7])
    path_end = int(aln_fields[8])
    
    # Extract CIGAR from optional fields
    cigar = None
    for field in aln_fields[12:]:
        if field.startswith('cg:Z:'):
            cigar = field[5:]
            break
    
    # Reverse complement if on minus strand
    if strand == '-':
        aligned_seq = ''.join(COMPLEMENT.get(b, 'N') for b in reversed(read_seq))
    else:
        aligned_seq = read_seq
    
    cigar_ops = parse_cigar(cigar)
    if cigar_ops:
        path_to_read, cigar_insertions = build_cigar_mappings(cigar_ops)
    else:
        path_to_read = {i: i for i in range(read_end - read_start)}
        cigar_insertions = {}
    
    # Parse the path string (e.g., ">12>34>56" or "<56<34<12")
    path_parts = re.findall(r'([><])([^><]+)', path_str)
    
    votes = {}
    ins_votes = []
    path_offset = 0
    
    for orient_char, node_id in path_parts:
        if node_id not in gfa_nodes:
            continue
        
        node_info = gfa_nodes[node_id]
        node_len = node_info.length
        node_orient = '+' if orient_char == '>' else '-'
        
        if node_id not in node_to_msa:
            path_offset += node_len
            continue
        
        # Determine overlap between alignment window and this node
        overlap_start = max(path_offset, path_start)
        overlap_end = min(path_offset + node_len, path_end)
        
        if overlap_start >= overlap_end:
            path_offset += node_len
            continue
        
        for path_pos in range(overlap_start, overlap_end):
            pos_in_node_orig = path_pos - path_offset
            rel_path_pos = path_pos - path_start
            
            # --- Regular base vote ---
            if rel_path_pos in path_to_read:
                read_pos = path_to_read[rel_path_pos]
                if 0 <= read_pos < len(aligned_seq):
                    nuc = aligned_seq[read_pos]
                    
                    if node_orient == '-':
                        pos_in_node = node_len - 1 - pos_in_node_orig
                        nuc = COMPLEMENT.get(nuc, 'N')
                    else:
                        pos_in_node = pos_in_node_orig
                    
                    if nuc in 'ATGC':
                        msa_pos = node_to_msa[node_id].get(pos_in_node)
                        if msa_pos is not None and msa_pos not in votes:
                            votes[msa_pos] = nuc
            
            # --- Insertion vote (v5b addition) ---
            if rel_path_pos in cigar_insertions:
                if node_orient == '-':
                    anchor_pos_in_node = node_len - 1 - pos_in_node_orig
                else:
                    anchor_pos_in_node = pos_in_node_orig
                
                anchor_msa = node_to_msa[node_id].get(anchor_pos_in_node)
                if anchor_msa is not None:
                    for sub_idx, ins_read_pos in enumerate(cigar_insertions[rel_path_pos]):
                        if 0 <= ins_read_pos < len(aligned_seq):
                            ins_nuc = aligned_seq[ins_read_pos]
                            if node_orient == '-':
                                ins_nuc = COMPLEMENT.get(ins_nuc, 'N')
                            if ins_nuc in 'ATGC':
                                ins_votes.append((anchor_msa, sub_idx, ins_nuc))
        
        path_offset += node_len
    
    return votes, ins_votes



#  Insertion collapse


def cluster_insertions(insertion_tracker: InsertionTracker, 
                       coverage: CoverageMatrix,
                       window: int = 5,
                       threshold: float = 0.5) -> Tuple[List[Tuple[int, str, int]], dict]:
    """
    Cluster nearby insertion positions and collapse to single events.
    
    Indel sliding causes the same biological insertion to appear at
    slightly different positions across reads. This function:
      1. Groups insertion anchors within `window` bp into clusters
      2. For each cluster, selects the anchor with highest total coverage
      3. Calls consensus at each sub-position of the best anchor
      4. Filters by coverage threshold relative to surrounding regular depth
    
    Args:
        insertion_tracker: recorded insertion events from PASS 2
        coverage: regular position coverage matrix
        window: max distance (bp) to group anchors into one cluster
        threshold: min insertion coverage as fraction of surrounding depth
    
    Returns:
        result: [(anchor_msa_pos, consensus_base, sub_index), ...]
        stats:  summary statistics dict
    """
    anchors = sorted(insertion_tracker.get_anchors())
    
    if not anchors:
        return [], {'total_clusters': 0, 'single': 0, 'multi': 0,
                    'collapsed_positions': 0, 'included': 0, 'filtered': 0}
    
    clusters = []
    current_cluster = [anchors[0]]
    
    for anchor in anchors[1:]:
        if anchor - current_cluster[-1] <= window:
            current_cluster.append(anchor)
        else:
            clusters.append(current_cluster)
            current_cluster = [anchor]
    clusters.append(current_cluster)
    
    result = []
    stats = {'total_clusters': len(clusters), 'single': 0, 'multi': 0,
             'collapsed_positions': 0, 'included': 0, 'filtered': 0}
    
    for cluster in clusters:
        if len(cluster) == 1:
            stats['single'] += 1
            best_anchor = cluster[0]
        else:
            stats['multi'] += 1
            stats['collapsed_positions'] += len(cluster) - 1
            # Pick anchor with highest total insertion coverage
            best_anchor = max(cluster, key=lambda a: sum(
                insertion_tracker.total_coverage_at(a, s)
                for s in insertion_tracker.get_sub_positions(a).keys()
            ))
        
        # Compute surrounding regular coverage for threshold
        surrounding_covs = []
        for offset in range(-2, 3):
            pos = best_anchor + offset
            c = coverage.total_coverage(pos)
            if c > 0:
                surrounding_covs.append(c)
        
        avg_surrounding = sum(surrounding_covs) / len(surrounding_covs) if surrounding_covs else 0
        min_ins_cov = avg_surrounding * threshold
        
        sub_positions = insertion_tracker.get_sub_positions(best_anchor)
        for sub_idx in sorted(sub_positions.keys()):
            base, cov = insertion_tracker.consensus_at(best_anchor, sub_idx)
            if base and cov >= max(min_ins_cov, 1):
                result.append((best_anchor, base, sub_idx))
                stats['included'] += 1
            else:
                stats['filtered'] += 1
    
    return result, stats



#  Main pipeline


def main():
    parser = argparse.ArgumentParser(
        description='Pangenome-based consensus generation v5b '
                    '(best alignment + insertion collapse)')
    
    # Required inputs
    parser.add_argument('--gfa', required=True,
                        help='Pangenome variation graph (GFA format)')
    parser.add_argument('--gaf', required=True,
                        help='Reads mapped to graph (GAF format)')
    parser.add_argument('--fastq1', required=True,
                        help='Forward reads (FASTQ, supports .gz)')
    parser.add_argument('--fastq2', default=None,
                        help='Reverse reads (FASTQ, supports .gz)')
    parser.add_argument('--correspondence', required=True,
                        help='Graph-to-MSA coordinate mapping (TSV)')
    
    # Output options
    parser.add_argument('--output-consensus', default='consensus_v5b.fasta',
                        help='Output consensus FASTA (default: consensus_v5b.fasta)')
    parser.add_argument('--output-name', default='consensus_v5b',
                        help='FASTA header name for consensus')
    parser.add_argument('--output-insertions', default=None,
                        help='Write insertion report TSV (optional)')
    
    # Consensus parameters
    parser.add_argument('--min-coverage', type=int, default=1,
                        help='Min read depth to call a base (default: 1; '
                             'note that 1 overshoots GT length by ~5-7%, '
                             '4-5 brings it down to 0.1-0.5%)')
    parser.add_argument('--no-trim', action='store_true',
                        help='Include all covered positions (no 5-prime trimming)')
    parser.add_argument('--start-codon-msa', type=int, default=None,
                        help='Force consensus start at this MSA position')
    
    # Insertion parameters
    parser.add_argument('--insertion-window', type=int, default=5,
                        help='Cluster insertions within this many bp (default: 5)')
    parser.add_argument('--insertion-threshold', type=float, default=0.5,
                        help='Min insertion coverage as fraction of surrounding (default: 0.5)')
    
    args = parser.parse_args()
    
    print("=" * 60)
    print("PANGENOME-BASED CONSENSUS GENERATION v5b")
    print("(Best alignment + insertion capture + sliding collapse)")
    print("=" * 60)
    
    print("\nLoading input files...")
    gfa_nodes, paths = parse_gfa(args.gfa)
    corr = load_correspondence(args.correspondence)
    node_to_msa = build_node_to_msa_lookup(gfa_nodes, paths, corr)
    
    fastq_files = [args.fastq1]
    if args.fastq2:
        fastq_files.append(args.fastq2)
    reads = parse_fastq(fastq_files)
    

    #  PASS 1: Find best alignment per read

    print(f"\n{'='*60}")
    print("PASS 1: Finding best alignment per read")
    print(f"{'='*60}")
    
    best_alignments = {}
    alignment_count = 0
    multi_mapped_reads = 0
    
    with open(args.gaf, 'r') as f:
        for i, line in enumerate(f):
            fields = line.strip().split('\t')
            if len(fields) < 12:
                continue
            
            alignment_count += 1
            read_name = fields[0].split()[0]
            as_score, mapq = parse_alignment_score(fields)
            
            if read_name not in best_alignments:
                best_alignments[read_name] = (as_score, mapq, fields)
            else:
                prev_as, prev_mapq, _ = best_alignments[read_name]
                if as_score > prev_as or (as_score == prev_as and mapq > prev_mapq):
                    best_alignments[read_name] = (as_score, mapq, fields)
                    multi_mapped_reads += 1
            
            if (i+1) % 5000 == 0:
                print(f"  {i+1:,} alignments processed...")
    
    print(f"\n  Total alignments: {alignment_count:,}")
    print(f"  Unique reads: {len(best_alignments):,}")
    print(f"  Multi-mapped (replaced): {multi_mapped_reads:,}")
    

    #  PASS 2: Extract votes (regular + insertions)

    print(f"\n{'='*60}")
    print("PASS 2: Extracting votes from best alignments")
    print(f"{'='*60}")
    
    coverage = CoverageMatrix()
    ins_tracker = InsertionTracker()
    total_bases = 0
    total_ins_bases = 0
    reads_with_votes = 0
    reads_with_insertions = 0
    
    for read_name, (as_score, mapq, fields) in best_alignments.items():
        read_seq = reads.get(read_name)
        if not read_seq:
            continue
        
        votes, ins_votes = extract_votes_from_alignment(
            fields, read_seq, gfa_nodes, node_to_msa)
        
        if votes:
            reads_with_votes += 1
            for msa_pos, base in votes.items():
                coverage.add(msa_pos, base)
                total_bases += 1
        
        if ins_votes:
            reads_with_insertions += 1
            for anchor_msa, sub_idx, base in ins_votes:
                ins_tracker.add(anchor_msa, sub_idx, base)
                total_ins_bases += 1
    
    print(f"  Reads with regular votes: {reads_with_votes:,}")
    print(f"  Total regular bases: {total_bases:,}")
    print(f"  Reads with insertions: {reads_with_insertions:,}")
    print(f"  Total insertion bases: {total_ins_bases:,}")
    print(f"  Unique insertion anchors: {len(ins_tracker.get_anchors()):,}")
    

    #  PASS 3: Cluster and collapse insertions

    print(f"\n{'='*60}")
    print(f"PASS 3: Clustering insertions "
          f"(window={args.insertion_window}bp, "
          f"threshold={args.insertion_threshold})")
    print(f"{'='*60}")
    
    collapsed_insertions, ins_stats = cluster_insertions(
        ins_tracker, coverage, 
        window=args.insertion_window,
        threshold=args.insertion_threshold
    )
    
    print(f"  Total insertion clusters: {ins_stats['total_clusters']}")
    print(f"  Single-position clusters (no sliding): {ins_stats['single']}")
    print(f"  Multi-position clusters (sliding collapsed): {ins_stats['multi']}")
    print(f"  Positions collapsed by clustering: {ins_stats['collapsed_positions']}")
    print(f"  Insertion bases included in consensus: {ins_stats['included']}")
    print(f"  Insertion bases filtered (low coverage): {ins_stats['filtered']}")
    
    # Build insertion lookup: anchor_msa_pos -> [(sub_idx, base), ...]
    insertion_bases = defaultdict(list)
    for anchor_msa, base, sub_idx in collapsed_insertions:
        insertion_bases[anchor_msa].append((sub_idx, base))
    for anchor in insertion_bases:
        insertion_bases[anchor].sort(key=lambda x: x[0])
    

    #  Build consensus sequence

    print(f"\n{'='*60}")
    print("BUILDING CONSENSUS SEQUENCE")
    print(f"{'='*60}")
    
    covered_positions = sorted(coverage.get_covered_positions())
    
    if not covered_positions:
        print("ERROR: No covered positions found!")
        return 1
    
    min_covered = min(covered_positions)
    max_covered = max(covered_positions)
    print(f"  Coverage range: MSA positions {min_covered} - {max_covered}")
    print(f"  Covered positions: {len(covered_positions)}")
    
    if args.no_trim:
        start_pos = min_covered
    elif args.start_codon_msa:
        start_pos = args.start_codon_msa
    else:
        # Auto-detect: find first position with decent coverage
        for pos in covered_positions:
            _, cov = coverage.consensus(pos, min_cov=1)
            if cov >= 10:
                start_pos = max(covered_positions[0], pos - 50)
                break
        else:
            start_pos = covered_positions[0]
    
    print(f"  Start position: {start_pos}")
    
    consensus_parts = []
    gap_count = 0
    base_count = 0
    ins_base_count = 0
    
    for msa_pos in covered_positions:
        if msa_pos < start_pos:
            continue
        
        # Regular base (majority vote)
        base, cov = coverage.consensus(msa_pos, args.min_coverage)
        
        if base is None:
            continue
        
        if base == '-':
            gap_count += 1
            continue
        
        consensus_parts.append(base)
        base_count += 1
        
        # Append any surviving insertion bases after this position
        if msa_pos in insertion_bases:
            for sub_idx, ins_base in insertion_bases[msa_pos]:
                consensus_parts.append(ins_base)
                ins_base_count += 1
    
    result = ''.join(consensus_parts)
    
    # ---- Statistics ----
    print(f"\n  --- Consensus Statistics ---")
    print(f"  Start MSA position: {start_pos}")
    print(f"  End MSA position:   {max_covered}")
    print(f"  Regular bases:      {base_count}")
    print(f"  Insertion bases:    {ins_base_count}")
    print(f"  Gaps skipped:       {gap_count}")
    print(f"  Final length:       {len(result)} bp")
    
    with open(args.output_consensus, 'w') as f:
        f.write(f">{args.output_name}\n")
        for i in range(0, len(result), 80):
            f.write(result[i:i+80] + '\n')
    
    print(f"\n  Wrote: {args.output_consensus}")
    
    if args.output_insertions:
        with open(args.output_insertions, 'w') as f:
            f.write("anchor_msa_pos\tsub_index\tconsensus_base\tcoverage\t"
                    "avg_surrounding\tcluster_size\tdecision\n")
            
            anchors = sorted(ins_tracker.get_anchors())
            if anchors:
                clusters = []
                current = [anchors[0]]
                for a in anchors[1:]:
                    if a - current[-1] <= args.insertion_window:
                        current.append(a)
                    else:
                        clusters.append(current)
                        current = [a]
                clusters.append(current)
                
                for cluster in clusters:
                    best_anchor = max(cluster, key=lambda a: sum(
                        ins_tracker.total_coverage_at(a, s)
                        for s in ins_tracker.get_sub_positions(a).keys()
                    ))
                    
                    surrounding_covs = []
                    for offset in range(-2, 3):
                        c = coverage.total_coverage(best_anchor + offset)
                        if c > 0:
                            surrounding_covs.append(c)
                    avg_surr = sum(surrounding_covs) / len(surrounding_covs) if surrounding_covs else 0
                    
                    for anchor in cluster:
                        subs = ins_tracker.get_sub_positions(anchor)
                        for sub_idx in sorted(subs.keys()):
                            base, cov = ins_tracker.consensus_at(anchor, sub_idx)
                            is_best = (anchor == best_anchor)
                            included = is_best and cov >= max(avg_surr * args.insertion_threshold, 1)
                            f.write(f"{anchor}\t{sub_idx}\t{base}\t{cov}\t"
                                    f"{avg_surr:.1f}\t{len(cluster)}\t"
                                    f"{'INCLUDED' if included else 'COLLAPSED' if not is_best else 'FILTERED'}\n")
        
        print(f"  Wrote: {args.output_insertions}")
    
    no_ins_file = args.output_consensus.replace('.fasta', '_no_ins.fasta')
    no_ins_result = []
    for msa_pos in covered_positions:
        if msa_pos < start_pos:
            continue
        base, cov = coverage.consensus(msa_pos, args.min_coverage)
        if base is None or base == '-':
            continue
        no_ins_result.append(base)
    
    no_ins_seq = ''.join(no_ins_result)
    with open(no_ins_file, 'w') as f:
        f.write(f">{args.output_name}_no_insertions\n")
        for i in range(0, len(no_ins_seq), 80):
            f.write(no_ins_seq[i:i+80] + '\n')
    
    print(f"  Wrote: {no_ins_file} ({len(no_ins_seq)} bp, for comparison)")
    
    print("\n" + "=" * 60)
    print("DONE!")
    print("=" * 60)
    return 0


if __name__ == '__main__':
    main()