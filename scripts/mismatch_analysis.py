#!/usr/bin/env python3
"""
Compare two consensus calls (a node-based "v5" and Shiver) against a
ground truth, position by position, and report where each one got things
wrong. Outputs a TSV of mismatch positions, a summary printout, and a
standalone HTML viz with a genome map and filterable table.

If the three sequences aren't already aligned, MAFFT gets called (tries
Docker first, falls back to a local install).

Originally set up for HIV-1, so the HTML viz uses HXB2-based gene
boundaries (approximate landmarks for other Group M subtypes, not exact).
"""

import argparse
import os
import sys
import subprocess
import tempfile
from collections import defaultdict



# FASTA I/O


def load_fasta(filepath):
    """Load a FASTA file, return dict {name: sequence}."""
    seqs = {}
    name = None
    seq_parts = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if name is not None:
                    seqs[name] = ''.join(seq_parts).upper()
                name = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)
        if name is not None:
            seqs[name] = ''.join(seq_parts).upper()
    return seqs


def write_fasta(seqs, filepath):
    """Write sequences dict to FASTA file."""
    with open(filepath, 'w') as f:
        for name, seq in seqs.items():
            f.write(f'>{name}\n')
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + '\n')



# ALIGNMENT


def align_sequences(seq_dict, output_path):
    """Align sequences using MAFFT (try Docker first, then local)."""
    # Use absolute paths to avoid Docker mount issues
    abs_output = os.path.abspath(output_path)
    abs_dir = os.path.dirname(abs_output)
    os.makedirs(abs_dir, exist_ok=True)

    tmp_input = os.path.join(abs_dir, 'tmp_mafft_input.fasta')
    write_fasta(seq_dict, tmp_input)

    # try Docker MAFFT first (the -v flag mounts the dir holding the temp file)
    cmd_docker = (
        f"docker run --rm -v {abs_dir}:/data biocontainers/mafft:v7.407-2-deb_cv1 "
        f"mafft --auto /data/tmp_mafft_input.fasta"
    )
    # Try local MAFFT
    cmd_local = f"mafft --auto {tmp_input}"

    aligned = None
    for cmd in [cmd_docker, cmd_local]:
        try:
            print(f"  Trying: {cmd[:80]}...")
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=120)
            if result.returncode == 0 and result.stdout.strip():
                with open(abs_output, 'w') as f:
                    f.write(result.stdout)
                aligned = load_fasta(abs_output)
                print(f"  MAFFT success! Aligned {len(aligned)} sequences.")
                break
            else:
                print(f"  Failed (rc={result.returncode}): {result.stderr[:200]}")
        except Exception as e:
            print(f"  Exception: {e}")
            continue

    # Cleanup
    if os.path.exists(tmp_input):
        os.remove(tmp_input)

    if aligned is None:
        print("ERROR: Could not run MAFFT. Install it or use Docker.", file=sys.stderr)
        sys.exit(1)

    return aligned


# MISMATCH CLASSIFICATION


PURINES = {'A', 'G'}
PYRIMIDINES = {'C', 'T'}

def classify_mutation(ref_base, alt_base):
    """Classify a substitution as transition or transversion."""
    if ref_base == alt_base:
        return 'match'
    if ref_base not in 'ATGC' or alt_base not in 'ATGC':
        return 'ambiguous'
    if {ref_base, alt_base} <= PURINES or {ref_base, alt_base} <= PYRIMIDINES:
        return 'transition'
    return 'transversion'


def get_mutation_label(ref_base, alt_base):
    """Get specific mutation label like A>G."""
    return f"{ref_base}>{alt_base}"



# CORE ANALYSIS


def analyze_mismatches(gt_seq, v5_seq, shiver_seq):
    """
    Walk through the three aligned sequences in lockstep and record every
    position where v5 or Shiver disagrees with the ground truth.

    Skips positions where GT itself is a gap (those are insertions in one
    of the other methods, not really "mismatches" in any useful sense).
    Returns a list of dicts, one per disagreement.
    """
    assert len(gt_seq) == len(v5_seq) == len(shiver_seq), \
        f"Aligned sequences must be same length: GT={len(gt_seq)}, v5={len(v5_seq)}, Shiver={len(shiver_seq)}"

    results = []
    genome_pos = 0  # ungapped position in ground truth

    for aln_pos in range(len(gt_seq)):
        gt_base = gt_seq[aln_pos]
        v5_base = v5_seq[aln_pos]
        sh_base = shiver_seq[aln_pos]

        # Track genome position (ungapped GT)
        if gt_base in 'ATGCN':
            genome_pos += 1

        # Skip positions where GT is a gap (insertion in one of the methods)
        if gt_base not in 'ATGC':
            continue

        v5_match = (v5_base == gt_base)
        sh_match = (sh_base == gt_base)

        # Only record positions where at least one method disagrees with GT
        if not v5_match or not sh_match:
            if not v5_match and not sh_match:
                if v5_base == sh_base:
                    category = 'shared_same_error'
                else:
                    category = 'shared_diff_error'
            elif not v5_match:
                category = 'v5_only'
            else:
                category = 'shiver_only'

            record = {
                'aln_pos': aln_pos + 1,
                'genome_pos': genome_pos,
                'gt_base': gt_base,
                'v5_base': v5_base if v5_base in 'ATGC-' else v5_base,
                'shiver_base': sh_base if sh_base in 'ATGC-' else sh_base,
                'v5_correct': v5_match,
                'shiver_correct': sh_match,
                'category': category,
            }

            if not v5_match and v5_base in 'ATGC':
                record['v5_mut_type'] = classify_mutation(gt_base, v5_base)
                record['v5_mut_label'] = get_mutation_label(gt_base, v5_base)
            else:
                record['v5_mut_type'] = 'gap/N' if not v5_match else 'correct'
                record['v5_mut_label'] = f"{gt_base}>{v5_base}" if not v5_match else '-'

            if not sh_match and sh_base in 'ATGC':
                record['shiver_mut_type'] = classify_mutation(gt_base, sh_base)
                record['shiver_mut_label'] = get_mutation_label(gt_base, sh_base)
            else:
                record['shiver_mut_type'] = 'gap/N' if not sh_match else 'correct'
                record['shiver_mut_label'] = f"{gt_base}>{sh_base}" if not sh_match else '-'

            results.append(record)

    return results



# OUTPUT: TSV TABLE


def write_tsv(results, filepath, label):
    """Write mismatch table to TSV."""
    header = [
        'dataset', 'genome_pos', 'aln_pos', 'gt_base',
        'v5_base', 'v5_correct', 'v5_mut_type', 'v5_mut_label',
        'shiver_base', 'shiver_correct', 'shiver_mut_type', 'shiver_mut_label',
        'category'
    ]
    with open(filepath, 'w') as f:
        f.write('\t'.join(header) + '\n')
        for r in results:
            row = [
                label,
                str(r['genome_pos']),
                str(r['aln_pos']),
                r['gt_base'],
                r['v5_base'],
                str(r['v5_correct']),
                r['v5_mut_type'],
                r['v5_mut_label'],
                r['shiver_base'],
                str(r['shiver_correct']),
                r['shiver_mut_type'],
                r['shiver_mut_label'],
                r['category'],
            ]
            f.write('\t'.join(row) + '\n')



# OUTPUT: SUMMARY STATISTICS


def print_summary(results, label, gt_length):
    """
    Print a summary block to stdout. Also returns a small dict with the
    headline numbers (accuracies, error counts) for downstream use.
    """
    v5_errors = [r for r in results if not r['v5_correct']]
    sh_errors = [r for r in results if not r['shiver_correct']]
    shared = [r for r in results if not r['v5_correct'] and not r['shiver_correct']]
    shared_same = [r for r in shared if r['category'] == 'shared_same_error']
    shared_diff = [r for r in shared if r['category'] == 'shared_diff_error']
    v5_only = [r for r in results if r['category'] == 'v5_only']
    sh_only = [r for r in results if r['category'] == 'shiver_only']

    v5_acc = (1 - len(v5_errors) / gt_length) * 100
    sh_acc = (1 - len(sh_errors) / gt_length) * 100

    print(f"\n{'='*60}")
    print(f"  MISMATCH ANALYSIS: {label}")
    print(f"{'='*60}")
    print(f"  Ground truth length: {gt_length} bp")
    print(f"")
    print(f"  v5 accuracy:     {v5_acc:.2f}%  ({len(v5_errors)} mismatches)")
    print(f"  Shiver accuracy: {sh_acc:.2f}%  ({len(sh_errors)} mismatches)")
    print(f"")
    print(f"  --- Error categories ---")
    print(f"  Shared errors (same call):  {len(shared_same)}")
    print(f"  Shared errors (diff call):  {len(shared_diff)}")
    print(f"  v5-only errors:             {len(v5_only)}")
    print(f"  Shiver-only errors:         {len(sh_only)}")
    print(f"  Total unique error positions: {len(results)}")

    # Mutation type breakdown
    for method_name, errors in [('v5', v5_errors), ('Shiver', sh_errors)]:
        ti = sum(1 for r in errors if r[f'{method_name.lower() if method_name == "v5" else "shiver"}_mut_type'] == 'transition')
        tv = sum(1 for r in errors if r[f'{method_name.lower() if method_name == "v5" else "shiver"}_mut_type'] == 'transversion')
        gap = sum(1 for r in errors if r[f'{method_name.lower() if method_name == "v5" else "shiver"}_mut_type'] == 'gap/N')
        key = 'v5' if method_name == 'v5' else 'shiver'
        ti = sum(1 for r in errors if r[f'{key}_mut_type'] == 'transition')
        tv = sum(1 for r in errors if r[f'{key}_mut_type'] == 'transversion')
        gap = sum(1 for r in errors if r[f'{key}_mut_type'] == 'gap/N')
        ti_ratio = ti / (ti + tv) if (ti + tv) > 0 else 0
        print(f"\n  {method_name} mutation breakdown:")
        print(f"    Transitions:    {ti}  ({ti_ratio*100:.0f}% of substitutions)")
        print(f"    Transversions:  {tv}")
        print(f"    Gaps/N:         {gap}")

    # List specific mutation types
    for method_name in ['v5', 'shiver']:
        key = method_name
        errors = [r for r in results if not r[f'{key}_correct']]
        mut_counts = defaultdict(int)
        for r in errors:
            lbl = r[f'{key}_mut_label']
            if lbl != '-':
                mut_counts[lbl] += 1
        if mut_counts:
            print(f"\n  {method_name} specific mutations:")
            for mut, cnt in sorted(mut_counts.items(), key=lambda x: -x[1]):
                print(f"    {mut}: {cnt}")

    print(f"\n{'='*60}")

    return {
        'label': label,
        'gt_length': gt_length,
        'v5_errors': len(v5_errors),
        'shiver_errors': len(sh_errors),
        'shared_same': len(shared_same),
        'shared_diff': len(shared_diff),
        'v5_only': len(v5_only),
        'shiver_only': len(sh_only),
        'v5_acc': v5_acc,
        'shiver_acc': sh_acc,
    }



# OUTPUT: HTML VISUALIZATION


def generate_html_visualization(results, label, gt_length, filepath):
    """Generate an interactive HTML visualization of mismatch positions."""

    v5_only = [r for r in results if r['category'] == 'v5_only']
    sh_only = [r for r in results if r['category'] == 'shiver_only']
    shared_same = [r for r in results if r['category'] == 'shared_same_error']
    shared_diff = [r for r in results if r['category'] == 'shared_diff_error']

    v5_errors = [r for r in results if not r['v5_correct']]
    sh_errors = [r for r in results if not r['shiver_correct']]
    v5_acc = (1 - len(v5_errors) / gt_length) * 100
    sh_acc = (1 - len(sh_errors) / gt_length) * 100

    # HIV genome regions (approximate positions based on HXB2)
    hiv_regions = [
        ("5'LTR", 1, 634),
        ("gag", 790, 2292),
        ("pol", 2085, 5096),
        ("vif", 5041, 5619),
        ("vpr", 5559, 5850),
        ("tat", 5831, 6045),
        ("rev", 5970, 6045),
        ("vpu", 6062, 6310),
        ("env", 6225, 8795),
        ("nef", 8797, 9417),
        ("3'LTR", 9086, 9719),
    ]

    # Build mismatch data as JS arrays
    def to_js_array(records, method_key=None):
        items = []
        for r in records:
            tooltip_parts = [
                f"Pos: {r['genome_pos']}",
                f"GT: {r['gt_base']}",
                f"v5: {r['v5_base']}",
                f"Shiver: {r['shiver_base']}",
            ]
            if method_key and r.get(f'{method_key}_mut_type'):
                tooltip_parts.append(f"Type: {r[f'{method_key}_mut_type']}")
            tooltip = ' | '.join(tooltip_parts)
            items.append(f"{{pos:{r['genome_pos']},tip:'{tooltip}'}}")
        return '[' + ','.join(items) + ']'

    regions_js = '[' + ','.join(
        f'{{name:"{name}",start:{s},end:{e}}}'
        for name, s, e in hiv_regions
    ) + ']'

    html = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="UTF-8">
<title>Mismatch Analysis: {label}</title>
<style>
  * {{ margin: 0; padding: 0; box-sizing: border-box; }}
  body {{ font-family: 'Segoe UI', Arial, sans-serif; background: #0f0f0f; color: #e0e0e0; padding: 20px; }}
  h1 {{ color: #fff; margin-bottom: 5px; font-size: 1.6em; }}
  .subtitle {{ color: #888; margin-bottom: 8px; font-size: 0.95em; }}
  .stats-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 12px; margin-bottom: 25px; }}
  .stat-card {{ background: #1a1a2e; border: 1px solid #333; border-radius: 8px; padding: 15px; text-align: center; }}
  .stat-card .value {{ font-size: 1.8em; font-weight: bold; margin: 5px 0; }}
  .stat-card .label {{ font-size: 0.85em; color: #999; }}
  .stat-card.v5 .value {{ color: #4fc3f7; }}
  .stat-card.shiver .value {{ color: #ff8a65; }}
  .stat-card.shared .value {{ color: #ce93d8; }}
  .stat-card.v5only .value {{ color: #4fc3f7; }}
  .stat-card.shonly .value {{ color: #ff8a65; }}

  .chart-container {{ background: #1a1a2e; border: 1px solid #333; border-radius: 8px; padding: 20px; margin-bottom: 20px; }}
  .chart-title {{ font-size: 1.1em; margin-bottom: 10px; color: #fff; }}
  canvas {{ width: 100%; }}
  
  .legend {{ display: flex; gap: 20px; flex-wrap: wrap; margin: 10px 0; }}
  .legend-item {{ display: flex; align-items: center; gap: 6px; font-size: 0.85em; }}
  .legend-dot {{ width: 12px; height: 12px; border-radius: 50%; }}

  .disclaimer-box {{ background: #1a1a2e; border: 1px solid #444; border-left: 3px solid #ffd54f; border-radius: 6px; padding: 14px 18px; margin-bottom: 20px; font-size: 0.85em; line-height: 1.5; }}
  .disclaimer-box .disc-title {{ color: #ffd54f; font-weight: bold; margin-bottom: 6px; font-size: 0.95em; }}
  .disclaimer-box a {{ color: #4fc3f7; text-decoration: none; }}
  .disclaimer-box a:hover {{ text-decoration: underline; }}
  .disclaimer-box .gene-table {{ margin-top: 10px; border-collapse: collapse; width: 100%; }}
  .disclaimer-box .gene-table th {{ background: #16213e; color: #aaa; padding: 4px 8px; text-align: left; font-weight: normal; font-size: 0.9em; }}
  .disclaimer-box .gene-table td {{ padding: 3px 8px; border-bottom: 1px solid #222; font-family: 'Consolas', 'Courier New', monospace; font-size: 0.9em; }}

  .table-container {{ background: #1a1a2e; border: 1px solid #333; border-radius: 8px; padding: 20px; margin-bottom: 20px; overflow-x: auto; }}
  table {{ border-collapse: collapse; width: 100%; font-size: 0.85em; }}
  th {{ background: #16213e; color: #ccc; padding: 8px 10px; text-align: left; position: sticky; top: 0; }}
  td {{ padding: 6px 10px; border-bottom: 1px solid #222; }}
  tr:hover td {{ background: #1f2940; }}
  .cat-shared_same_error {{ color: #ce93d8; }}
  .cat-shared_diff_error {{ color: #f48fb1; }}
  .cat-v5_only {{ color: #4fc3f7; }}
  .cat-shiver_only {{ color: #ff8a65; }}
  .mut-transition {{ color: #81c784; }}
  .mut-transversion {{ color: #e57373; }}

  .filter-bar {{ margin-bottom: 10px; display: flex; gap: 10px; flex-wrap: wrap; }}
  .filter-btn {{ background: #16213e; border: 1px solid #444; color: #ccc; padding: 5px 12px; border-radius: 4px; cursor: pointer; font-size: 0.85em; }}
  .filter-btn.active {{ border-color: #4fc3f7; color: #4fc3f7; }}
  .filter-btn:hover {{ border-color: #888; }}

  .mutation-grid {{ display: grid; grid-template-columns: 1fr 1fr; gap: 15px; margin-bottom: 20px; }}
  .mutation-box {{ background: #1a1a2e; border: 1px solid #333; border-radius: 8px; padding: 15px; }}
  .mutation-box h3 {{ font-size: 1em; margin-bottom: 8px; }}
  .mut-bar {{ display: flex; align-items: center; gap: 8px; margin: 4px 0; font-size: 0.85em; }}
  .mut-bar-fill {{ height: 16px; border-radius: 3px; min-width: 2px; }}
  .mut-bar-label {{ min-width: 50px; }}
</style>
</head>
<body>

<h1>Mismatch Analysis: {label}</h1>
<p class="subtitle">v5 (node-based consensus) vs Shiver vs Ground Truth</p>

<!-- Summary stats -->
<div class="stats-grid">
  <div class="stat-card v5">
    <div class="label">v5 Accuracy</div>
    <div class="value">{v5_acc:.2f}%</div>
    <div class="label">{len(v5_errors)} mismatches</div>
  </div>
  <div class="stat-card shiver">
    <div class="label">Shiver Accuracy</div>
    <div class="value">{sh_acc:.2f}%</div>
    <div class="label">{len(sh_errors)} mismatches</div>
  </div>
  <div class="stat-card shared">
    <div class="label">Shared Errors</div>
    <div class="value">{len(shared_same) + len(shared_diff)}</div>
    <div class="label">{len(shared_same)} same call, {len(shared_diff)} different</div>
  </div>
  <div class="stat-card v5only">
    <div class="label">v5-Only Errors</div>
    <div class="value">{len(v5_only)}</div>
    <div class="label">errors unique to v5</div>
  </div>
  <div class="stat-card shonly">
    <div class="label">Shiver-Only Errors</div>
    <div class="value">{len(sh_only)}</div>
    <div class="label">errors unique to Shiver</div>
  </div>
</div>

<!-- Genome-wide mismatch map -->
<div class="chart-container">
  <div class="chart-title">Genome-Wide Mismatch Map</div>
  <div class="legend">
    <div class="legend-item"><div class="legend-dot" style="background:#ce93d8"></div>Shared error</div>
    <div class="legend-item"><div class="legend-dot" style="background:#4fc3f7"></div>v5-only error</div>
    <div class="legend-item"><div class="legend-dot" style="background:#ff8a65"></div>Shiver-only error</div>
  </div>
  <canvas id="genomeMap" height="120"></canvas>
</div>

<!-- Gene Region Source Disclaimer -->
<div class="disclaimer-box">
  <div class="disc-title">&#9432; HIV-1 Genome Region Annotation Source</div>
  Gene region boundaries displayed in the genome map and the "HIV Region" column are based on the
  <strong>HXB2 reference genome</strong> (GenBank accession:
  <a href="https://www.ncbi.nlm.nih.gov/nuccore/K03455.1" target="_blank">K03455.1</a>),
  the standard HIV-1 coordinate reference used by the
  <a href="https://www.hiv.lanl.gov/content/sequence/HIV/MAP/landmark.html" target="_blank">Los Alamos HIV Sequence Database</a>.
  Coordinates were verified against the GenBank FEATURES table for K03455 and cross-referenced with
  Korber et al., "Numbering Positions in HIV Relative to HXB2CG",
  <em>Human Retroviruses and AIDS</em>, 1998
  (<a href="https://www.hiv.lanl.gov/content/sequence/HIV/REVIEWS/HXB2.html" target="_blank">link</a>).
  <br><br>
  <strong>Note:</strong> Since the ground truth sequence in this analysis is <strong>{label}</strong>
  (not HXB2 subtype B), the gene boundaries shown are <strong>approximate</strong>.
  Insertions and deletions between subtypes can shift exact gene start/stop positions by several
  nucleotides. The positions below are accurate for HXB2 and serve as a reliable approximate
  reference for other HIV-1 Group M subtypes.

  <table class="gene-table">
    <thead>
      <tr><th>Region</th><th>HXB2 Start</th><th>HXB2 End</th><th>Source</th></tr>
    </thead>
    <tbody>
      <tr><td>5&prime; LTR</td><td>1</td><td>634</td><td>K03455 repeat_region</td></tr>
      <tr><td>gag</td><td>790</td><td>2292</td><td>K03455 CDS</td></tr>
      <tr><td>pol</td><td>2085</td><td>5096</td><td>K03455 CDS</td></tr>
      <tr><td>vif</td><td>5041</td><td>5619</td><td>K03455 CDS</td></tr>
      <tr><td>vpr</td><td>5559</td><td>5850</td><td>K03455 CDS</td></tr>
      <tr><td>tat (exon 1)</td><td>5831</td><td>6045</td><td>K03455 CDS/exon</td></tr>
      <tr><td>rev (exon 1)</td><td>5970</td><td>6045</td><td>K03455 CDS/exon</td></tr>
      <tr><td>vpu</td><td>6062</td><td>6310</td><td>K03455 CDS</td></tr>
      <tr><td>env</td><td>6225</td><td>8795</td><td>K03455 CDS</td></tr>
      <tr><td>nef</td><td>8797</td><td>9417</td><td>K03455 CDS</td></tr>
      <tr><td>3&prime; LTR</td><td>9086</td><td>9719</td><td>K03455 repeat_region</td></tr>
    </tbody>
  </table>
</div>

<!-- Mutation type breakdown -->
<div class="mutation-grid" id="mutationGrid"></div>

<!-- Mismatch table -->
<div class="table-container">
  <div class="chart-title">Mismatch Details</div>
  <div class="filter-bar">
    <button class="filter-btn active" onclick="filterTable('all')">All</button>
    <button class="filter-btn" onclick="filterTable('shared_same_error')">Shared (same)</button>
    <button class="filter-btn" onclick="filterTable('shared_diff_error')">Shared (diff)</button>
    <button class="filter-btn" onclick="filterTable('v5_only')">v5-only</button>
    <button class="filter-btn" onclick="filterTable('shiver_only')">Shiver-only</button>
  </div>
  <table>
    <thead>
      <tr>
        <th>Genome Pos</th>
        <th>GT</th>
        <th>v5 Call</th>
        <th>v5 Type</th>
        <th>Shiver Call</th>
        <th>Shiver Type</th>
        <th>Category</th>
        <th>HIV Region (HXB2)</th>
      </tr>
    </thead>
    <tbody id="tableBody"></tbody>
  </table>
</div>

<script>
const genomeLen = {gt_length};
const regions = {regions_js};

const sharedSame = {to_js_array(shared_same)};
const sharedDiff = {to_js_array(shared_diff)};
const v5Only = {to_js_array(v5_only)};
const shOnly = {to_js_array(sh_only)};

// Full data for table
const allData = {generate_table_json(results, hiv_regions)};

// ---- Genome map ----
function drawGenomeMap() {{
  const canvas = document.getElementById('genomeMap');
  const rect = canvas.parentElement.getBoundingClientRect();
  canvas.width = rect.width - 40;
  canvas.height = 120;
  const ctx = canvas.getContext('2d');
  const W = canvas.width, H = canvas.height;
  const margin = {{left: 50, right: 20, top: 10, bottom: 35}};
  const plotW = W - margin.left - margin.right;
  const plotH = H - margin.top - margin.bottom;

  const x = pos => margin.left + (pos / genomeLen) * plotW;

  // Background
  ctx.fillStyle = '#1a1a2e';
  ctx.fillRect(0, 0, W, H);

  // Gene regions
  const regionH = 14;
  const regionY = margin.top + plotH - regionH - 2;
  const colors = ['#2a2a4a','#1f2f4a','#2a2a4a','#1f2f4a','#2a2a4a','#1f2f4a','#2a2a4a','#1f2f4a','#2a2a4a','#1f2f4a','#2a2a4a'];
  regions.forEach((r, i) => {{
    const x1 = x(r.start), x2 = x(r.end);
    ctx.fillStyle = colors[i % colors.length];
    ctx.fillRect(x1, regionY, x2-x1, regionH);
    ctx.fillStyle = '#888';
    ctx.font = '10px sans-serif';
    ctx.textAlign = 'center';
    ctx.fillText(r.name, (x1+x2)/2, regionY + regionH + 12);
  }});

  // Genome line
  ctx.strokeStyle = '#555';
  ctx.lineWidth = 1;
  const lineY = margin.top + 10;
  ctx.beginPath();
  ctx.moveTo(margin.left, lineY);
  ctx.lineTo(margin.left + plotW, lineY);
  ctx.stroke();

  // Ticks
  ctx.fillStyle = '#888';
  ctx.font = '10px sans-serif';
  ctx.textAlign = 'center';
  for (let p = 0; p <= genomeLen; p += 1000) {{
    const xp = x(p);
    ctx.beginPath();
    ctx.moveTo(xp, lineY);
    ctx.lineTo(xp, lineY + 4);
    ctx.stroke();
    ctx.fillText(p.toString(), xp, lineY + 14);
  }}

  // Draw mismatches as vertical lines
  function drawMarkers(data, color, yOff) {{
    ctx.strokeStyle = color;
    ctx.lineWidth = 1.5;
    data.forEach(d => {{
      const xp = x(d.pos);
      ctx.beginPath();
      ctx.moveTo(xp, lineY + 18 + yOff);
      ctx.lineTo(xp, lineY + 30 + yOff);
      ctx.stroke();
    }});
  }}

  // Labels for rows
  ctx.fillStyle = '#ce93d8'; ctx.textAlign = 'right'; ctx.font = '10px sans-serif';
  ctx.fillText('Shared', margin.left - 5, lineY + 27);
  ctx.fillStyle = '#4fc3f7';
  ctx.fillText('v5', margin.left - 5, lineY + 42);
  ctx.fillStyle = '#ff8a65';
  ctx.fillText('Shiver', margin.left - 5, lineY + 57);

  drawMarkers([...sharedSame, ...sharedDiff], '#ce93d8', 0);
  drawMarkers(v5Only, '#4fc3f7', 15);
  drawMarkers(shOnly, '#ff8a65', 30);
}}

drawGenomeMap();
window.addEventListener('resize', drawGenomeMap);

// ---- Mutation breakdown bars ----
function drawMutationBreakdown() {{
  const grid = document.getElementById('mutationGrid');
  
  // Count mutations by type for each method
  const v5Muts = {{}};
  const shMuts = {{}};
  allData.forEach(d => {{
    if (d.v5_mut_label !== '-' && d.v5_mut_label) {{
      v5Muts[d.v5_mut_label] = (v5Muts[d.v5_mut_label] || 0) + 1;
    }}
    if (d.sh_mut_label !== '-' && d.sh_mut_label) {{
      shMuts[d.sh_mut_label] = (shMuts[d.sh_mut_label] || 0) + 1;
    }}
  }});

  function makeBars(muts, title, color) {{
    const sorted = Object.entries(muts).sort((a,b) => b[1]-a[1]);
    const max = sorted.length > 0 ? sorted[0][1] : 1;
    let html = '<div class="mutation-box"><h3 style="color:' + color + '">' + title + '</h3>';
    sorted.forEach(([label, count]) => {{
      const w = Math.max(2, (count/max)*100);
      html += '<div class="mut-bar">';
      html += '<span class="mut-bar-label">' + label + '</span>';
      html += '<div class="mut-bar-fill" style="width:' + w + '%;background:' + color + '"></div>';
      html += '<span>' + count + '</span></div>';
    }});
    html += '</div>';
    return html;
  }}

  grid.innerHTML = makeBars(v5Muts, 'v5 Mutations', '#4fc3f7') + makeBars(shMuts, 'Shiver Mutations', '#ff8a65');
}}
drawMutationBreakdown();

// ---- Table ----
function getRegion(pos) {{
  for (const r of regions) {{
    if (pos >= r.start && pos <= r.end) return r.name;
  }}
  return '-';
}}

function renderTable(filter) {{
  const tbody = document.getElementById('tableBody');
  const filtered = filter === 'all' ? allData : allData.filter(d => d.category === filter);
  let html = '';
  filtered.forEach(d => {{
    const catClass = 'cat-' + d.category;
    const v5Class = d.v5_mut_type === 'transition' ? 'mut-transition' : d.v5_mut_type === 'transversion' ? 'mut-transversion' : '';
    const shClass = d.sh_mut_type === 'transition' ? 'mut-transition' : d.sh_mut_type === 'transversion' ? 'mut-transversion' : '';
    html += '<tr>';
    html += '<td>' + d.genome_pos + '</td>';
    html += '<td>' + d.gt_base + '</td>';
    html += '<td style="color:' + (d.v5_correct ? '#81c784' : '#e57373') + '">' + d.v5_base + '</td>';
    html += '<td class="' + v5Class + '">' + (d.v5_mut_label || '-') + '</td>';
    html += '<td style="color:' + (d.sh_correct ? '#81c784' : '#e57373') + '">' + d.sh_base + '</td>';
    html += '<td class="' + shClass + '">' + (d.sh_mut_label || '-') + '</td>';
    html += '<td class="' + catClass + '">' + d.category.replace(/_/g, ' ') + '</td>';
    html += '<td>' + getRegion(d.genome_pos) + '</td>';
    html += '</tr>';
  }});
  tbody.innerHTML = html;
}}

function filterTable(filter) {{
  document.querySelectorAll('.filter-btn').forEach(b => b.classList.remove('active'));
  event.target.classList.add('active');
  renderTable(filter);
}}

renderTable('all');
</script>
</body>
</html>"""
    
    with open(filepath, 'w') as f:
        f.write(html)


def generate_table_json(results, hiv_regions):
    """
    Build the JS array of mismatch objects embedded in the HTML table.

    Despite the name this isn't real JSON (unquoted keys, single-quoted
    strings); it's a JS object literal for inline embedding only. Don't
    try to JSON.parse it, i just kept the function name through the evolution
    of the script (worked with JSON in some tries).
    """
    items = []
    for r in results:
        items.append(
            '{' +
            f"genome_pos:{r['genome_pos']},"
            f"gt_base:'{r['gt_base']}',"
            f"v5_base:'{r['v5_base']}',"
            f"v5_correct:{'true' if r['v5_correct'] else 'false'},"
            f"v5_mut_type:'{r['v5_mut_type']}',"
            f"v5_mut_label:'{r['v5_mut_label']}',"
            f"sh_base:'{r['shiver_base']}',"
            f"sh_correct:{'true' if r['shiver_correct'] else 'false'},"
            f"sh_mut_type:'{r['shiver_mut_type']}',"
            f"sh_mut_label:'{r['shiver_mut_label']}',"
            f"category:'{r['category']}'"
            + '}'
        )
    return '[' + ','.join(items) + ']'



# MAIN


def find_sequence_by_hint(seqs, hints):
    """Find a sequence name containing any of the hint strings."""
    for name in seqs:
        for hint in hints:
            if hint.lower() in name.lower():
                return name
    return None


def main():
    parser = argparse.ArgumentParser(description='Mismatch analysis: v5 vs Shiver vs Ground Truth')
    parser.add_argument('--v5', required=True, help='v5 consensus FASTA')
    parser.add_argument('--shiver', required=True, help='Shiver consensus FASTA')
    parser.add_argument('--gt', required=True, help='Ground truth FASTA')
    parser.add_argument('--label', required=True, help='Dataset label (e.g., CON_A1)')
    parser.add_argument('--output-dir', default='.', help='Output directory')
    parser.add_argument('--aligned', default=None,
                        help='Pre-aligned FASTA with all 3 sequences (skip MAFFT)')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    if args.aligned:
        print(f"Loading pre-aligned file: {args.aligned}")
        aligned = load_fasta(args.aligned)
    else:
        v5_seqs = load_fasta(args.v5)
        sh_seqs = load_fasta(args.shiver)
        gt_seqs = load_fasta(args.gt)

        combined = {}
        v5_name = list(v5_seqs.keys())[0]
        sh_name = list(sh_seqs.keys())[0]
        gt_name = list(gt_seqs.keys())[0]

        # Rename to avoid conflicts
        combined['V5_' + v5_name] = v5_seqs[v5_name]
        combined['SHIVER_' + sh_name] = sh_seqs[sh_name]
        combined['GT_' + gt_name] = gt_seqs[gt_name]

        aln_path = os.path.join(args.output_dir, f'{args.label}_aligned.fasta')
        print(f"Aligning sequences with MAFFT...")
        aligned = align_sequences(combined, aln_path)

    # Identify which sequence is which
    gt_name = find_sequence_by_hint(aligned, ['GT_', 'EMBOSS', 'ground_truth', 'simref', 'consensus_fixed'])
    v5_name = find_sequence_by_hint(aligned, ['V5_', 'v5', 'node_v5', 'node_based_v5'])
    sh_name = find_sequence_by_hint(aligned, ['SHIVER_', 'shiver', 'Shiver'])

    if not gt_name or not v5_name or not sh_name:
        print(f"\nAvailable sequence names: {list(aligned.keys())}")
        print(f"  Detected GT:     {gt_name}")
        print(f"  Detected v5:     {v5_name}")
        print(f"  Detected Shiver: {sh_name}")
        if not gt_name:
            print("ERROR: Could not identify ground truth sequence.")
        if not v5_name:
            print("ERROR: Could not identify v5 consensus sequence.")
        if not sh_name:
            print("ERROR: Could not identify Shiver consensus sequence.")
        print("\nPlease check sequence names in your aligned FASTA.")
        sys.exit(1)

    print(f"\n  GT:     {gt_name}")
    print(f"  v5:     {v5_name}")
    print(f"  Shiver: {sh_name}")

    gt_seq = aligned[gt_name]
    v5_seq = aligned[v5_name]
    sh_seq = aligned[sh_name]

    # ground truth ungapped length
    gt_length = sum(1 for b in gt_seq if b in 'ATGCN')

    results = analyze_mismatches(gt_seq, v5_seq, sh_seq)

    tsv_path = os.path.join(args.output_dir, f'{args.label}_mismatches.tsv')
    write_tsv(results, tsv_path, args.label)
    print(f"\n  TSV saved: {tsv_path}")

    summary = print_summary(results, args.label, gt_length)

    html_path = os.path.join(args.output_dir, f'{args.label}_mismatch_viz.html')
    generate_html_visualization(results, args.label, gt_length, html_path)
    print(f"  HTML saved: {html_path}")

    return summary


if __name__ == '__main__':
    main()