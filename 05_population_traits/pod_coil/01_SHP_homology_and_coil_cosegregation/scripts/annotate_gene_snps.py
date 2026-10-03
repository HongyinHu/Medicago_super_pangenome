#!/usr/bin/env python3
import argparse
import csv


CODON = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}
COMP = str.maketrans('ACGTN', 'TGCAN')


def revcomp(seq):
    return seq.upper().translate(COMP)[::-1]


def read_region_fasta(path):
    header = ''
    seq = []
    with open(path) as handle:
        for line in handle:
            if line.startswith('>'):
                header = line[1:].strip().split()[0]
            else:
                seq.append(line.strip())
    chrom, coords = header.split(':')
    start, end = map(int, coords.split('-'))
    return chrom, start, end, ''.join(seq).upper()


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--association', required=True)
    p.add_argument('--gff', required=True)
    p.add_argument('--gene', required=True)
    p.add_argument('--region-fasta', required=True)
    p.add_argument('--output', required=True)
    args = p.parse_args()

    gene_row = None
    cds = []
    with open(args.gff) as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            f = line.rstrip().split('\t')
            if len(f) != 9:
                continue
            attrs = f[8]
            if f[2] == 'gene' and ('ID=' + args.gene) in attrs:
                gene_row = f
            if f[2] == 'CDS' and ('Parent=' + args.gene + '.1') in attrs:
                cds.append((int(f[3]), int(f[4])))
    if gene_row is None or not cds:
        raise SystemExit('Gene/CDS not found: ' + args.gene)
    chrom = gene_row[0]
    gene_start, gene_end, strand = int(gene_row[3]), int(gene_row[4]), gene_row[6]
    rchrom, rstart, rend, region_seq = read_region_fasta(args.region_fasta)
    if chrom != rchrom:
        raise SystemExit('Chromosome mismatch')

    ordered = sorted(cds, reverse=(strand == '-'))
    full_cds = ''
    segment_offsets = []
    for start, end in ordered:
        fragment = region_seq[start - rstart:end - rstart + 1]
        if strand == '-':
            fragment = revcomp(fragment)
        segment_offsets.append((start, end, len(full_cds)))
        full_cds += fragment

    with open(args.association) as inp, open(args.output, 'w', newline='') as out:
        reader = csv.DictReader(inp, delimiter='\t')
        fields = reader.fieldnames + ['feature', 'cds_nt_position', 'aa_position', 'ref_codon', 'alt_codon', 'aa_change']
        writer = csv.DictWriter(out, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        for row in reader:
            pos = int(row['pos'])
            ref, alt = row['ref'].upper(), row['alt'].upper()
            feature = 'intergenic'
            cds_offset = None
            if gene_start <= pos <= gene_end:
                feature = 'intron'
                for start, end, before in segment_offsets:
                    if start <= pos <= end:
                        feature = 'CDS'
                        cds_offset = before + ((pos - start) if strand == '+' else (end - pos))
                        break
            elif pos < gene_start:
                feature = 'upstream' if strand == '+' else 'downstream'
            elif pos > gene_end:
                feature = 'downstream' if strand == '+' else 'upstream'

            row.update({'feature': feature, 'cds_nt_position': '', 'aa_position': '', 'ref_codon': '', 'alt_codon': '', 'aa_change': ''})
            if cds_offset is not None and len(ref) == 1 and len(alt) == 1:
                transcript_ref = ref if strand == '+' else revcomp(ref)
                transcript_alt = alt if strand == '+' else revcomp(alt)
                codon_start = (cds_offset // 3) * 3
                ref_codon = full_cds[codon_start:codon_start + 3]
                alt_chars = list(ref_codon)
                alt_chars[cds_offset % 3] = transcript_alt
                alt_codon = ''.join(alt_chars)
                ref_aa = CODON.get(ref_codon, 'X')
                alt_aa = CODON.get(alt_codon, 'X')
                row.update({
                    'cds_nt_position': cds_offset + 1,
                    'aa_position': cds_offset // 3 + 1,
                    'ref_codon': ref_codon,
                    'alt_codon': alt_codon,
                    'aa_change': '{}{}{}'.format(ref_aa, cds_offset // 3 + 1, alt_aa)
                })
            writer.writerow(row)


if __name__ == '__main__':
    main()
