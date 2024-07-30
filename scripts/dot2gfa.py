#!/usr/bin/env python3

import click
import bisect
import re
from Bio import SeqIO


def remove_overlapping_prefix(vertexLabels, edges, k):
    """Remove overlapping prefix from selected vertices.

    Algorithm: if a vertex has in-degree > 0, then
               remove first (k-1) characters from its label
    """
    trimVertexLabel = [0] * len(vertexLabels)
    for (u, v) in edges:
        assert u >= 0
        assert u < len(vertexLabels)
        assert v >= 0
        assert v < len(vertexLabels)
        trimVertexLabel[v] = 1

    for i in range(len(vertexLabels)):
        if trimVertexLabel[i] == 1:
            vertexLabels[i] = vertexLabels[i][k - 1:]


def get_out_neighbours(edges):
    """Compute out-neighbours of each vertex."""
    out_index = dict()
    for (u, v) in edges:
        out_index[u] = list()

    for (u, v) in edges:
        out_index[u].append(v)

    return out_index


@click.command()
@click.option('-k', type=int, default=25, show_default=True, help='k-mer size')
@click.option('--fasta', '-f', type=click.Path(exists=True), help='FASTA file')
@click.option('--output', '-o', type=click.Path(), default="output.gfa",
              show_default=True, help='Output GFA file')
@click.argument('input', type=click.Path(exists=True))
def dot2gfa(input, k, fasta, output):
    """Convert output of splitMEM in DOT format to GFA1 format."""
    lpos = list()
    fasta_pos_file = input + "fastaPos.txt"
    with open(fasta_pos_file, 'r') as f:
        lpos = [int(l.strip()) for l in f]

    lrec = [str(r.seq) for r in SeqIO.parse(fasta, "fasta")]

    vertexLabels = list()
    edges = list()

    with open(input, 'r') as f:
        currentVertexId = -1
        for line in f:
            tokens = line.strip().split()
            if len(tokens) == 2:  # must be a vertex label
                currentVertexId += 1
                assert int(tokens[0]) == currentVertexId
                # remove leading '[label="' and trailing '"]", and
                # choose the last 'pos:len'
                label_range = tokens[1][8:-2].split(',')[-1].split(':')
                pos = int(label_range[0])
                length = int(label_range[1])
                idx = bisect.bisect_left(lpos, position + 1) - 1
                offset = position - lpos[idx]
                # NOTE: replace N with A, as splitMEM does
                seq = re.sub(r'[N]', 'A', lrec[idx])
                # NOTE: splitMEM add a "N" at the end of each FASTA entry
                seq += "N"
                subseq = seq[offset:offset + length]
                subseq = re.sub(r'[^ACGTN\$]', '', subseq)  #get rid of ambiguous character
                subseq = re.sub(r'[\$]', 'N', subseq)       # get rid of ambiguous characters
                vertexLabels.append(subseq)
            elif len(tokens) == 3:
                if tokens[1] == "->":  # must be an edge
                    assert int(tokens[0]) == currentVertexId
                    edges.append((int(tokens[0]), int(tokens[2])))

    # Remove duplicates
    edges = list(dict.fromkeys(edges))

    remove_overlapping_prefix(vertexLabels, edges, k)
    out_index = get_out_neighbours(edges)

    with open(output, 'w') as f:
        f.write('H\tVN:Z:1.0\n')
        for u in range(len(vertexLabels)):
            f.write(f'S\t{u+1}\t{vertexLabels[u]}\n')
            if u in out_index:
                for v in out_index[u]:
                    f.write(f'L\t{u+1}\t+\t{v+1}\t+\t0M\n')


if __name__ == '__main__':
    dot2gfa()
