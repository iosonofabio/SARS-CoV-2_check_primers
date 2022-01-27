# vim: fdm=indent
'''
author:     Fabio Zanini
date:       27/01/22
content:    Run in-silico PCR for SARS-CoV-2
'''
import os
import sys
from Bio import SeqIO
from Bio.Seq import reverse_complement as rc
sys.path.append('../libraries')
import seqanpy


if __name__ == '__main__':

    refs = list(SeqIO.parse('../data/references/sars.fasta', 'fasta'))

    primers = {
        '1': {'FWD': 'TGCTATCTCTGGGACCAATG',
              'REV': 'TGGAAGCAAAATAAACACCATC'},
        '2': {'FWD': 'AGAACTCAATTACCCCCTGCATAC',
              'REV': 'GGAAAAGAAAGGTAAGAACAAGTCC'},
    }

    products_exp = {
        '1': 'TGCTATCTCTGGGACCAATGGTACTAAGAGGTTTGATAACCCTGTCCTACCATTTAATGATGGTGTTTATTTTGCTTCCA',
        '2': 'AGAACTCAATTACCCCCTGCATACACTAATTCTTTCACACGTGGTGTTTATTACCCTGACAAAGTTTTCAGATCCTCAGTTTTACATTCAACTCAGGACTTGTTCTTACCTTTCTTTTCC'
    }


    def find_amplicons(primer_pair):
        prfwd, prrev = primer_pair['FWD'], primer_pair['REV']
        prrev_rc = rc(prrev)

        amplicons = {}
        for ref in refs:
            seq = ref.seq
            scoref, alifr, alifp = seqanpy.align_overlap(seq, prfwd)
            scorer, alirr, alirp = seqanpy.align_overlap(seq, prrev_rc)
            start_fwd = len(alifp) - len(alifp.lstrip('-'))
            end_fwd = len(alifp.rstrip('-'))
            end_rev = len(alirp.rstrip('-'))

            print(ref.name, 'fwd', scoref, 'of', 3 * len(prfwd))
            print(ref.name, 'rev', scorer, 'of', 3 * len(prrev))

            if scoref < 3 * (len(prfwd) - 1):
                continue
            if scorer < 3 * (len(prrev) - 1):
                continue

            # Borderline case: check that the 3' end is perfect
            threeplen = 5
            if alifr[end_fwd - threeplen: end_fwd] != alifp[end_fwd - threeplen: end_fwd]:
                continue
            if alirr[end_rev - threeplen: end_rev] != alirp[end_rev - threeplen: end_rev]:
                continue

            if end_rev <= start_fwd:
                continue

            amp = str(seq[start_fwd: end_rev])
            amplicons[ref.name] = amp
        return amplicons

    amplicons = {key: find_amplicons(value) for key, value in primers.items()}
