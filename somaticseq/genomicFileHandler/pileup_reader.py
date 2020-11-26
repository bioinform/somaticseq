#!/usr/bin/env python3

# This class will be deprecated

import sys, os, gzip, re, math

nan = float('nan')
inf = float('inf')

def seq(reads):

    '''A function to convert a string into a character generator.'''
    for read_i in reads:
        yield read_i


class Pileup_line:

    def __init__(self, pileup_line):

        '''Argument is a line in pileup file.'''
        self.pileup_line = pileup_line.rstrip('\n')

        try:
            self.chromosome, self.position, self.refbase, self.dp, self.reads, self.qualities = pileup_line.rstrip('\n').split('\t')
            self.position = int(self.position)
            self.dp = int(self.dp)

        except ValueError:

            try:
                self.chromosome, self.position, self.refbase, self.dp = pileup_line.rstrip('\n').split('\t')
                self.position = int(self.position)
                self.dp = int(self.dp)
                self.reads = self.qualities = ''

            except ValueError:
                self.chromosome = self.refbase = self.reads = self.qualities = ''
                self.position = None
                self.dp = 0


    def alt_read_count(self, pattern):
        alt_reads = re.findall(pattern, self.reads, re.I)
        total_count = len(alt_reads)
        return total_count

    def base_reads(self):

        # Make the base_calls into a generator:
        base_calls = seq(self.reads)

        ref_forward_count = ref_reverse_count = 0
        alt_forward = []
        alt_reverse = []
        del_forward = []
        del_reverse = []
        ins_forward = []
        ins_reverse = []
        n_count = N_count = 0

        for base_i in base_calls:

            # Indicate beginning of a read:
            if base_i == '^':

                # Skip the MQ ascii:
                next(base_calls)

            # Reference forward reads
            elif base_i == '.':
                ref_forward_count = ref_forward_count + 1

            # Reference reverse reads
            elif base_i == ',':
                ref_reverse_count = ref_reverse_count + 1

            # Random alt base
            elif base_i == 'n':
                n_count = n_count + 1
            elif base_i == 'N':
                N_count = N_count + 1

            # SNP
            elif re.match(r'[GCTAU]', base_i):
                alt_forward.append(base_i)
            elif re.match(r'[gctau]', base_i):
                alt_reverse.append(base_i)

            # Deletion
            elif base_i == '-':

                num = ''
                base_i = next(base_calls)

                while base_i.isnumeric():
                    num = num + base_i
                    base_i = next(base_calls)

                num = int(num)
                deleted_seq = base_i
                for i in range(num-1):
                    deleted_seq = deleted_seq + next(base_calls)

                if deleted_seq.isupper():
                    del_forward.append(deleted_seq)
                elif deleted_seq.islower():
                    del_reverse.append(deleted_seq)

            # Insertion
            elif base_i == '+':
                num = ''
                base_i = next(base_calls)

                while base_i.isnumeric():
                    num = num + base_i
                    base_i = next(base_calls)

                num = int(num)
                inserted_seq = base_i
                for i in range(num-1):
                    inserted_seq = inserted_seq + next(base_calls)

                if inserted_seq.isupper():
                    ins_forward.append(inserted_seq)
                elif inserted_seq.islower():
                    ins_reverse.append(inserted_seq)

        return ref_forward_count, ref_reverse_count, alt_forward, alt_reverse, del_forward, del_reverse, ins_forward, ins_reverse, n_count, N_count


    def total_insertion_calls(self):
        ins = re.findall(r'\+[0-9]+[ACGTNacgtn]+', self.reads)
        total_count = len(ins)
        return total_count

    def total_deletion_calls(self):
        dels = re.findall(r'-[0-9]+[ACGTNacgtn]+', self.reads)
        total_count = len(dels)
        return total_count

    def indel_fraction(self):
        ins = self.total_insertion_calls()
        dels = self.total_deletion_calls()

        try:
            fraction = (ins + dels ) / self.dp
        except ZeroDivisionError:
            fraction = 0
        return fraction





class Base_calls(Pileup_line):

    def __init__(self, pileup_line):

        Pileup_line.__init__(self, pileup_line)

        # Make the base_calls into a generator:
        base_calls = seq(self.reads)

        ref_forward_count = ref_reverse_count = 0
        n_count = a_count = c_count = g_count = t_count = N_count = A_count = C_count = G_count = T_count = 0
        del_forward = []
        del_reverse = []
        ins_forward = []
        ins_reverse = []

        for base_i in base_calls:

            # Indicate beginning of a read:
            if base_i == '^':

                # Skip the MQ ascii:
                next(base_calls)

            # Reference forward reads
            elif base_i == '.':
                ref_forward_count = ref_forward_count + 1

            # Reference reverse reads
            elif base_i == ',':
                ref_reverse_count = ref_reverse_count + 1

            # Random alt base
            elif base_i == 'n':
                n_count = n_count + 1
            elif base_i == 'N':
                N_count = N_count + 1

            # SNP
            elif re.match(r'a', base_i):
                a_count += 1
            elif re.match(r'c', base_i):
                c_count += 1
            elif re.match(r'g', base_i):
                g_count += 1
            elif re.match(r't', base_i):
                t_count += 1
            elif re.match(r'A', base_i):
                A_count += 1
            elif re.match(r'C', base_i):
                C_count += 1
            elif re.match(r'G', base_i):
                G_count += 1
            elif re.match(r'T', base_i):
                T_count += 1

            # Deletion
            elif base_i == '-':

                num = ''
                base_i = next(base_calls)

                while base_i.isnumeric():
                    num = num + base_i
                    base_i = next(base_calls)

                num = int(num)
                deleted_seq = base_i
                for i in range(num-1):
                    deleted_seq = deleted_seq + next(base_calls)

                if deleted_seq.isupper():
                    del_forward.append(deleted_seq)
                elif deleted_seq.islower():
                    del_reverse.append(deleted_seq)

            # Insertion
            elif base_i == '+':
                num = ''
                base_i = next(base_calls)

                while base_i.isnumeric():
                    num = num + base_i
                    base_i = next(base_calls)

                num = int(num)
                inserted_seq = base_i
                for i in range(num-1):
                    inserted_seq = inserted_seq + next(base_calls)

                if inserted_seq.isupper():
                    ins_forward.append(inserted_seq)
                elif inserted_seq.islower():
                    ins_reverse.append(inserted_seq)

        # Before going there, find out what is the ref base and re-assign accordingly:
        if self.refbase.upper() == 'A':
            a_count = ref_reverse_count
            A_count = ref_forward_count
        elif self.refbase.upper() == 'C':
            c_count = ref_reverse_count
            C_count = ref_forward_count
        elif self.refbase.upper() == 'G':
            g_count = ref_reverse_count
            G_count = ref_forward_count
        elif self.refbase.upper() == 'T':
            t_count = ref_reverse_count
            T_count = ref_forward_count
        elif self.refbase.upper() == 'N':
            n_count = ref_reverse_count
            N_count = ref_forward_count

        deletion_calls = {}
        for call_i in del_forward:
            try:
                deletion_calls[call_i] += 1
            except KeyError:
                deletion_calls[call_i] = 1

        for call_i in del_reverse:
            try:
                deletion_calls[call_i.upper()] += 1
            except KeyError:
                deletion_calls[call_i.upper()] = 1

        insertion_calls = {}
        for call_i in ins_forward:
            try:
                insertion_calls[call_i] += 1
            except KeyError:
                insertion_calls[call_i] = 1

        for call_i in ins_reverse:
            try:
                insertion_calls[call_i.upper()] += 1
            except KeyError:
                insertion_calls[call_i.upper()] = 1

        # Returns these:
        self.A   = A_count, a_count
        self.C   = C_count, c_count
        self.G   = G_count, g_count
        self.T   = T_count, t_count
        self.DEL = len(del_forward), len(del_reverse)
        self.INS = len(ins_forward), len(ins_reverse)
        self.N   = N_count, n_count
        self.ref = ref_forward_count, ref_reverse_count

        self.deletions       = del_forward, del_reverse  # list of deletions
        self.insertions      = ins_forward, ins_reverse  # list of insertions
        self.deletion_calls  = deletion_calls  # dictionary of deletion calls (strand agnostic) and their occurrence
        self.insertion_calls = insertion_calls # dictionary of insertion calls (strand agnostic) and their occurrence

