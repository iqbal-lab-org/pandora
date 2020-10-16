/*
 * Adapted from https://github.com/iqbal-lab/cortex/blob/master/include/basic/binary_kmer.h
 *
 * Copyright 2009-2011 Zamin Iqbal and Mario Caccamo
 *
 * CORTEX project contacts:
 * 		M. Caccamo (mario.caccamo@bbsrc.ac.uk) and
 * 		Z. Iqbal (zam@well.ox.ac.uk)
 *
 * **********************************************************************
 *
 * This file is part of CORTEX.
 *
 * CORTEX is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CORTEX is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CORTEX.  If not, see <http://www.gnu.org/licenses/>.
 *
 * **********************************************************************
 */
/*
  binary_kmer.h - basic manipulation of binary kmers
*/

#ifndef BINARY_KMER_H_
#define BINARY_KMER_H_

#include <inttypes.h>
#include "cortex/ns.h"

#include "global.h"
#include "event_encoding.h"

#define BINVERSION 6

// Think of NUMBER_OF_BITFIELDS_IN_BINARY_KMER as the number of uint64_t we
// encode the kmer in
constexpr int NB_BITFIELDS_IN_BINARY_KMER {1};
constexpr int BINARY_KMER_BYTES {(NB_BITFIELDS_IN_BINARY_KMER * sizeof(uint64_t))};
using bitfield_bits = uint64_t;


class cortex::BinaryKmer {
private:
    bitfield_bits m_bitfield[NB_BITFIELDS_IN_BINARY_KMER];

public:
    const short kmer_size;
    // initialise to zero and kmer size as max kmer size
    BinaryKmer();
    // initialiase to zero
    BinaryKmer(short kmer_size);
    // set all bitfields
    void set_all(bitfield_bits val);

    friend bool operator== (const BinaryKmer& left, const BinaryKmer& right);
    friend bool operator!= (const BinaryKmer& left, const BinaryKmer& right);
    friend bool operator<= (const BinaryKmer& left, const BinaryKmer& right);
    // right shift one base
    void right_shift();
    // left shift one base
    void left_shift();
    // left shirt one base and insert new base at right end
    void left_shift_add_right(Nucleotide base);
};


typedef struct
{
    int nkmers;
    BinaryKmer * kmer;
} KmerSlidingWindow;


//a set of KmerArrays
typedef struct
{
    int nwindows;
    int max_nwindows;
    KmerSlidingWindow * window;
} KmerSlidingWindowSet;


// basic BinaryKmer operations
//void    binary_kmer_initialise_to_zero(BinaryKmer* bkmer);
//void    binary_kmer_assignment_operator(BinaryKmer left, BinaryKmer right);
//void    binary_kmer_set_all_bitfields(BinaryKmer assignee, bitfield_of_64bits val);
//boolean binary_kmer_comparison_operator(const BinaryKmer const left, const BinaryKmer const right);
//boolean binary_kmer_less_than(const BinaryKmer const left, const BinaryKmer const right,short kmer_size);
//void    binary_kmer_right_shift_one_base(BinaryKmer kmer);
//void    binary_kmer_left_shift_one_base(BinaryKmer kmer, short kmer_size);
//void    binary_kmer_left_shift_one_base_and_insert_new_base_at_right_end(BinaryKmer* bkmer, Nucleotide n, short kmer_size);


char * nucleotides_to_string(Nucleotide * nucleotides, int length, char * string);

//get overlapping kmers from sequence
//int get_sliding_windows_from_sequence(char * sequence,char * qualities, int length, char quality_cutoff, short kmer_size, KmerSlidingWindowSet * windows, int max_windows, int max_kmers);
int get_sliding_windows_from_sequence(char * seq,  char * qualities, int length, char quality_cut_off, short kmer_size, KmerSlidingWindowSet * windows,
    int max_windows, int max_kmers, boolean break_homopolymers, int homopolymer_cutoff);


char * binary_kmer_to_seq(BinaryKmer* kmer, short kmer_size, char * seq);

BinaryKmer* seq_to_binary_kmer(char * seq, short kmer_size, BinaryKmer* prealloced_kmer);

char * seq_reverse_complement(char * in, int length, char * out);

BinaryKmer* binary_kmer_reverse_complement(BinaryKmer* kmer, short kmer_size, BinaryKmer* prealloc_reverse_kmer);

Nucleotide binary_kmer_get_first_nucleotide(BinaryKmer* kmer, short kmer_size);
Nucleotide binary_kmer_get_last_nucleotide(BinaryKmer* kmer);

char reverse_char_nucleotide(char c);

void binary_kmer_alloc_kmers_set(KmerSlidingWindowSet * windows, int max_windows, int max_kmers);
void binary_kmer_free_kmers(KmerSlidingWindow * *);
void binary_kmer_free_kmers_set(KmerSlidingWindowSet * *);

void nucleotide_iterator(void (*f)(Nucleotide));
void nucleotide_iterator_orientation(void (*f)(Nucleotide n,Orientation o));

#endif /* BINARY_KMER_H_ */
