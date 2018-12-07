// Copyright (c) 2018, Nicolas Aravena, Diego Arroyuelo.  All rights reserved.
// Use of this source code is governed
// by a BSD license that can be found in the LICENSE file.
/*! \file s9_vector.hpp
   \brief s9_vector.hpp contains the sdsl::s9_vector class, and
          classes which support rank and select for s9_vector.
   \author Nicolás Aravena, Diego Arroyuelo
*/
#ifndef INCLUDED_SDSL_S9_VECTOR
#define INCLUDED_SDSL_S9_VECTOR

#include <sdsl/int_vector.hpp>
#include <sdsl/util.hpp>

using namespace std; 
//! Namespace for the succinct data structure library
namespace sdsl
{

// forward declaration needed for friend declaration
template<uint8_t t_b=1, uint16_t t_bs=128, typename t_int_vector = int_vector<>>
class rank_support_s9;

// forward declaration needed for friend declaration
template<uint8_t t_b=1, uint16_t t_bs=128, typename t_int_vector = int_vector<>>
class select_support_s9;

//! Bit vector gap encoding using S9 with rank and select support.
/*!
 *   TODO: Detailed description
 *   \tparam    t_bs            Block size (default = 128).
 *   \tparam    t_int_vector    Int Vector type for rank and select directories (default = int_vector<>).
 *
 *   References:
 *    - Jorge Dinator, Diego Arroyuelo
 *      ESTRUCTURAS COMPRIMIDAS PARA RANK, SELECT Y NEXT BASADAS EN GAP ENCODING
 *      MEMORIA DE TITULACIÓN 2016.
 *
 *    Block size can be adjusted.
 */
template<uint16_t t_bs=128, typename t_int_vector = int_vector<>>
class s9_vector{
    public:
        typedef typename t_int_vector::size_type    size_type;
        typedef typename t_int_vector::value_type   value_type;
        typedef int_vector<32>::value_type          word_type;

        typedef rank_support_s9<1, t_bs>    rank_1_type;
        typedef rank_support_s9<0, t_bs>    rank_0_type;
        typedef select_support_s9<1, t_bs>  select_1_type;
        //typedef select_support_s9<0, t_bs>  select_0_type;

        friend class rank_support_s9<0, t_bs>;
        friend class rank_support_s9<1, t_bs>;
        friend class select_support_s9<1, t_bs>;
        friend class select_support_s9<0, t_bs>;

    private:
        size_type       m;          //Number of ones in the original bitsequence
        size_type       m_size;     //Length of the original bit vector (n)
        int_vector<32>  m_seq;      //vector with the S9 encoded representation of the ones (seq)
        t_int_vector    m_blocks;   //Array with the starting positions of the S9 blocks (blockS)
        t_int_vector    m_absolute; //Array with the absolute starting positions of the S9 blocks (absolute)
        t_int_vector    m_rank;     //One-level rank directory (rankD)
        
        void copy(const s9_vector& s9){
            m = s9.m;
            m_size = s9.m_size;
            m_seq = s9.m_seq;
            m_blocks = s9.m_blocks;
            m_absolute = s9.m_absolute;
            m_rank = s9.m_rank;
        }

        //The following private members might be moved to a helper class.

        // The following table is used to determine, given a number of integers, 
        // the maximum amount of these integers that can be encoded within a
        // single S9 word.
        uint8_t maxNumbers[29] = {1 ,1 ,2 ,3 ,4 ,5 ,5 ,7 ,7 ,9 ,9 ,9 ,9 ,9 ,14,14,14,14,14,14,14,14,14,14,14,14,14,14,28};
        uint8_t Bits[29]       = {28,28,14,9 ,7 ,5 ,5 ,4 ,4 ,3 ,3 ,3 ,3 ,3 ,2 ,2 ,2 ,2 ,2 ,2 ,2 ,2 ,2 ,2 ,2 ,2 ,2 ,2 ,1 };
        // Shows max size of a integer in a single s9 word for a given size.
        uint8_t T[33] = {1,1,2,3,4,5,7,7,9,9,14,14,14,14,14,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28};
        // Gets the case for encoding a given amount of integers
        uint8_t Case[29] = {0, 1, 2, 3, 4, 5, 0 , 6, 0 , 7, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9};
        // Given a bit count gets the amount of integers that can be encoded 
        uint8_t Numbers[29] = {0, 28, 14, 9, 7, 5, 0 , 4, 0 , 3, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};

        void add_s9_word(word_type &word, int_vector<64> &gaps, size_type i, size_type cant){
            size_type max_bits = Bits[cant]; 

            i += cant - 1;
            word |= gaps[i];
            i--;
            for(uint8_t k = 1; k < cant; k++, i--){
                word <<= max_bits;
                word |= gaps[i];
            }

            word_type header = Case[cant];
            header <<= 28;
            word |= header;
        }

        /** Stores a given bit sequence into 32 bit word
         * @param A 32 bit word
         * @param ini Starting position
         * @param fin Store until end-1
         * @param x Value to be stored
         */
        inline void write_in_word(word_type &A, const size_type ini, const size_type fin, const word_type x) {
            uint8_t W = 32;
            if(ini==fin+1) return;
            uint i=ini/W, j=ini-i*W;
            uint len = (fin-ini+1);
            uint mask = ((j+len) < W ? ~0u << (j+len) : 0)
                | ((W-j) < W ? ~0u >> (W-j) : 0);
            A = (A & mask) | x << j;
        }
        /** Packs numbers in a block
         * @param gaps gaps int vector
         * @param gaps_i gaps starting position
         * @param seq s9 word int vector
         * @param seq_i seq starting position
         * @param block_size block size
         * @return the number of S9 words used to do the packing
         */
        value_type pack_in_block(t_int_vector &gaps, size_type gaps_i, int_vector<32> &seq, size_type seq_i, uint16_t block_size){
            size_type nbits = 0, max_bits = 0, cur_cant = 0;
            size_type words = 0, i = 0;
            while (i < block_size){
                //Bits needed to encode the current integer
                if (gaps[gaps_i + i] == 0) nbits = 1;
                else nbits = sdsl::bits::hi(gaps[gaps_i + i]) + 1; //comparar con libcds
                if (nbits >= max_bits) max_bits = T[nbits];

                // curCant+1 integers can be packed into a single word,
                if (cur_cant + 1 == Bits[max_bits]){
                    //add_s9_word(m_seq[seq_i + words], gaps, gaps_i + i - cur_cant, cur_cant + 1);
                    // using maxBits each 
                    write_in_word(m_seq[seq_i + words], 28, 31, Case[cur_cant+1]);     
                    for (uint k = 0; k < cur_cant + 1; k++)
                        write_in_word(m_seq[seq_i + words], max_bits*k, max_bits*k+(max_bits-1), gaps[gaps_i + i - cur_cant + k]);

                    words++;                   
                    cur_cant = 0;
                    max_bits = 0;
                    i++;
                }
                else 
                    if (cur_cant + 1 > Bits[max_bits]){
                        // adding a new integer exceeds the capacity of a word.
                        // so we try to pack as much of the unpacked integers as we can
                        // into a single word, and ressume the algorithm from there.
                        i -= cur_cant;
                        write_in_word(m_seq[seq_i + words], 28, 31, Case[maxNumbers[cur_cant]]);     
                        max_bits = Bits[cur_cant];
                        for (uint k=0; k < maxNumbers[cur_cant]; k++, i++)
                            write_in_word(m_seq[seq_i + words], max_bits * k, max_bits * k + (max_bits - 1), gaps[gaps_i + i]);
                        words++;
                        cur_cant = 0;
                        max_bits = 0;
                    }
                    else { i++; cur_cant++; }
            }
            // now pack the remaining "cur_cant" integers
            i-= cur_cant;
            while (cur_cant > 0) {
                write_in_word(m_seq[seq_i + words], 28, 31, Case[maxNumbers[cur_cant]]);
                max_bits = Bits[cur_cant];
                for (uint j = 0; j < maxNumbers[cur_cant]; j++, i++)
                    write_in_word(m_seq[seq_i + words], max_bits * j, max_bits * j + (max_bits - 1), gaps[gaps_i + i]);
                words++;
                cur_cant -= maxNumbers[cur_cant];
            }
            return words;
        }

        //! return the sum of a s9 encoded word and assigns r to the amount of words in regAux
        /*! \param regAux 32 bit word.
         *  \param *r Pointer to word counter
         */
        inline word_type sum_s9_word(register word_type regAux, word_type *r) const{
            switch (regAux>>28) {
                    case 1: *r = 1; return (regAux & 0x0fffffff); 
                            
                            break;
                    case 2: *r = 2; return (regAux & 0x00003fff) + ((regAux>>14) & 0x00003fff) + 1;
                            break;
                    case 3: *r = 3; return (regAux & 0x000001ff) + ((regAux>>9) & 0x000001ff) + ((regAux>>18) & 0x000001ff) + 2;
                    break;
                    case 4: *r = 4; return (regAux & 0x0000007f) + ((regAux>>7) & 0x0000007f) + ((regAux>>14) & 0x0000007f) 
                                    + ((regAux>>21) & 0x0000007f) + 3;
                    break;
                    case 5: *r = 5; return (regAux & 0x0000001f) + ((regAux>>5) & 0x0000001f) + ((regAux>>10) & 0x0000001f) 
                                    + ((regAux>>15) & 0x0000001f) + ((regAux>>20) & 0x0000001f) + 4;
                    break;
                    case 6: *r = 7; return (regAux & 0x0000000f) + ((regAux>>4) & 0x0000000f) + ((regAux>>8) & 0x0000000f) 
                                    + ((regAux>>12) & 0x0000000f) + ((regAux>>16) & 0x0000000f) + ((regAux>>20) & 0x0000000f) 
                                    + ((regAux>>24) & 0x0000000f) + 6;
                    break;
                    case 7: *r = 9; return (regAux & 0x00000007) + ((regAux>>3) & 0x00000007) + ((regAux>>6) & 0x00000007) 
                                    + ((regAux>>9) & 0x00000007) + ((regAux>>12) & 0x00000007) + ((regAux>>15) & 0x00000007) 
                                    + ((regAux>>18) & 0x00000007) + ((regAux>>21) & 0x00000007) + ((regAux>>24) & 0x00000007) + 8;
                    break;
                    case 8: *r = 14; return (regAux & 0x00000003) + ((regAux>>2) & 0x00000003) + ((regAux>>4) & 0x00000003) 
                                    + ((regAux>>6) & 0x00000003) + ((regAux>>8) & 0x00000003) + ((regAux>>10) & 0x00000003) 
                                    + ((regAux>>12) & 0x00000003) + ((regAux>>14) & 0x00000003) + ((regAux>>16) & 0x00000003) 
                                    + ((regAux>>18) & 0x00000003) + ((regAux>>20) & 0x00000003) + ((regAux>>22) & 0x00000003)
                                    + ((regAux>>24) & 0x00000003) + ((regAux>>26) & 0x00000003) + 13;
                    break;
                    case 9: *r = 28; return (regAux & 0x00000001) + ((regAux>>1) & 0x00000001) + ((regAux>>2) & 0x00000001) 
                                    + ((regAux>>3) & 0x00000001) + ((regAux>>4) & 0x00000001) + ((regAux>>5) & 0x00000001) 
                                    + ((regAux>>6) & 0x00000001) + ((regAux>>7) & 0x00000001) + ((regAux>>8) & 0x00000001) 
                                    + ((regAux>>9) & 0x00000001) + ((regAux>>10) & 0x00000001) + ((regAux>>11) & 0x00000001)
                                    + ((regAux>>12) & 0x00000001) + ((regAux>>13) & 0x00000001) + ((regAux>>14) & 0x00000001) 
                                    + ((regAux>>15) & 0x00000001) + ((regAux>>16) & 0x00000001) + ((regAux>>17) & 0x00000001) 
                                    + ((regAux>>18) & 0x00000001) + ((regAux>>19) & 0x00000001) + ((regAux>>20) & 0x00000001) 
                                    + ((regAux>>21) & 0x00000001) + ((regAux>>22) & 0x00000001) + ((regAux>>23) & 0x00000001)
                                    + ((regAux>>24) & 0x00000001) + ((regAux>>25) & 0x00000001) + ((regAux>>26) & 0x00000001) 
                                    + ((regAux>>27) & 0x00000001) + 27;
                    break;
            }
            return 0;
        }

    public:
        //! Default constructor
        s9_vector() {};

        //! Copy constructor
        s9_vector(const s9_vector& s9){
            copy(s9);
        }

        //! Move constructor
        s9_vector(s9_vector&& s9){
            *this = std::move(s9);
        }

        //! Constructor
        /*!
        *  \param bv  Uncompressed bit_vector.
        */
        s9_vector(const bit_vector& bv){
            m_size = bv.size();
            m = util::cnt_one_bits(bv);
            
            size_type last_block = m % t_bs;
            size_type blocks_size = (m / t_bs) + ((last_block)? 1 : 0);
            m_blocks = t_int_vector(blocks_size, 0);
            m_absolute = t_int_vector(blocks_size, 0);
            
            t_int_vector m_gaps = t_int_vector (m, 0);
            // storing absolute positions of 1s in bv
            for(size_type k = 0, j = 0; k < m_size; k++){
               if (bv[k] == 1) m_gaps[j++] = k;
            }
            // now stores positions in differential form
            // leaves in absolute form the first position in each block
            for(size_type k = m - 1, j = blocks_size - 1; k > 0; k--) {
               if ((k%t_bs)!= 0) {
                  m_gaps[k] -= m_gaps[k-1] + 1; //-+1 reemplaza m_gaps[k]--;
               }
               else {
                  m_absolute[j--] = m_gaps[k];
                  m_gaps[k] = 0;
               }
            }
            m_absolute[0] = m_gaps[0];
            m_gaps[0] = 0;

            m_seq = int_vector<32>(m,0);
            value_type words = 0;
            size_type i = 0;
            for(; i < blocks_size - ((last_block)? 1 : 0); i++){
                m_blocks[i] = words;
                words += pack_in_block(m_gaps, i*t_bs, m_seq, words, t_bs);
            }
            //Storing elements of the last non completed block (temporary)
            if(last_block){
                m_blocks[i] = words;
                words += pack_in_block(m_gaps, i*t_bs, m_seq, words, last_block);
            }

            size_type seq_size = words;
            m_seq.resize(seq_size);

            //filling the rank directory
            size_type block_rank = 64 * t_bs;
            size_type rank_size = m_size/block_rank + ((m_size % block_rank)!=0);
            m_rank = t_int_vector(rank_size);
            m_rank[0] = 0;
            for (size_type i = 1, j = 1; i < rank_size; i++){
                for (;(j < blocks_size - 1) && (m_absolute[j] < (block_rank * i)); j++);
                if (m_absolute[j] == block_rank * i)
                    m_rank[i] = j;
                else
                    m_rank[i] = j-1;
            }
        }

        //! Swap method
        void swap(s9_vector& s9){
            if (this != &s9) {
                std::swap(m, s9.m);
                std::swap(m_size, s9.m_size);
                m_seq.swap(s9.m_seq);
                m_blocks.swap(s9.m_blocks);
                m_absolute.swap(s9.m_absolute);
                m_rank.swap(s9.m_rank);
            }
        }

        //! Accessing the i-th element of the original bit_vector
        /*! \param i An index i with \f$ 0 \leq i < size()  \f$.
           \return The i-th bit of the original bit_vector
        */
        bool operator[](size_type i) {
            assert(i < m_size);
            if(i+1==0) return 0;

            word_type block_rank = 64 * t_bs;
            word_type j = m_rank[i/block_rank]+1;
            word_type buffer[28];

            register word_type regAux;

            word_type w, s, sPrev, r=0, rPrev;

            while(j < m_absolute.size() && i >= m_absolute[j])
                ++j;
            --j;
            
            if (m_absolute[j] > i) return j&&t_bs;
            else if (i==m_absolute[j]) return 1+j*t_bs;
            else {
                word_type Limit;
                if (j==m_absolute.size()-1) Limit = m_seq.size();
                else Limit = m_blocks[j+1];
                w = m_blocks[j];
                regAux = m_seq[w++];
                s = sPrev = m_absolute[j];
                s+= sum_s9_word(regAux, &r);
                for (rPrev = 0; w < Limit && s < i;) {
                    sPrev = s;
                    rPrev = r;
                    regAux = m_seq[w++];
                    switch(regAux >> 28) {
                        case 1: s+= (regAux & 0x0fffffff)+1;
                                ++r;
                                break;
                        case 2: s+= (regAux & 0x00003fff) + ((regAux>>14) & 0x00003fff) + 2;
                                r+= 2;
                                break;
                        case 3: s+= (regAux & 0x000001ff) + ((regAux>>9) & 0x000001ff) + ((regAux>>18) & 0x000001ff) + 3;
                                r+= 3;
                                break;
                        case 4: s+= (regAux & 0x0000007f) + ((regAux>>7) & 0x0000007f) + ((regAux>>14) & 0x0000007f) + ((regAux>>21) & 0x0000007f) + 4;
                                r+= 4;
                                break;
                        case 5: s+= (regAux & 0x0000001f) + ((regAux>>5) & 0x0000001f) + ((regAux>>10) & 0x0000001f) + ((regAux>>15) & 0x0000001f)
                                    + ((regAux>>20) & 0x0000001f) + 5;
                                r+= 5;
                                break;
                        case 6: s+= (regAux & 0x0000000f) + ((regAux>>4) & 0x0000000f) + ((regAux>>8) & 0x0000000f) + ((regAux>>12) & 0x0000000f)
                                    + ((regAux>>16) & 0x0000000f) + ((regAux>>20) & 0x0000000f) + ((regAux>>24) & 0x0000000f) + 7;
                                r+= 7;
                                break;
                        case 7: s+= (regAux & 0x00000007) + ((regAux>>3) & 0x00000007) + ((regAux>>6) & 0x00000007) + ((regAux>>9) & 0x00000007)
                                    + ((regAux>>12) & 0x00000007) + ((regAux>>15) & 0x00000007) + ((regAux>>18) & 0x00000007) + ((regAux>>21) & 0x00000007)
                                    + ((regAux>>24) & 0x00000007) + 9;
                                r+= 9;
                                break;
                        case 8: s+= (regAux & 0x00000003) + ((regAux>>2) & 0x00000003) + ((regAux>>4) & 0x00000003) + ((regAux>>6) & 0x00000003)
                                    + ((regAux>>8) & 0x00000003) + ((regAux>>10) & 0x00000003) + ((regAux>>12) & 0x00000003) + ((regAux>>14) & 0x00000003)
                                    + ((regAux>>16) & 0x00000003) + ((regAux>>18) & 0x00000003) + ((regAux>>20) & 0x00000003) + ((regAux>>22) & 0x00000003)
                                    + ((regAux>>24) & 0x00000003) + ((regAux>>26) & 0x00000003) + 14;
                                r+= 14;
                                break;
                        case 9: s+= (regAux & 0x00000001) + ((regAux>>1) & 0x00000001) + ((regAux>>2) & 0x00000001) + ((regAux>>3) & 0x00000001)
                                    + ((regAux>>4) & 0x00000001) + ((regAux>>5) & 0x00000001) + ((regAux>>6) & 0x00000001) + ((regAux>>7) & 0x00000001)
                                    + ((regAux>>8) & 0x00000001) + ((regAux>>9) & 0x00000001) + ((regAux>>10) & 0x00000001) + ((regAux>>11) & 0x00000001)
                                    + ((regAux>>12) & 0x00000001) + ((regAux>>13) & 0x00000001) + ((regAux>>14) & 0x00000001) + ((regAux>>15) & 0x00000001)
                                    + ((regAux>>16) & 0x00000001) + ((regAux>>17) & 0x00000001) + ((regAux>>18) & 0x00000001) + ((regAux>>19) & 0x00000001)
                                    + ((regAux>>20) & 0x00000001) + ((regAux>>21) & 0x00000001) + ((regAux>>22) & 0x00000001) + ((regAux>>23) & 0x00000001)
                                    + ((regAux>>24) & 0x00000001) + ((regAux>>25) & 0x00000001) + ((regAux>>26) & 0x00000001) + ((regAux>>27) & 0x00000001) + 28;
                                r+= 28;
                                break;
                    }
                }
                if (s > i) {
                    s = sPrev;
                    r = rPrev;
                    buffer[0] = 0;
                    switch (regAux >> 28) {
                        case 1: buffer[0] = regAux & 0x0fffffff;
                                break;
                        case 2: buffer[0] = regAux & 0x00003fff;        buffer[1] = (regAux>>14) & 0x00003fff;
                                break;
                        case 3: buffer[0] = regAux & 0x000001ff;        buffer[1] = (regAux>>9) & 0x000001ff;
                                buffer[2] = (regAux>>18) & 0x000001ff;
                                break;
                        case 4: buffer[0] = regAux & 0x0000007f;        buffer[1] = (regAux>>7) & 0x0000007f;
                                buffer[2] = (regAux>>14) & 0x0000007f;  buffer[3] = (regAux>>21) & 0x0000007f;
                                break;
                        case 5: buffer[0] = regAux & 0x0000001f;        buffer[1] = (regAux>>5) & 0x0000001f;
                                buffer[2] = (regAux>>10) & 0x0000001f;  buffer[3] = (regAux>>15) & 0x0000001f;
                                buffer[4] = (regAux>>20) & 0x0000001f;
                                break;
                        case 6: buffer[0] = regAux & 0x0000000f;        buffer[1] = (regAux>>4) & 0x0000000f;
                                buffer[2] = (regAux>>8) & 0x0000000f;   buffer[3] = (regAux>>12) & 0x0000000f;
                                buffer[4] = (regAux>>16) & 0x0000000f;  buffer[5] = (regAux>>20) & 0x0000000f;
                                buffer[6] = (regAux>>24) & 0x0000000f;
                                break;
                        case 7: buffer[0] = regAux & 0x00000007;        buffer[1] = (regAux>>3) & 0x00000007;
                                buffer[2] = (regAux>>6) & 0x00000007;   buffer[3] = (regAux>>9) & 0x00000007;
                                buffer[4] = (regAux>>12) & 0x00000007;  buffer[5] = (regAux>>15) & 0x00000007;
                                buffer[6] = (regAux>>18) & 0x00000007;  buffer[7] = (regAux>>21) & 0x00000007;
                                buffer[8] = (regAux>>24) & 0x00000007;
                                break;
                        case 8: buffer[0] = regAux & 0x00000003;        buffer[1] = (regAux>>2) & 0x00000003;
                                buffer[2] = (regAux>>4) & 0x00000003;   buffer[3] = (regAux>>6) & 0x00000003;
                                buffer[4] = (regAux>>8) & 0x00000003;   buffer[5] = (regAux>>10) & 0x00000003;
                                buffer[6] = (regAux>>12) & 0x00000003;  buffer[7] = (regAux>>14) & 0x00000003;
                                buffer[8] = (regAux>>16) & 0x00000003;  buffer[9] = (regAux>>18) & 0x00000003;
                                buffer[10] = (regAux>>20) & 0x00000003; buffer[11] = (regAux>>22) & 0x00000003;
                                buffer[12] = (regAux>>24) & 0x00000003; buffer[13] = (regAux>>26) & 0x00000003;
                                break;
                        case 9: buffer[0] = regAux & 0x00000001;        buffer[1] = (regAux>>1) & 0x00000001;
                                buffer[2] = (regAux>>2) & 0x00000001;   buffer[3] = (regAux>>3) & 0x00000001;
                                buffer[4] = (regAux>>4) & 0x00000001;   buffer[5] = (regAux>>5) & 0x00000001;
                                buffer[6] = (regAux>>6) & 0x00000001;   buffer[7] = (regAux>>7) & 0x00000001;
                                buffer[8] = (regAux>>8) & 0x00000001;   buffer[9] = (regAux>>9) & 0x00000001;
                                buffer[10] = (regAux>>10) & 0x00000001; buffer[11] = (regAux>>11) & 0x00000001;
                                buffer[12] = (regAux>>12) & 0x00000001; buffer[13] = (regAux>>13) & 0x00000001;
                                buffer[14] = (regAux>>14) & 0x00000001; buffer[15] = (regAux>>15) & 0x00000001;
                                buffer[16] = (regAux>>16) & 0x00000001; buffer[17] = (regAux>>17) & 0x00000001;
                                buffer[18] = (regAux>>18) & 0x00000001; buffer[19] = (regAux>>19) & 0x00000001;
                                buffer[20] = (regAux>>20) & 0x00000001; buffer[21] = (regAux>>21) & 0x00000001;
                                buffer[22] = (regAux>>22) & 0x00000001; buffer[23] = (regAux>>23) & 0x00000001;
                                buffer[24] = (regAux>>24) & 0x00000001; buffer[25] = (regAux>>25) & 0x00000001;
                                buffer[26] = (regAux>>26) & 0x00000001; buffer[27] = (regAux>>27) & 0x00000001;
                                break;
                    }

                    word_type k;
                    if (r==0) { k = 1; r = 1;}
                    else {k = 0;}

                    while (s < i) {
                        s+= (buffer[k++]+1);
                        r += (s<=i);
                    }
                }
                return (s==i);
            }
        }

        //! Assignment operator
        s9_vector& operator=(const s9_vector& s9){
            if (this != &s9) {
                copy(s9);
            }
            return *this;
        }

        //! Move assignment operator
        s9_vector& operator=(s9_vector&& s9){
            swap(s9);
            return *this;
        }

        //! Returns the size of the original bit vector.
        size_type size()const{
            return m_size;
        }

        //! Serializes the data structure into the given ostream
        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const{
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += write_member(m_size, out, child, "m_size");
            written_bytes += write_member(m, out, child, "m");
            written_bytes += m_seq.serialize(out, child, "m_seq");
            written_bytes += m_blocks.serialize(out, child, "m_blocks");
            written_bytes += m_absolute.serialize(out, child, "m_absolute");
            written_bytes += m_rank.serialize(out, child, "m_rank");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Loads the data structure from the given istream.
        void load(std::istream& in){
            read_member(m_size, in);
            read_member(m, in);
            m_seq.load(in);
            m_blocks.load(in);
            m_absolute.load(in);
            m_rank.load(in);
        }

        #ifdef S9_DEBUG
        //Temp functions for debugging the compressed sequence.
        int_vector<32> get_seq(){
            return m_seq;
        }
        int_vector<64> get_absolute(){
            return m_absolute;
        }
        int_vector<64> get_blocks(){
            return m_blocks;
        }
        int_vector<64> get_rank(){
            return m_rank;
        }
        size_type get_m(){
            return m;
        }
        size_type get_n(){
            return m_size;
        }
        #endif
};
//! rank_support for the s9_vector class
/*!
* \tparam t_b           The bit pattern of size one. (so `0` or `1`)
* \tparam t_bs          The block size of the corresponding s9_vector
* \tparam t_int_vector  s9_vector parameter
*/
template<uint8_t t_b, uint16_t t_bs, typename t_int_vector>
class rank_support_s9{
        static_assert(t_b == 1u or t_b == 0u , "rank_support_s9: bit pattern must be `0` or `1`");
    public:
        typedef s9_vector<t_bs, t_int_vector>              bit_vector_type;
        typedef typename bit_vector_type::size_type     size_type;
        typedef typename bit_vector_type::value_type    value_type;
        typedef typename bit_vector_type::word_type     word_type;
        enum { bit_pat = t_b };
        enum { bit_pat_len = (uint8_t)1 };

    private:
        const bit_vector_type* m_v; //!< Pointer to the rank supported s9_vector


        //rank0 implementation
        size_type rank0(size_type i) const{
            assert(i < m_v->m_size);
            if(i+1==0) return 0;
            return 1+i-rank1(i);
        }

        //implementation for rank1
        size_type rank1(size_type i) const{
            i--;
            //assert(i <= m_v->m_size);
            if(i+1==0) return 0;

            word_type block_rank = 64 * t_bs;
            word_type j = m_v->m_rank[i/block_rank]+1;
            word_type buffer[28];

            register word_type regAux;

            word_type w, s, sPrev, r=0, rPrev;

            while(j < m_v->m_absolute.size() && m_v->m_absolute[j] <= i)
                ++j;
            --j;
            
            if (m_v->m_absolute[j] > i) return j*t_bs;
            else if (i==m_v->m_absolute[j]) return 1+j*t_bs;
            else {
                word_type Limit;
                if (j==m_v->m_absolute.size()-1) Limit = m_v->m_seq.size();
                else Limit = m_v->m_blocks[j+1];
                w = m_v->m_blocks[j];
                regAux = m_v->m_seq[w++];
                s = sPrev = m_v->m_absolute[j];
                s+= m_v->sum_s9_word(regAux, &r);
                for (rPrev = 0; w < Limit && s < i;) {
                    sPrev = s;
                    rPrev = r;
                    regAux = m_v->m_seq[w++];
                    switch(regAux >> 28) {
                        case 1: s+= (regAux & 0x0fffffff)+1;
                                ++r;
                                break;
                        case 2: s+= (regAux & 0x00003fff) + ((regAux>>14) & 0x00003fff) + 2;
                                r+= 2;
                                break;
                        case 3: s+= (regAux & 0x000001ff) + ((regAux>>9) & 0x000001ff) + ((regAux>>18) & 0x000001ff) + 3;
                                r+= 3;
                                break;
                        case 4: s+= (regAux & 0x0000007f) + ((regAux>>7) & 0x0000007f) + ((regAux>>14) & 0x0000007f) + ((regAux>>21) & 0x0000007f) + 4;
                                r+= 4;
                                break;
                        case 5: s+= (regAux & 0x0000001f) + ((regAux>>5) & 0x0000001f) + ((regAux>>10) & 0x0000001f) + ((regAux>>15) & 0x0000001f)
                                    + ((regAux>>20) & 0x0000001f) + 5;
                                r+= 5;
                                break;
                        case 6: s+= (regAux & 0x0000000f) + ((regAux>>4) & 0x0000000f) + ((regAux>>8) & 0x0000000f) + ((regAux>>12) & 0x0000000f)
                                    + ((regAux>>16) & 0x0000000f) + ((regAux>>20) & 0x0000000f) + ((regAux>>24) & 0x0000000f) + 7;
                                r+= 7;
                                break;
                        case 7: s+= (regAux & 0x00000007) + ((regAux>>3) & 0x00000007) + ((regAux>>6) & 0x00000007) + ((regAux>>9) & 0x00000007)
                                    + ((regAux>>12) & 0x00000007) + ((regAux>>15) & 0x00000007) + ((regAux>>18) & 0x00000007) + ((regAux>>21) & 0x00000007)
                                    + ((regAux>>24) & 0x00000007) + 9;
                                r+= 9;
                                break;
                        case 8: s+= (regAux & 0x00000003) + ((regAux>>2) & 0x00000003) + ((regAux>>4) & 0x00000003) + ((regAux>>6) & 0x00000003)
                                    + ((regAux>>8) & 0x00000003) + ((regAux>>10) & 0x00000003) + ((regAux>>12) & 0x00000003) + ((regAux>>14) & 0x00000003)
                                    + ((regAux>>16) & 0x00000003) + ((regAux>>18) & 0x00000003) + ((regAux>>20) & 0x00000003) + ((regAux>>22) & 0x00000003)
                                    + ((regAux>>24) & 0x00000003) + ((regAux>>26) & 0x00000003) + 14;
                                r+= 14;
                                break;
                        case 9: s+= (regAux & 0x00000001) + ((regAux>>1) & 0x00000001) + ((regAux>>2) & 0x00000001) + ((regAux>>3) & 0x00000001)
                                    + ((regAux>>4) & 0x00000001) + ((regAux>>5) & 0x00000001) + ((regAux>>6) & 0x00000001) + ((regAux>>7) & 0x00000001)
                                    + ((regAux>>8) & 0x00000001) + ((regAux>>9) & 0x00000001) + ((regAux>>10) & 0x00000001) + ((regAux>>11) & 0x00000001)
                                    + ((regAux>>12) & 0x00000001) + ((regAux>>13) & 0x00000001) + ((regAux>>14) & 0x00000001) + ((regAux>>15) & 0x00000001)
                                    + ((regAux>>16) & 0x00000001) + ((regAux>>17) & 0x00000001) + ((regAux>>18) & 0x00000001) + ((regAux>>19) & 0x00000001)
                                    + ((regAux>>20) & 0x00000001) + ((regAux>>21) & 0x00000001) + ((regAux>>22) & 0x00000001) + ((regAux>>23) & 0x00000001)
                                    + ((regAux>>24) & 0x00000001) + ((regAux>>25) & 0x00000001) + ((regAux>>26) & 0x00000001) + ((regAux>>27) & 0x00000001) + 28;
                                r+= 28;
                                break;
                    }
                }
                if (s > i) {
                    s = sPrev;
                    r = rPrev;
                    buffer[0] = 0;
                    switch (regAux >> 28) {
                        case 1: buffer[0] = regAux & 0x0fffffff;
                                break;
                        case 2: buffer[0] = regAux & 0x00003fff;        buffer[1] = (regAux>>14) & 0x00003fff;
                                break;
                        case 3: buffer[0] = regAux & 0x000001ff;        buffer[1] = (regAux>>9) & 0x000001ff;
                                buffer[2] = (regAux>>18) & 0x000001ff;
                                break;
                        case 4: buffer[0] = regAux & 0x0000007f;        buffer[1] = (regAux>>7) & 0x0000007f;
                                buffer[2] = (regAux>>14) & 0x0000007f;  buffer[3] = (regAux>>21) & 0x0000007f;
                                break;
                        case 5: buffer[0] = regAux & 0x0000001f;        buffer[1] = (regAux>>5) & 0x0000001f;
                                buffer[2] = (regAux>>10) & 0x0000001f;  buffer[3] = (regAux>>15) & 0x0000001f;
                                buffer[4] = (regAux>>20) & 0x0000001f;
                                break;
                        case 6: buffer[0] = regAux & 0x0000000f;        buffer[1] = (regAux>>4) & 0x0000000f;
                                buffer[2] = (regAux>>8) & 0x0000000f;   buffer[3] = (regAux>>12) & 0x0000000f;
                                buffer[4] = (regAux>>16) & 0x0000000f;  buffer[5] = (regAux>>20) & 0x0000000f;
                                buffer[6] = (regAux>>24) & 0x0000000f;
                                break;
                        case 7: buffer[0] = regAux & 0x00000007;        buffer[1] = (regAux>>3) & 0x00000007;
                                buffer[2] = (regAux>>6) & 0x00000007;   buffer[3] = (regAux>>9) & 0x00000007;
                                buffer[4] = (regAux>>12) & 0x00000007;  buffer[5] = (regAux>>15) & 0x00000007;
                                buffer[6] = (regAux>>18) & 0x00000007;  buffer[7] = (regAux>>21) & 0x00000007;
                                buffer[8] = (regAux>>24) & 0x00000007;
                                break;
                        case 8: buffer[0] = regAux & 0x00000003;        buffer[1] = (regAux>>2) & 0x00000003;
                                buffer[2] = (regAux>>4) & 0x00000003;   buffer[3] = (regAux>>6) & 0x00000003;
                                buffer[4] = (regAux>>8) & 0x00000003;   buffer[5] = (regAux>>10) & 0x00000003;
                                buffer[6] = (regAux>>12) & 0x00000003;  buffer[7] = (regAux>>14) & 0x00000003;
                                buffer[8] = (regAux>>16) & 0x00000003;  buffer[9] = (regAux>>18) & 0x00000003;
                                buffer[10] = (regAux>>20) & 0x00000003; buffer[11] = (regAux>>22) & 0x00000003;
                                buffer[12] = (regAux>>24) & 0x00000003; buffer[13] = (regAux>>26) & 0x00000003;
                                break;
                        case 9: buffer[0] = regAux & 0x00000001;        buffer[1] = (regAux>>1) & 0x00000001;
                                buffer[2] = (regAux>>2) & 0x00000001;   buffer[3] = (regAux>>3) & 0x00000001;
                                buffer[4] = (regAux>>4) & 0x00000001;   buffer[5] = (regAux>>5) & 0x00000001;
                                buffer[6] = (regAux>>6) & 0x00000001;   buffer[7] = (regAux>>7) & 0x00000001;
                                buffer[8] = (regAux>>8) & 0x00000001;   buffer[9] = (regAux>>9) & 0x00000001;
                                buffer[10] = (regAux>>10) & 0x00000001; buffer[11] = (regAux>>11) & 0x00000001;
                                buffer[12] = (regAux>>12) & 0x00000001; buffer[13] = (regAux>>13) & 0x00000001;
                                buffer[14] = (regAux>>14) & 0x00000001; buffer[15] = (regAux>>15) & 0x00000001;
                                buffer[16] = (regAux>>16) & 0x00000001; buffer[17] = (regAux>>17) & 0x00000001;
                                buffer[18] = (regAux>>18) & 0x00000001; buffer[19] = (regAux>>19) & 0x00000001;
                                buffer[20] = (regAux>>20) & 0x00000001; buffer[21] = (regAux>>21) & 0x00000001;
                                buffer[22] = (regAux>>22) & 0x00000001; buffer[23] = (regAux>>23) & 0x00000001;
                                buffer[24] = (regAux>>24) & 0x00000001; buffer[25] = (regAux>>25) & 0x00000001;
                                buffer[26] = (regAux>>26) & 0x00000001; buffer[27] = (regAux>>27) & 0x00000001;
                                break;
                    }

                    word_type k;
                    if (r==0) { k = 1; r = 1;}
                    else {k = 0;}

                    while (s < i) {
                        s+= (buffer[k++]+1);
                        r += (s<=i);
                    }
                }
            return r + j*t_bs;
            }
        }


    public:
        //! Standard constructor
        /*! \param v Pointer to the s9_vector, which should be supported
         */
        explicit rank_support_s9(const bit_vector_type* v=nullptr){
            set_vector(v);
        }
        //! Answers rank queries
        /*! \param i Argument for the length of the prefix v[0..i-1], with \f$0\leq i \leq size()\f$.
           \returns Number of 1-bits in the prefix [0..i-1] of the original bit_vector.
           \par Time complexity
                \f$ \Order{ sample\_rate of the rrr\_vector} \f$
        */
        const size_type rank(size_type i)const{
            return  t_b ? rank1(i) : rank0(i);
        }
        
        //! Short hand for rank(i)
        const size_type operator()(size_type i)const{
            return rank(i);
        }
        //! Returns the size of the original vector
        const size_type size()const{
            return m_v->size();
        }

        //! Set the supported vector.
        void set_vector(const bit_vector_type* v=nullptr){
            m_v = v;
        }

        rank_support_s9& operator=(const rank_support_s9& rs){
            if (this != &rs) {
                set_vector(rs.m_v);
            }
            return *this;
        }

        void swap(rank_support_s9&) { }

        //! Load the data structure from a stream and set the supported vector.
        void load(std::istream&, const bit_vector_type* v=nullptr)
        {
            set_vector(v);
        }

        //! Serializes the data structure into a stream.
        size_type serialize(std::ostream&, structure_tree_node* v=nullptr, std::string name="")const
        {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            structure_tree::add_size(child, 0);
            return 0;
        }

};

//! Select support for the s9_vector class. 
/*
* \tparam t_b           Bit pattern of size one. (so `0` or `1`). Only implemented for 1 at the moment.
* \tparam t_bs          Block size of the corresponding s9_vector
* \tparam t_int_vector  s9_vector parameter
*/
template<uint8_t t_b, uint16_t t_bs, typename t_int_vector>
class select_support_s9{
    static_assert(t_b == 1u or t_b == 0u , "select_support_s9: bit pattern must be `0` or `1`");
    public:
        typedef s9_vector<t_bs, t_int_vector>              bit_vector_type;
        typedef typename bit_vector_type::size_type     size_type;
        typedef typename bit_vector_type::value_type    value_type;
        typedef typename bit_vector_type::word_type     word_type;
        enum { bit_pat = t_b };
        enum { bit_pat_len = (uint8_t)1 };

    private:
        const bit_vector_type* m_v; //!< Pointer to the rank supported s9_vector

        //TODO: select0 implementation
        size_type select0(size_type i) const
        {
            if(i==0) return (size_type)-1;
            if(i>m_v->m_size-m_v->m) return (size_type)-1;
            return i;
        }

        //select1 implementation
        size_type select1(size_type i) const {
            if(i==0) return (size_type)-1;
            if(i>m_v->m) return (size_type)-1;
            
            word_type buffer[t_bs], *buff, *p;
            
            word_type j = ((i-1)%t_bs), s;

            if (j==0) return m_v->m_absolute[(i-1)/t_bs]; //es -1 porq por ejemplo si t_bs=128 e i=128, i pertenece al primer bloque, y no al segundo.
                
            register word_type regAux;
            
            p = buffer;
            
            word_type w = m_v->m_blocks[(i-1)/t_bs]; 
            for (int k = (int)j+1;  k > 0;  ) { 
                regAux = m_v->m_seq[w++]; 
                    
                switch(regAux >> 28) {
                    case 1: *p= regAux & 0x0fffffff; 
                        ++p; --k;
                        break;
                    case 2: *p = regAux & 0x00003fff; 
                            *(p + 1) = (regAux>>14) & 0x00003fff;
                            p += 2; k -= 2;
                            break;
                    case 3: *p = regAux & 0x000001ff;
                            *(p + 1) = (regAux>>9) & 0x000001ff;
                            *(p + 2) = (regAux>>18) & 0x000001ff;
                            p += 3; k -= 3;
                            break;
                    case 4: *p = regAux & 0x0000007f;
                            *(p + 1) = (regAux>>7) & 0x0000007f;
                            *(p + 2) = (regAux>>14) & 0x0000007f; 
                            *(p + 3) = (regAux>>21) & 0x0000007f; 
                            p += 4; k -= 4;
                            break;
                    case 5: *p = regAux & 0x0000001f;
                            *(p + 1) = (regAux>>5) & 0x0000001f;
                            *(p + 2) = (regAux>>10) & 0x0000001f; 
                            *(p + 3) = (regAux>>15) & 0x0000001f; 
                            *(p + 4) = (regAux>>20) & 0x0000001f; 
                            p += 5; k -= 5;
                            break;
                    case 6: *p = regAux & 0x0000000f; 
                            *(p + 1) = (regAux>>4) & 0x0000000f; 
                            *(p + 2) = (regAux>>8) & 0x0000000f; 
                            *(p + 3) = (regAux>>12) & 0x0000000f; 
                            *(p + 4) = (regAux>>16) & 0x0000000f; 
                            *(p + 5) = (regAux>>20) & 0x0000000f; 
                            *(p + 6) = (regAux>>24) & 0x0000000f;
                            p += 7; k -= 7;
                            break;
                    case 7: *p = regAux & 0x00000007;
                            *(p  + 1) = (regAux>>3) & 0x00000007;
                            *(p  + 2) = (regAux>>6) & 0x00000007; 
                            *(p  + 3) = (regAux>>9) & 0x00000007; 
                            *(p  + 4) = (regAux>>12) & 0x00000007;
                            *(p  + 5) = (regAux>>15) & 0x00000007; 
                            *(p  + 6) = (regAux>>18) & 0x00000007;
                            *(p  + 7) = (regAux>>21) & 0x00000007;
                            *(p  + 8) = (regAux>>24) & 0x00000007; 
                            p += 9; k -= 9;
                            break;
                    case 8: *p = regAux & 0x00000003; 
                            *(p  + 1) = (regAux>>2) & 0x00000003;
                            *(p  + 2) = (regAux>>4) & 0x00000003;
                            *(p  + 3) = (regAux>>6) & 0x00000003;
                            *(p  + 4) = (regAux>>8) & 0x00000003;
                            *(p  + 5) = (regAux>>10) & 0x00000003;
                            *(p  + 6) = (regAux>>12) & 0x00000003;
                            *(p  + 7) = (regAux>>14) & 0x00000003;
                            *(p  + 8) = (regAux>>16) & 0x00000003;
                            *(p  + 9) = (regAux>>18) & 0x00000003;
                            *(p  + 10) = (regAux>>20) & 0x00000003;
                            *(p  + 11) = (regAux>>22) & 0x00000003;
                            *(p  + 12) = (regAux>>24) & 0x00000003;
                            *(p  + 13) = (regAux>>26) & 0x00000003; 
                            p += 14; k -= 14;
                            break;
                    case 9: *p = regAux & 0x00000001; *(p + 1) = (regAux>>1) & 0x00000001; 
                            *(p + 2) = (regAux>>2) & 0x00000001; *(p + 3) = (regAux>>3) & 0x00000001; 
                            *(p + 4) = (regAux>>4) & 0x00000001; *(p + 5) = (regAux>>5) & 0x00000001;                   
                            *(p + 6) = (regAux>>6) & 0x00000001; *(p + 7) = (regAux>>7) & 0x00000001; 
                            *(p + 8) = (regAux>>8) & 0x00000001; *(p + 9) = (regAux>>9) & 0x00000001; 
                            *(p + 10) = (regAux>>10) & 0x00000001; *(p + 11) = (regAux>>11) & 0x00000001; 
                            *(p + 12) = (regAux>>12) & 0x00000001; *(p + 13) = (regAux>>13) & 0x00000001; 
                            *(p + 14) = (regAux>>14) & 0x00000001; *(p + 15) = (regAux>>15) & 0x00000001;                   
                            *(p + 16) = (regAux>>16) & 0x00000001; *(p + 17) = (regAux>>17) & 0x00000001; 
                            *(p + 18) = (regAux>>18) & 0x00000001; *(p + 19) = (regAux>>19) & 0x00000001; 
                            *(p + 20) = (regAux>>20) & 0x00000001; *(p + 21) = (regAux>>21) & 0x00000001; 
                            *(p + 22) = (regAux>>22) & 0x00000001; *(p + 23) = (regAux>>23) & 0x00000001; 
                            *(p + 24) = (regAux>>24) & 0x00000001; *(p + 25) = (regAux>>25) & 0x00000001; 
                            *(p + 26) = (regAux>>26) & 0x00000001; *(p + 27) = (regAux>>27) & 0x00000001; 
                            p += 28; k -= 28;
                            break;
                }     
            }  
                
            s = m_v->m_absolute[(i-1)/t_bs];
            buff = buffer;
            

            for(word_type k = 0; k <= j; k+=8, buff+=8) {
            buff[0] += s;
            buff[1]+= buff[0];//+1;
            buff[2]+= buff[1];//+1;
            buff[3]+= buff[2];//+1;
            buff[4]+= buff[3];//+1;
            buff[5]+= buff[4];//+1;
            buff[6]+= buff[5];//+1;
            buff[7]+= buff[6];//+1;
            s = buff[7];//+1;
            }
                
            return buffer[j]+j; //Se suma j para compensar los -1 de la codificacion por gaps.5
        }

    public:
        explicit select_support_s9(const bit_vector_type* v=nullptr){
            set_vector(v);
        }

        //! Answers select queries
        size_type select(size_type i)const {
            return  t_b ? select1(i) : select0(i);
        }

        const size_type operator()(size_type i)const{
            return select(i);
        }

        const size_type size()const{
            return m_v->size();
        }

        void set_vector(const bit_vector_type* v=nullptr){
            m_v = v;
        }

        select_support_s9& operator=(const select_support_s9& rs){
            if (this != &rs) {
                set_vector(rs.m_v);
            }
            return *this;
        }

        void swap(select_support_s9&) { }

        void load(std::istream&, const bit_vector_type* v=nullptr){
            set_vector(v);
        }

        size_type serialize(std::ostream&, structure_tree_node* v=nullptr, std::string name="")const{
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            structure_tree::add_size(child, 0);
            return 0;
        }
};

} // end namespace
#endif
