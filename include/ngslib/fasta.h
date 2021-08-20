// The C++ code for reading FASTA.
// Author: Shujia Huang
// Date: 2021-08-18
#ifndef __INCLUDE_NGSLIB_FASTA_H__
#define __INCLUDE_NGSLIB_FASTA_H__

#include <iostream>
#include <string>

#include <htslib/faidx.h>


namespace ngslib {

    // Identify the FASTA format
    class Fasta {

    private:
        typedef unsigned long ulong;

        std::string fname;
        faidx_t *fai;

        void _load_data(const char* name);

    public:
        // default constructor
        Fasta() : fai(NULL) {}
        Fasta(const char *file_name) { this->_load_data(file_name); }
        Fasta(const std::string &file_name) { this->_load_data(file_name.c_str()); }
        Fasta(const Fasta &);  // copy constructor

        // Destroy the malloc'ed faidx_t index inside object
        ~Fasta() { if (fai) fai_destroy(fai); }

        Fasta & operator=(const char *s);
        Fasta & operator=(const std::string & s) {
            return *this = s.c_str();    // inline definition
        }
        Fasta & operator=(const Fasta &s) { return *this = s.fname; }
        friend std::ostream & operator<<(std::ostream & os, const Fasta & fa);

        int sequence_length(const char *seq_id) const { return faidx_seq_len(fai, seq_id); }
        int sequence_length(const std::string seq_id) const { return sequence_length(seq_id.c_str()); }

        /** fetch a string from the fasta sequence
         * @param chromosome name of the reference to query
         * @param start position.  1. Zero-based
         * @param end position.    2. Zero-based
         *
         * @exception Throws an invalid_argument if start > end, chromosome not found, or seq not found
         * @note This is currently NOT thread safe
         */
        std::string fetch(const char * chromosome, const ulong start, const ulong end) const;
        std::string fetch(const std::string & chromosome, const ulong start, const ulong end) const {
            return fetch(chromosome.c_str(), start, end);
        }

        std::string fetch(const std::string & chromosome, const ulong start) const {
            return fetch(chromosome, start, sequence_length(chromosome));
        }
        std::string fetch(const std::string & chromosome) const {
            return fetch(chromosome, 0, sequence_length(chromosome));
        }
    };  // class Fasta

}  // namespace ngslib

#endif  // #ifndef __INCLUDE_NGSLIB_FASTA_H__
