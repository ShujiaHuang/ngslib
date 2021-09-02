// C++ codes for dealing BAM data
// Author: Shujia Huang
// Date: 2021-08-22

#ifndef __INCLUDE_NGSLIB_BAM_RECORD_H__
#define __INCLUDE_NGSLIB_BAM_RECORD_H__

#include <iostream>
#include <string>

#include <htslib/sam.h>
#include "ngslib/bam_header.h"


namespace ngslib {

    /*! BASES is defined according to information of bam_get_seq() in sam.h,
     * detail for bam_get_seq() is bellow:
     * @abstract  Get query sequence
     * @param  b  pointer to an alignment
     * @return    pointer to sequence
     *
     * @discussion Each base is encoded in 4 bits: 1 for A, 2 for C, 4 for G,
     * 8 for T and 15 for N. Two bases are packed in one byte with the base
     * at the higher 4 bits having smaller coordinate on the read. It is
     * recommended to use bam_seqi() macro to get the base.
     **/
    static const char _BASES[16] = {' ', 'A', 'C', ' ',
                                    'G', ' ', ' ', ' ',
                                    'T', ' ', ' ', ' ',
                                    ' ', ' ', ' ', 'N'};

    class BamRecord {

    private:
        bam1_t *_b;  // bam record

        /* get the max size of Op in cigar */
        unsigned int _max_cigar_Opsize(const char op) const;

    public:
        BamRecord() : _b(NULL) {}  // initial to be NULL.
        ~BamRecord() { destroy(); }

        BamRecord(const BamRecord &b);

        BamRecord(const bam1_t *b);

        BamRecord &operator=(const BamRecord &b);

        BamRecord &operator=(const bam1_t *b);

        void init();

        void destroy();

        // return the bam record pointer.
        bam1_t *b() const { return _b; }

        // conversion function
        operator bool() const { return bool(_b != NULL); }

        friend std::ostream &operator<<(std::ostream &os, const BamRecord &b);

        /// Inline functions for dealing the FLAG of BAM alignment record
        /* The reads are Pair-end in sequencing, no matter whether it is mapped in a pair */
        bool is_paired() const { return _b && (_b->core.flag & BAM_FPAIRED); }

        /* This is read1 */
        bool is_read1() const { return _b && (_b->core.flag & BAM_FREAD1); }

        /* This is read2 */
        bool is_read2() const { return _b && (_b->core.flag & BAM_FREAD2); }

        /* The read itself is mapped */
        bool is_mapped() const { return _b && ((_b->core.flag & BAM_FUNMAP) == 0); }

        /* The meta read is mapped */
        bool is_mate_mapped() const { return _b && ((_b->core.flag & BAM_FMUNMAP) == 0); }

        /* The read is mapped to the reverse strand (-) */
        bool is_mapped_reverse() const { return _b && (_b->core.flag & BAM_FREVERSE); }

        /* The mate read is mapped to the reverse strand (-) */
        bool is_mate_mapped_reverse() const { return _b && (_b->core.flag & BAM_FMREVERSE); }

        /*  The read is mapped in a proper pair */
        bool is_proper_pair() const { return _b && (_b->core.flag & BAM_FPROPER_PAIR); }

        /* The read is a secondary alignment (not primary) */
        bool is_secondary() const { return _b && (_b->core.flag & BAM_FSECONDARY); }

        /* The mapped read is failed QC */
        bool qc_fail() const { return _b && (_b->core.flag & BAM_FQCFAIL); }

        /* The read is a duplicate */
        bool is_duplicate() const { return _b && (_b->core.flag & BAM_FDUP); }

        /* The read is a supplementary alignment */
        bool is_supplementary() const { return _b && (_b->core.flag & BAM_FSUPPLEMENTARY); }

        // -- End the FLAG inline functions --

        /// functions for CIGAR field

        // -- End the CIGAR functions --

        /// functions need Bam header
        // hts_pos_t is a alisa name of int64_t defined in sam.h.
        hts_pos_t tid_length(const BamHeader &hdr) const {
            return _b ? hdr.seq_length(_b->core.tid) : -1;
        }

        hts_pos_t mate_tid_length(const BamHeader &hdr) const {
            return _b ? hdr.seq_length(_b->core.mtid) : -1;
        }

        /* Get the alignment chromosome of this read */
        std::string tid_name(const BamHeader &hdr) const {
            return _b ? hdr.seq_name(_b->core.tid) : "";
        }

        /* Get the alignment chromosome of mate read */
        std::string mate_tid_name(const BamHeader &hdr) const {
            return is_paired() ? hdr.seq_name(_b->core.mtid) : "";
        }

        /// inline functions for the alignment reference information

        /* Get the full alignment flag of this read */
        uint16_t flag() const { return _b->core.flag; }

        /* Get the id of alignment chromosome, defined by sam_hdr_t */
        // Use `target_name` in sam.h to get the name of chromosome.
        int32_t tid() const { return _b ? _b->core.tid : -1; }

        /* chromosome ID of next read in template, defined by sam_hdr_t */
        int32_t mate_tid() const { return _b ? _b->core.mtid : -2; }

        /* Get the alignment strand of this read, Should be one of '*', '-', '+' */
        char map_strand() const {
            return _b ? (is_mapped_reverse() ? '-' : '+') : '*';
        }

        /* Get the alignment strand of mate read, Should be one of '*', '-', '+' */
        char mate_map_strand() const {
            return _b ? (is_mate_mapped_reverse() ? '-' : '+') : '*';
        }

        /* Get the begin mapped position on the reference genome, 0-based
         * return the hts_pos_t on success, -1 on NULL.
         * */
        hts_pos_t reference_start_pos() const {
            // Return the leftmost position of an alignment on the reference genome, 0-based
            return _b ? _b->core.pos : -1;
        }

        /* Get the end mapped position of the alignment on the reference.
          Calculate the rightmost base position of an alignment on the reference genome.

          @return  The coordinate of the first base after the alignment, 0-based.

          For a mapped read, this is just b->core.pos + bam_cigar2rlen.
          For an unmapped read (either according to its flags or if it has no cigar
          string) or a read whose cigar string consumes no reference bases at all,
          we return b->core.pos + 1 by convention.
         */
        hts_pos_t reference_end_pos() const {
            return _b ? bam_endpos(_b) : -1;
        }

        /* Get the begin mapped position of mate on the reference genome, 0-based
         * return the hts_pos_t on success, -1 on NULL.
         * */
        hts_pos_t mate_reference_start_pos() const {
            // 0-based leftmost coordinate of next read in template
            return _b ? _b->core.mpos : -1;
        }

        /* Get mapping quality */
        int mapq() const { return _b ? _b->core.qual : 0; }

        /* convert CIGAR to a string */
        std::string cigar() const;

        /** Return the number of "aligned bases" exclude the base mark as I, D, N, S, H, and P in CIGAR.
         *
         * BamTools reports AlignedBases, which for example returns the literal strings:
         * 3S5M - CTG
         * 5M - CTAGC
         * 3M1D3M - ATG-TGA
         * 3M1I3M - ATGCTGA
         *
         */
        unsigned int align_length() const;

        /* Get the total number of matched base ('M') in this alignment */
        unsigned int match_size() const;

        /* Get max insertion size of this alignment */
        unsigned int max_insertion_size() const;

        /* Get max deletion size of this alignment */
        unsigned int max_deletion_size() const;

        /* Get insert size */
        hts_pos_t insert_size() const { return is_mapped() ? _b->core.isize : 0; }

        /// Functions for the alignment query information
        /* Get the qname of this read as a string */
        std::string qname() const { return std::string(bam_get_qname(_b)); }

        /* get the length of query sequence */
        int query_length() const { return _b ? _b->core.l_qseq : -1; }

        /* Retrieve the sequencing bases of this read as a string (ACTGN) */
        std::string query_sequence() const;

        /* Retrieve the sequencing qualities of this read as a string
         *
         * @param offset Encoding offset for Phred quality scores. Default 33
         * @return quality scores after converting offset. If first char is empty, returns empty string.
         *
         * */
        std::string query_qual(int offset=33) const;

        /* Calculate the mean sequencing quality of the whole read */
        double mean_qqual() const;

        /* Get the alignment start position on this read, by removing soft-clips. 0-based coordinate
         *
         * @return 0-base position on the read on success, -1 on NULL.
         *
         * */
        int32_t query_start_pos() const;

        /* Get the alignment start position on this read, by removing soft-clips. 0-based coordinate
         * Do it in the reverse orientation.
         * @return 0-base position on the read on success, -1 on NULL.
         * */
        int32_t query_start_pos_reverse() const;

        /* Get the end of the alignment on the read, by removing soft-clips, 1-base
         * @return The last position in 1-base on the read on success, -1 on NULL.
         * */
        int32_t query_end_pos() const;

        /* Get the end of the alignment on the read, by removing soft-clips, 1-base
         * Do it in the reverse orientation.
         * @return The last position in 1-base on the read on success, -1 on NULL.
         * */
        int32_t query_end_pos_reverse() const;

        /// Other useful functions

        /* BamRecord has proper orientation (FR): lower position read is mapped to
         * forward strand(+) the higher one mapped to reverse strand(-).
         * */
        bool is_proper_orientation() const;

        /* has a specific TAG in the alignment or not */
        bool has_tag(const std::string tag) const;

        /** Get a string (Z) tag
         * @param tag Name of the tag. eg "XP"
         * @param s The string to be filled in with the tag information
         * @return the tag value as a string is present, even if empty.
         * */
        std::string get_Z_tag(const std::string tag) const;

        /** Get an int (i) tag
         * @param tag Name of the tag. eg "XP"
         * @param t Value to be filled in with the tag value.
         * @return the tag value as a string is present, even if empty.
         * */
        std::string get_Int_tag(const std::string tag) const;

        /** Get a float (f) tag
         * @param tag Name of the tag. eg "AS"
         * @param t Value to be filled in with the tag value.
         * @return the tag value as a string is present, even if empty.
         * */
        std::string get_Float_tag(const std::string tag) const;

        /** Get a string of either Z, f or i type. Useful if tag type not known
         * at compile time.
         *
         * @param tag Name of the tag. eg "XP"
         * @param s The string to be filled in with the tag information
         * @return the tag value as a string is present, even if empty.
         *
         * */
        std::string get_tag(const std::string tag) const;

        /** Get the read group, first from RG tag , then by qname.
         * @return empty string if no read group found
         */
        std::string read_group() const;

        /// Set data to the alignment to change the status of alignment record
        /* Set QC fail for this alignment read */
        void set_fail() {
            if (is_mapped())
                _b->core.flag |= BAM_FQCFAIL;

            return;
        }
    };  // class BamRecord

}  // namespace ngslib

#endif
