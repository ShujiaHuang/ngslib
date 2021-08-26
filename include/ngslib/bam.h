// The C++ codes for BAM/SAM/CRAM file
// Author: Shujia Huang
// Date: 2021-08-20

#ifndef __INCLUDE_NGSLIB_BAM_H__
#define __INCLUDE_NGSLIB_BAM_H__

#include <iostream>
#include <string>

#include "ngslib/bam_header.h"
#include "ngslib/bam_record.h"

namespace ngslib {

    // A Bam file I/O class
    class Bam {
    private:
        std::string fname;  // input file name
        std::string _mode;  // Mode matching / [rwa][bcefFguxz0-9]* /
        int _io_status;     // I/O status code in read() and write() function

        samFile *_fp;       // sam file pointer, samFile is as the same as htsFile in sam.h
        BamHeader _hdr;     // the sam/bam/cram header
        // BamRecord _brc;     // the sam/bam/cram record (可以考虑不要添加这个类成员)
        hts_itr_t *_itr;    // a SAM/BAM/CRAM iterator for specify region
        hts_idx_t *_idx;  // BAM or CRAM index pointer.

        // call `hts_open` function to open file.
        /*!
          @abstract       Open a sequence data (SAM/BAM/CRAM) or variant data (VCF/BCF)
                          or possibly-compressed textual line-orientated file
          @param fn       The file name or "-" for stdin/stdout. For indexed files
                          with a non-standard naming, the file name can include the
                          name of the index file delimited with HTS_IDX_DELIM
          @param mode     Mode matching / [rwa][bcefFguxz0-9]* /
          @discussion
              With 'r' opens for reading; any further format mode letters are ignored
              as the format is detected by checking the first few bytes or BGZF blocks
              of the file.  With 'w' or 'a' opens for writing or appending, with format
              specifier letters:
                b  binary format (BAM, BCF, etc) rather than text (SAM, VCF, etc)
                c  CRAM format
                f  FASTQ format
                F  FASTA format
                g  gzip compressed
                u  uncompressed
                z  bgzf compressed
                [0-9]  zlib compression level
              and with non-format option letters (for any of 'r'/'w'/'a'):
                e  close the file on exec(2) (opens with O_CLOEXEC, where supported)
                x  create the file exclusively (opens with O_EXCL, where supported)
              Note that there is a distinction between 'u' and '0': the first yields
              plain uncompressed output whereas the latter outputs uncompressed data
              wrapped in the zlib format.
          @example
              [rw]b  .. compressed BCF, BAM, FAI
              [rw]bu .. uncompressed BCF
              [rw]z  .. compressed VCF
              [rw]   .. uncompressed VCF
        */
        void _open(const char *fn, const char *mode);

        Bam(const Bam &b) = delete;  // reject using copy constructor (C++11 style).
        Bam &operator=(const Bam &b) = delete;  // reject using copy/assignment operator (C++11 style).

    public:
        Bam() : _fp(NULL), _itr(NULL), _idx(NULL), _io_status(-1) {}
        ~Bam();

        // @mode matching: [rwa]
        Bam(const char *fn, const char *mode) { _open(fn, mode); }

        Bam(const std::string &fn, const std::string mode) {
            _open(fn.c_str(), mode.c_str());
        }

        // set a new header from outside BamHeader
        // void set_header(const BamHeader &h) { _hdr = h; }

        // return the read-only BAM header
        const BamHeader &header();  // return the Bam file header.

        /// Generate and save an index file
        /** @param fn        Input BAM/etc filename, to which .csi/etc will be added
            @param min_shift Positive to generate CSI, or 0 to generate BAI
            @return  0 if successful, or negative if an error occurred (usually -1; or
                     -2: opening fn failed; -3: format not indexable; -4:
                     failed to create and/or save the index)
        */
        int index_build(int min_shift = 0) {
            return sam_index_build(fname.c_str(), min_shift);
        }

        // load index of BAM or CRAM
        void index_load();

        /// call sam_read1 - Read a record from a file
        /** @param fp   Pointer to the source file
         *  @param h    Pointer to the header previously read (fully or partially)
         *  @param b    Pointer to the record placeholder
         *  @return >= 0 on successfully reading a new record, -1 on end of stream, < -1 on error
         **/
        int read(BamRecord &b);

        /// Create a SAM/BAM/CRAM iterator pointer (hts_itr_t*) for one region.
        /** @param idx     Index
            @param hdr     Header
            @param region  Region specification
            @return An iterator on success; NULL on failure

         Regions are parsed by hts_parse_reg(), and take one of the following forms:

            region          | Outputs
            --------------- | -------------
            REF             | All reads with RNAME REF
            REF:            | All reads with RNAME REF
            REF:START       | Reads with RNAME REF overlapping START to end of REF
            REF:-END        | Reads with RNAME REF overlapping start of REF to END
            REF:START-END   | Reads with RNAME REF overlapping START to END
            .               | All reads from the start of the file
            *               | Unmapped reads at the end of the file (RNAME '*' in SAM)

         The form `REF:` should be used when the reference name itself contains a colon.
         Note that SAM files must be bgzf-compressed for iterators to work.
        **/
        bool set_itr_region(const std::string &region);

        /// call sam_write1 - Write a record to a file
        /** @param fp    Pointer to the destination file
         *  @param h     Pointer to the header structure previously read
         *  @param b     Pointer to the record to be written
         *  @return >= 0 on successfully writing the record, -1 on error
         **/
        int write(const BamRecord &b);

        // For reading: >= 0 on successfully reading a new record,
        //              -1 on end of stream, < -1 on error;
        // For writing: >= 0 on successfully writing the record, -1 on error.
        int io_status() { return _io_status; }

        operator bool() const { return _io_status >= 0; }

        friend std::ostream &operator<<(std::ostream &os, const Bam &b);

    };

}  // namespace ngslib

#endif