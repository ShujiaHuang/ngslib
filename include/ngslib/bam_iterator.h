// Code for iterator BAM/CRAM
// Author: Shujia Huang
// Date: 2021-09-04

#ifndef __INCLUDE_NGSLIB_BAM_ITERATOR_H__
#define __INCLUDE_NGSLIB_BAM_ITERATOR_H__

#include <iostream>
#include <string>

#include <htslib/sam.h>
#include "ngslib/bam.h"
#include "ngslib/bam_record.h"

namespace ngslib {

    // This BamIterator is designed to be called only by Bam object.
    class BamIterator {

    private:
        hts_itr_t *_itr;  // A SAM/BAM/CRAM iterator for a specify region

        samFile *_fp;     // A pointer of SAM/BAM/CRAM
        hts_idx_t *_idx;  // Index
        sam_hdr_t *_hdr;  // A pointer to the header of SAM/BAM/CRAM

    public:

        BamIterator() : _fp(NULL), _idx(NULL), _hdr(NULL), _itr(NULL) {}

        /// Create a SAM/BAM/CRAM iterator pointer (hts_itr_t*) for one region.
        /** @param region  Region specification
            @return 1 on success; 0 on failure

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
        explicit BamIterator(samFile *fp, hts_idx_t *idx, sam_hdr_t *hdr);
        explicit BamIterator(samFile *fp, hts_idx_t *idx, sam_hdr_t *hdr,
                             const std::string &region);

        BamIterator(const BamIterator &bi);   // copy constructor
        BamIterator &operator=(const BamIterator &bi);

        ~BamIterator() { destroy(); }

        // set a SAM/BAM/CRAM iterator for to a specific region.
        /*
         * @param   region  Region specification
         *
         * Regions are parsed by hts_parse_reg(), and take one of the following forms:

            region          | Outputs
            --------------- | -------------
            REF             | All reads with RNAME REF
            REF:            | All reads with RNAME REF
            REF:START       | Reads with RNAME REF overlapping START to end of REF
            REF:-END        | Reads with RNAME REF overlapping start of REF to END
            REF:START-END   | Reads with RNAME REF overlapping START to END
            .               | All reads from the start of the file
            *               | Unmapped reads at the end of the file (RNAME '*' in SAM)
         *
         * */
        bool fetch(const std::string &region);

        /// Read a record from a file
        /** @param fp   Pointer to the source file
         *  @param h    Pointer to the header previously read (fully or partially)
         *  @param b    Pointer to the record placeholder
         *  @return >= 0 on successfully reading a new record, -1 on end of stream,
         *          < -1 on error
         **/
        int next(BamRecord &br);

        void destroy();
    };
}

#endif