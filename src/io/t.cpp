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

    class BamIterator {

    private:

        Bam *_bam;
        hts_itr_t *_itr;     // A SAM/BAM/CRAM iterator for a specify region

        BamIterator(const BamIterator &) = delete;             // reject using copy constructor (C++11 style).
        BamIterator &operator=(const BamIterator &) = delete;  // reject using copy/assignment operator (C++11 style).

    public:

        BamIterator() : _bam(NULL), _itr(NULL) {}

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
        explicit BamIterator(Bam *b, const std::string &region);
        explicit BamIterator(Bam *b);
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

        operator bool() const { return bool(*_bam); }
    };
}

#endif


#include <stdexcept>

#include "ngslib/bam_iterator.h"

namespace ngslib {

    BamIterator::BamIterator(Bam *b) : _itr(NULL) {
        _bam = b;
    }

    BamIterator::BamIterator(Bam *b, const std::string &region) : _itr(NULL) {
        _bam = b;
        this->fetch(region);
    }

    bool BamIterator::fetch(const std::string &region) {

        // Reset a iterator, An iterator on success; NULL on failure
        if (_itr)
            sam_itr_destroy(_itr);

        _itr = sam_itr_querys(_bam->idx(), _bam->header().h(), region.c_str());
        if (!_itr) {
            throw std::invalid_argument("[bam_iterator.cpp::BamIterator:fetch] "
                                        "Fail to fetch the alignment data in "
                                        "region: " + region);
        }

        return _itr != NULL;
    }

    int BamIterator::next(BamRecord &br) {

        int io_status;
        if (!_itr) {
            io_status = br.load_read(_bam->fp(), _bam->header().h());
        } else {
            io_status = br.next_read(_bam->fp(), _itr);
        }

        // Destroy BamRecord and set br to be NULL if fail to read data
        if (io_status < 0)
            br.destroy();

        return io_status;
    }

    void BamIterator::destroy() {

        if (_itr) sam_itr_destroy(_itr);
        _itr = NULL;

        // Set _bam to NULL instead of destructing _bam, which is just a pointer
        // of Bam object.
        _bam = NULL;
        return;
    }
}

