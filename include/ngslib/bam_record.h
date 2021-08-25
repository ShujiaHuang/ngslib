// C++ codes for dealing BAM data
// Author: Shujia Huang
// Date: 2021-08-22

#ifndef __INCLUDE_NGSLIB_BAM_RECORD_H__
#define __INCLUDE_NGSLIB_BAM_RECORD_H__

#include <iostream>
#include <htslib/sam.h>

namespace ngslib {

    class BamRecord {

    private:
        bam1_t *_b;  // bam record

    public:
        BamRecord(): _b(NULL) {};
        ~BamRecord() { if(_b) {bam_destroy1(_b);} }

    };
}

#endif
