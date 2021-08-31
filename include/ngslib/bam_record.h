// C++ codes for dealing BAM data
// Author: Shujia Huang
// Date: 2021-08-22

#ifndef __INCLUDE_NGSLIB_BAM_RECORD_H__
#define __INCLUDE_NGSLIB_BAM_RECORD_H__

#include <iostream>
#include <string>

#include <htslib/sam.h>

namespace ngslib {

    class BamRecord {

    private:
        bam1_t *_b;  // bam record

    public:
        BamRecord(): _b(NULL) {}  // initial to be NULL.
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
    };
}

#endif
