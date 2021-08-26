#include <stdexcept>

#include <htslib/hts.h>
#include "ngslib/bam_record.h"
#include "ngslib/utils.h"


namespace ngslib {

    BamRecord::BamRecord(const BamRecord &b) {
        *this = b;
    }

    BamRecord::BamRecord(const bam1_t *b) {
        *this = b;
    }

    BamRecord &BamRecord::operator=(const BamRecord &b) {

        if (_b) bam_destroy1(_b);
        _b = bam_dup1(b._b);
        return *this;
    }

    BamRecord &BamRecord::operator=(const bam1_t *b) {

        if (_b) bam_destroy1(_b);
        _b = bam_dup1(b);
        return *this;
    }

    void BamRecord::init() {
        if (_b) destroy();
        _b = bam_init1();
    }

    void BamRecord::destroy() {
        bam_destroy1(_b);
        _b = NULL;
    }

    std::ostream &operator<<(std::ostream &os, const BamRecord &b) {
        return os;
    }
}
