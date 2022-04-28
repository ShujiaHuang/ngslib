#include <stdexcept>

#include "ngslib/bam_iterator.h"

namespace ngslib {

    BamIterator::BamIterator(samFile *fp, hts_idx_t *idx,
                             sam_hdr_t *hdr) : _itr(NULL) {
        _fp = fp;
        _idx = idx;
        _hdr = hdr;
    }

    BamIterator::BamIterator(samFile *fp, hts_idx_t *idx, sam_hdr_t *hdr,
                             const std::string &region) : _itr(NULL) {
        _fp = fp;
        _idx = idx;
        _hdr = hdr;
        this->fetch(region);
    }

    BamIterator::BamIterator(const BamIterator &bi) {

        _fp = bi._fp;
        _idx = bi._idx;
        _hdr = bi._hdr;
        _itr = bi._itr;
    }

    BamIterator &BamIterator::operator=(const BamIterator &bi) {

        _fp = bi._fp;
        _idx = bi._idx;
        _hdr = bi._hdr;
        _itr = bi._itr;

        return *this;
    }

    bool BamIterator::fetch(const std::string &region) {

        // Reset a iterator, An iterator on success; NULL on failure
        if (_itr)
            sam_itr_destroy(_itr);

        _itr = sam_itr_querys(_idx, _hdr, region.c_str());
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
            io_status = br.load_read(_fp, _hdr);
        } else {
            io_status = br.next_read(_fp, _itr);
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
        _fp = NULL;
        _idx = NULL;
        _hdr = NULL;
        return;
    }
}