#include <stdexcept>

#include <htslib/hts.h>
#include "ngslib/bam.h"
#include "ngslib/utils.h"

namespace ngslib {

    void Bam::_open(const std::string fn, const std::string mode) {

        _fp = NULL;
        _itr = NULL;
        _idx = NULL;

        _fname = fn;
        _mode = mode;
        _io_status = 0;  // Everything is OK by default.

        if ((mode[0] == 'r') && (!is_readable(fn))) {
            throw std::invalid_argument("[bam.cpp::Bam:_open] file not found - " + _fname);
        }

        _fp = sam_open(fn.c_str(), mode.c_str()); // Open a Sam/Bam/Cram file
        if (!_fp) {
            throw std::invalid_argument("[bam.cpp::Bam:_open] file open failure.");
        }

        return;
    }

    Bam::~Bam() {
        // sam_close function is an alias name of hts_close.
        if (_fp) sam_close(_fp);
        if (_idx) hts_idx_destroy(_idx);
        if (_itr) sam_itr_destroy(_itr);

        _io_status = -1;
    }

    const BamHeader &Bam::header() {
        if (!_hdr) {
            _hdr = BamHeader(_fp);  // call copy constructor of BamHeader.
        }
        return _hdr;
    }

    void Bam::index_load() {

        if (_idx)
            hts_idx_destroy(_idx);

        _idx = sam_index_load(_fp, _fname.c_str());
        if (!_idx) {
            throw std::invalid_argument(
                    "[bam.cpp::Bam:index_load] Failed to load index BAM/CRAM "
                    "file or the index file is not available. Rebuild by "
                    "samtools index please."
            );
        }
    }

    // fetch 这个函数在使用多线程的时候会不会发生问题？尝试多区间处理方式？
    // 特别是类成员参数, _fp/_idx 在并行处理时是否存在问题? (htslib/thread_pool.h 参考一下)
    // 最好不要在一份文件中做多线程，而是以文件为单位跑多线程，从而在根上避免？
    // Create a SAM/BAM/CRAM iterator for one region.
    bool Bam::fetch(const std::string &region) {

        if (!_idx) index_load();  // May not be thread safety?
        if (!_hdr) _hdr = BamHeader(_fp);  // If NULL, set BAM header to _hdr.

        // Reset a iterator, An iterator on success; NULL on failure
        if (_itr) sam_itr_destroy(_itr);
        _itr = sam_itr_querys(_idx, _hdr.h(), region.c_str());

        if (!_itr) {
            throw std::invalid_argument("[bam.cpp::Bam:set_itr_region] Fail to fetch "
                                        "the alignment data in region: " + region);
        }

        return _itr != NULL;
    }

    int Bam::read(BamRecord &br) {

        // If NULL, initial the BAM header by _fp.
        if (!_hdr.h()) _hdr = BamHeader(_fp);

        if (!_itr) {
            _io_status = br.load_read(_fp, _hdr.h());
        } else {
            _io_status = br.next_read(_fp, _itr);
        }

        // Destroy BamRecord and set br to be NULL if fail to read data
        if (_io_status < 0) br.destroy();

        return _io_status;
    }

    std::ostream &operator<<(std::ostream &os, const Bam &b) {

        if (b) {
            os << b._fname;
        }

        return os;
    }

}  // namespace ngslib

