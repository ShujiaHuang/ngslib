#include <stdexcept>

#include <htslib/hts.h>
#include "ngslib/bam.h"
#include "ngslib/utils.h"

namespace ngslib {

    void Bam::_open(const char *fn, const char *mode) {

        _fp = NULL;
        _itr = NULL;
        _idx = NULL;
        _io_status = 0;  // Everything is OK by default.
        fname = tostring(fn);
        _mode = tostring(mode);

        if ((mode[0] == 'r') && (!is_readable(fn))) {
            throw std::invalid_argument("[bam::Bam:_open] file not found - " + fname);
        }

        _fp = sam_open(fn, mode); // Open a Sam/Bam/Cram file
        if (!_fp) {
            throw std::invalid_argument("[bam::Bam:_open] file open failure.");
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

        if (_idx) hts_idx_destroy(_idx);
        _idx = sam_index_load(_fp, fname.c_str());
        if (!_idx) {
            throw std::invalid_argument("[bam::Bam:index_load] Failed to load index BAM/CRAM "
                                        "file or the index file is not available. Rebuild by "
                                        "samtools index please.");
        }
    }

    // 这个函数是否有线程问题，因为 br 并非局部变量，线程内数据变化，它也会变化？
    int Bam::read(BamRecord &br) {
        if (!_hdr.h()) _hdr = BamHeader(_fp);  // If NULL, inital the BAM header by _fp

        br.init();
        if (_itr) {
            _io_status = sam_itr_next(_fp, _itr, br.b());
        } else {
            _io_status = sam_read1(_fp, _hdr.h(), br.b());
        }

        // set BamRecord to NULL if fail to read data
        if (_io_status < 0) br.set_null();

        return _io_status;
    }

    // set_itr_region 这个函数在使用多线程的时候会不会发生问题？
    // 特别是类成员参数, _fp/_idx 在并行处理时是否存在问题? (htslib/thread_pool.h 参考一下)
    // 最好不要在一份文件中做多线程，而是以文件为单位跑多线程.
    // Create a SAM/BAM/CRAM iterator for one region.
    bool Bam::set_itr_region(const std::string &region) {

        if (!_idx) index_load();  // May not be thread safety?
        if (!_hdr) _hdr = BamHeader(_fp);  // If NULL, initial the BAM header by _fp

        if (_itr) sam_itr_destroy(_itr);

        // An iterator on success; NULL on failure
        _itr = sam_itr_querys(_idx, _hdr.h(), region.c_str());
        if (!_itr) {
            throw std::invalid_argument("[bam::Bam:set_itr_region] Fail "
                                        "to create iterator.");
        }

        return _itr != NULL;
    }

    int Bam::write(const BamRecord &br) {
        _io_status = sam_write1(_fp, _hdr.h(), br.b());
        return _io_status;
    }

    std::ostream &operator<<(std::ostream &os, const Bam &b) {

        if (b) {
            os << b.fname;
        }

        return os;
    }

}  // namespace ngslib

