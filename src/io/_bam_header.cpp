#include <stdexcept>

#include <htslib/hts.h>
#include "ngslib/_bam_header.h"
#include "ngslib/utils.h"

namespace ngslib {

    BamHeader::BamHeader(const std::string &fn) {

        if (!is_readable(fn)) {
            throw std::invalid_argument("_bam_header::BamHeader: " + fn + "not found.");
        }

        samFile *fp = hts_open(fn.c_str(), "r");
        _h = sam_hdr_read(fp);  // get a BAM header pointer on success, NULL on failure.
        sam_close(fp);
    }
    BamHeader::BamHeader(samFile *fp) {
        samFile *f = hts_open(fp->fn, "r");
        _h = sam_hdr_read(f);
        sam_close(f);
    }

    BamHeader &BamHeader::operator=(const BamHeader &bh) {

        if (_h) {
            sam_hdr_destroy(_h);
        }

        _h = sam_hdr_dup(bh._h);
        return *this;
    }

    BamHeader &BamHeader::operator=(samFile *fp) {

        if (_h) {
            sam_hdr_destroy(_h);
        }

        samFile *f = hts_open(fp->fn, "r");
        _h = sam_hdr_read(f);
        sam_close(f);
        return *this;
    }

    BamHeader &BamHeader::operator=(const sam_hdr_t *hdr) {

        if (_h) {
            sam_hdr_destroy(_h);
        }

        _h = sam_hdr_dup(hdr);
        return *this;
    }

    BamHeader &BamHeader::operator=(const char *fn) {

        if (_h) {
            sam_hdr_destroy(_h);
        }

        samFile *fp = hts_open(fn, "r");
        _h = sam_hdr_read(fp);  // get a BAM header pointer on success, NULL on failure.
        sam_close(fp);

        return *this;
    }

    std::ostream &operator<<(std::ostream &os, const BamHeader &hd) {

        if (hd._h)
            os << sam_hdr_str(hd._h);

        return os;
    }

}  // namespace ngslib