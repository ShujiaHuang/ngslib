#include <stdexcept>

#include <htslib/hts.h>
#include "ngslib/bam_header.h"
#include "ngslib/utils.h"

namespace ngslib {

    BamHeader::BamHeader(samFile *fp) {
        _h = sam_hdr_read(fp);
    }

    BamHeader::BamHeader(const std::string &fn) {

        if (!is_readable(fn)) {
            throw std::invalid_argument("_bam_header::BamHeader: " + fn + "not found.");
        }

        samFile *fp = hts_open(fn.c_str(), "r");
        _h = sam_hdr_read(fp);  // get a BAM header pointer on success, NULL on failure.
        sam_close(fp);
    }

    BamHeader &BamHeader::operator=(const BamHeader &bh) {

        // release _h pointer if _h is not NULL
        sam_hdr_destroy(_h);
        _h = sam_hdr_dup(bh._h);
        return *this;
    }

    BamHeader &BamHeader::operator=(const sam_hdr_t *hdr) {

        // release _h pointer if _h is not NULL
        sam_hdr_destroy(_h);
        _h = sam_hdr_dup(hdr);
        return *this;
    }

    std::ostream &operator<<(std::ostream &os, const BamHeader &hd) {

        if (hd._h)
            os << sam_hdr_str(hd._h);

        return os;
    }

    void BamHeader::destroy() {
        sam_hdr_destroy(_h);
        _h = NULL;
    }
}  // namespace ngslib