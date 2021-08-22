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
        h = sam_hdr_read(fp);  // get a BAM header pointer on success, NULL on failure.
        sam_close(fp);
    }

    BamHeader &BamHeader::operator=(const BamHeader &bh) {

        if (h) {
            sam_hdr_destroy(h);
        }

        h = sam_hdr_dup(bh.h);
        return *this;
    }

    BamHeader &BamHeader::operator=(const sam_hdr_t *hdr) {

        if (h) {
            sam_hdr_destroy(h);
        }

        h = sam_hdr_dup(hdr);
        return *this;
    }

    BamHeader &BamHeader::operator=(const char *fn) {

        if (h) {
            sam_hdr_destroy(h);
        }

        samFile *fp = hts_open(fn, "r");
        h = sam_hdr_read(fp);  // get a BAM header pointer on success, NULL on failure.
        sam_close(fp);

        return *this;
    }

    std::ostream &operator<<(std::ostream &os, const BamHeader &hd) {

        if (hd.h)
            os << sam_hdr_str(hd.h);

        return os;
    }

}  // namespace ngslib