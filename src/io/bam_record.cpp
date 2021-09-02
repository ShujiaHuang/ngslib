#include <stdexcept>
#include <sstream>

#include <htslib/hts.h>
#include "ngslib/bam_record.h"
#include "ngslib/utils.h"


namespace ngslib {

    BamRecord::BamRecord(const BamRecord &b) {
        _b = bam_dup1(b._b);
    }

    BamRecord::BamRecord(const bam1_t *b) {
        _b = bam_dup1(b);
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

    std::ostream &operator<<(std::ostream &os, const BamRecord &r) {
        if (!r._b) return os;

        os << r.qname() << "\t"
           << r.flag() << "\t"
           << r.tid() + 1 << "\t"
           << r.reference_start_pos() + 1 << "\t"  // mapping position +1 to make 1-base coordinate
           << r.mapq() << "\t"
           << r.cigar() << "\t"
           << r.mate_tid() + 1 << "\t"
           << r.mate_reference_start_pos() + 1 << "\t"  // mapping position +1 to make 1-base coordinate
           << r.insert_size() << "\t"
           << r.query_sequence() << "\t"
           << r.query_qual();

        return os;
    }

    std::string BamRecord::cigar() const {

        if (!_b) return "*";  // empty

        std::stringstream cig;

        uint32_t *c = bam_get_cigar(_b);
        for (size_t i = 0; i < _b->core.n_cigar; ++i) {
            // Fetch cigar type BAM_CIGAR_STR by and cigar len by bam_cigar_oplen()
            cig << bam_cigar_oplen(c[i]) << BAM_CIGAR_STR[c[i] & BAM_CIGAR_MASK];
        }

        return cig.str();
    }

    unsigned int BamRecord::align_length() const {

        if (!_b) return 0;

        int length = 0;
        char cg;
        uint32_t *c = bam_get_cigar(_b);
        for (size_t i = 0; i < _b->core.n_cigar; ++i) {
            cg = bam_cigar_opchr(c[i]);
            if (cg == 'M' ||
                cg == '=' ||
                cg == 'X') {
                length += bam_cigar_oplen(c[i]);
            }
        }

        return length;
    }

    unsigned int BamRecord::match_size() const {

        if (!_b) return 0;

        unsigned int m_size = 0;
        uint32_t *c = bam_get_cigar(_b);
        for (size_t i = 0; i < _b->core.n_cigar; i++) {
            if (bam_cigar_opchr(c[i]) == 'M')
                m_size += bam_cigar_oplen(c[i]);
        }

        return m_size;
    }

    unsigned int BamRecord::_max_cigar_Opsize(const char op) const {

        if (!_b) return 0;

        unsigned int max_size = 0;
        uint32_t *c = bam_get_cigar(_b);
        for (size_t i = 0; i < _b->core.n_cigar; i++) {
            if (bam_cigar_opchr(c[i]) == op)
                max_size = std::max(bam_cigar_oplen(c[i]), max_size);
        }

        return max_size;
    }

    unsigned int BamRecord::max_insertion_size() const {
        return _max_cigar_Opsize('I');
    }

    unsigned int BamRecord::max_deletion_size() const {
        return _max_cigar_Opsize('D');
    }

    std::string BamRecord::query_sequence() const {

        if (!_b) return "";

        uint8_t *p = bam_get_seq(_b);
        std::string seq(_b->core.l_qseq, 'N');  // initial by a batch of 'N'
        for (size_t i = 0; i < _b->core.l_qseq; ++i)
            seq[i] = _BASES[bam_seqi(p, i)];

        return seq;

    }

    std::string BamRecord::query_qual(int offset) const {

        if (!_b) return "";

        uint8_t *p = bam_get_qual(_b);
        if (!p) return "";

        std::string qual(_b->core.l_qseq, ' ');
        for (size_t i = 0; i < _b->core.l_qseq; ++i)
            qual[i] = (char) (p[i] + offset);

        return qual;
    }

    double BamRecord::mean_qqual() const {

        if (!_b || (_b->core.l_qseq <= 0))
            return -1;

        double total_phred_score = 0;
        uint8_t *p = bam_get_qual(_b);
        for (size_t i = 0; i < _b->core.l_qseq; ++i)
            total_phred_score += p[i];

        return total_phred_score / _b->core.l_qseq;
    }

    int32_t BamRecord::query_start_pos() const {

        if (!_b) return -1;

        uint32_t *c = bam_get_cigar(_b);
        int32_t p = 0;
        for (size_t i = 0; i < _b->core.n_cigar; ++i) {
            if (bam_cigar_opchr(c[i]) == 'S') {
                p += bam_cigar_oplen(c[i]);
            } else {
                break;
            }
        }
        return p;
    }

    int32_t BamRecord::query_start_pos_reverse() const {

        if (!_b) return -1;

        uint32_t *c = bam_get_cigar(_b);
        int32_t p = 0;

        // loop from the end
        for (size_t i = _b->core.n_cigar - 1; i >= 0; --i) {
            if (bam_cigar_opchr(c[i]) == 'S') {
                p += bam_cigar_oplen(c[i]);
            } else { // not a clip, stop counting
                break;
            }
        }
        return p;
    }

    int32_t BamRecord::query_end_pos() const {

        if (!_b) return -1;

        uint32_t *c = bam_get_cigar(_b);
        int32_t p = 0;
        for (int32_t i = _b->core.n_cigar - 1; i >= 0; --i) { // loop from the end
            if (bam_cigar_opchr(c[i]) == 'S') {
                p += bam_cigar_oplen(c[i]);
            } else { // not a clip, so stop counting
                break;
            }
        }
        return (_b->core.l_qseq - p);

    }

    int32_t BamRecord::query_end_pos_reverse() const {

        if (!_b) return -1;

        uint32_t *c = bam_get_cigar(_b);
        int32_t p = 0;
        for (size_t i = 0; i < _b->core.n_cigar; ++i) {
            if (bam_cigar_opchr(c[i]) == 'S') {
                p += bam_cigar_oplen(c[i]);
            } else { // not a clip, so stop counting
                break;
            }
        }
        return (_b->core.l_qseq - p);
    }

    bool BamRecord::is_proper_orientation() const {

        // _b is NULL
        if (!_b) return false;

        // Get false if mate read mapped on different chromosome
        if (_b->core.tid != _b->core.mtid) return false;

        // Get true if FR: Consider read1 must map in front of read2
        if (_b->core.pos < _b->core.mpos) {
            // present read is mapped in front of meta => read1
            // Return true if read1 is mapped to the forward strand (+) and the mate (read2)
            // is mapped to the reverse one (-).
            return (_b->core.flag & BAM_FREVERSE) == 0 && (_b->core.flag & BAM_FMREVERSE) != 0 ? true : false;
        } else {
            // present read is mapped behind meta => read2
            // Return false if read2 is mapped to forward strand and the mate (read1)
            // is mapped to the reverse one.
            return (_b->core.flag & BAM_FREVERSE) == 0 && (_b->core.flag & BAM_FMREVERSE) != 0 ? false : true;
        }
    }

    bool BamRecord::has_tag(const std::string tag) const {

        if (!_b) return false;

        uint8_t *p = bam_aux_get(_b, tag.c_str());
        return bool(p);
    }

    std::string BamRecord::get_Z_tag(const std::string tag) const {

        std::string tag_str;
        if (has_tag(tag)) {

            uint8_t *p = bam_aux_get(_b, tag.c_str());
            if (*p == 'Z') {
                char *pp = bam_aux2Z(p);
                if (pp) {
                    tag_str = std::string(pp);
                }
            }
        }

        return tag_str;
    }

    std::string BamRecord::get_Int_tag(const std::string tag) const {

        std::string tag_str;
        if (has_tag(tag)) {

            uint8_t *p = bam_aux_get(_b, tag.c_str());
            if (*p == 'i' || *p == 'I' ||
                *p == 'S' || *p == 's' ||
                *p == 'C' || *p == 'c') {

                tag_str = tostring(bam_aux2i(p));
            }
        }

        return tag_str;
    }

    std::string BamRecord::get_Float_tag(const std::string tag) const {

        std::string tag_str;
        if (has_tag(tag)) {

            uint8_t *p = bam_aux_get(_b, tag.c_str());
            if (*p == 'f' || *p == 'd') {
                tag_str = tostring(bam_aux2f(p));
            }
        }

        return tag_str;
    }

    std::string BamRecord::get_tag(const std::string tag) const {

        std::string tag_str;

        tag_str = get_Z_tag(tag);
        if (tag_str.size()) {
            return tag_str;
        }

        tag_str = get_Int_tag(tag);
        if (tag_str.size()) {
            return tag_str;
        }

        tag_str = get_Float_tag(tag);
        if (tag_str.size()) {
            return tag_str;
        }

        return tag_str;
    }

    std::string BamRecord::read_group() const {

        std::string rg;
        if (has_tag("RG")) {
            // try to get from RG tag first
            rg = get_Z_tag("RG");

        } else {
            // try to get the read group tag from qname.
            std::string qn = qname();
            size_t pos = qn.find(":", 0);

            rg = (pos != std::string::npos) ? qn.substr(0, pos) : "";
        }

        // Get the read group, return empty string if no read group found.
        return rg;
    }

}  // namespace ngslib
