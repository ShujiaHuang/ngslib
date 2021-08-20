#include <stdexcept>

#include "ngslib/fasta.h"
#include "ngslib/utils.h"


namespace ngslib {

    void Fasta::_load_data(const char *file_name) {
        fai = NULL;

        if (!is_readable(file_name)) {
            throw std::invalid_argument("fasta::Fasta: file not found - " + tostring(file_name));
        }

        fname = tostring(file_name);
        fai = fai_load(file_name);  // load data

        if (!fai) {
            throw std::invalid_argument("fasta::Fasta: index not loaded.");
        }
    }

    Fasta::Fasta(const Fasta &ft) {  // copy constructor
        if (fai) {
            fai_destroy(fai);
        }

        this->_load_data(ft.fname.c_str());   // re-open the FASTA file.
    }

    Fasta & Fasta::operator=(const char *file_name) {
        if (fai) {
            fname.clear();
            fai_destroy(fai);
        }

        this->_load_data(file_name);
        return *this;
    }

    std::string Fasta::fetch(const char *chromosome,
                             const ulong start,   // start: 0-base, end: 0-base.
                             const ulong end) const {
        // check if we have loaded the fasta index
        if (!fai) throw std::invalid_argument("Fasta::fetch index not loaded");
        if (start > end) throw std::invalid_argument("Fasta::fetch the start position must be <= end.");
        if (start < 0) throw std::invalid_argument("Fasta::fetch the start position must be >= 0");

        int length;
        char *f = faidx_fetch_seq(fai, chromosome, start, end, &length);

        if (!f) {
            throw std::invalid_argument("Fasta::fetch - Could not find valid sequence");
        }

        std::string sub_seq(f);
        free(f);

        if (sub_seq.empty()) {
            throw std::invalid_argument("Fasta::fetch - Returning empty query on " + tostring(chromosome) +
                                        ":" + tostring(start) + "-" + tostring(end));
        }
        return sub_seq;
    }

    // Output the filename of FASTA
    std::ostream & operator<<(std::ostream & os, const Fasta & fa) {
        os << fa.fname;
        return os;
    }

}  // namespace ngslib


