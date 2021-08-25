// The C++ codes for dealing BAM header
// Author: Shujia Huang
// Date: 2021-08-20

#ifndef __INCLUDE_NGSLIB_BAM_HEADER_H__
#define __INCLUDE_NGSLIB_BAM_HEADER_H__

#include <iostream>
#include <string>

#include <htslib/sam.h>


namespace ngslib {

    // Store a header to a BAM file, which also acts as a dictionary of
    // reference sequences with names and lengths.
    class BamHeader {

    private:

        /*! `sam_hdr_t` is defined in htslib/sam.h. The data structure of `sam_hdr_t` is:
         *
         * @abstract Structure for the alignment header.
         *
         * @field n_targets   number of reference sequences
         * @field l_text      length of the plain text in the header (may be zero if
         *                    the header has been edited)
         * @field target_len  lengths of the reference sequences
         * @field target_name names of the reference sequences
         * @field text        plain text (may be NULL if the header has been edited)
         * @field sdict       header dictionary
         * @field hrecs       pointer to the extended header struct (internal use only)
         * @field ref_count   reference count
         * @note The text and l_text fields are included for backwards compatibility.

         These fields may be set to NULL and zero respectively as a side-effect
         of calling some header API functions.  New code that needs to access the
         header text should use the sam_hdr_str() and sam_hdr_length() functions
         instead of these fields.
         */

        sam_hdr_t *_h;   // `bam_hdr_t` is an alisa name of `sam_hdr_t`, so I keep using `sam_hdr_t`.

    public:

        // Initializes a new empty BamHeader with no data.
        BamHeader() : _h(NULL) {}

        /** Read the header from a BAM compressed file.
         *
         * @param fp  File pointer
         * @return    A valid pointer to new header on success, NULL on failure
         *
         * This function works on SAM, BAM and CRAM files.
         *
         */
        BamHeader(const std::string &fn);

        explicit BamHeader(samFile *fp);  // Only explicit conversions allowed.

        // Create a new BamHeader from a raw htslib header, rarely use.
        BamHeader(const sam_hdr_t *hdr) { _h = sam_hdr_dup(hdr); }

        BamHeader(const BamHeader &bh) { _h = sam_hdr_dup(bh._h); }  // copy constructor

        ~BamHeader() { sam_hdr_destroy(_h); }

        BamHeader &operator=(const BamHeader &bh);

        BamHeader &operator=(const sam_hdr_t *hdr);

        BamHeader &operator=(const char *fn);

        BamHeader &operator=(const std::string &fn) { return *this = fn.c_str(); }

        friend std::ostream &operator<<(std::ostream &os, const BamHeader &hd);

        void init() { if (_h) _h = sam_hdr_init(); }

        // Free the memory of set Bam file header pointer to be NULL to save memory.
        void destroy();
        void set_null() { return destroy(); }

        // conversion the Bamheader to be a bool type by determine the _h is NULL or nor.
        operator bool() const { return bool(_h != NULL); }

        // Write BAM header to a BAM file.
        int write(samFile *fp) {
            // samFile is an alias of htsFile which define in sam.h by: `typedef htsFile samFile;`
            return sam_hdr_write(fp, _h);
        }

        // return the `sam_hdr_t` pointer of BAM file header.
        sam_hdr_t *h() const { return _h; }
    };
}

#endif