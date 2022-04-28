// Author: Shujia Huang
// Date: 2021-09-06
#include <iostream>
#include <string>

#include <ngslib/bam.h>
#include <ngslib/bam_record.h>
#include <ngslib/bam_iterator.h>

using ngslib::Bam;
using ngslib::BamRecord;
using ngslib::BamIterator;

void read_br(BamIterator &b) {

    // reading alignment record
    BamRecord al;
    int read_count = 0;

    // fetch alignments
    while (b.next(al) >= 0) {  // -1 => on end of the file.
        ++read_count;
        std::cout << "* Read count: " << read_count << "\t" << al << "\n";
    }
}

BamIterator ret_bi(Bam &b, std::string region) {

    BamIterator bi;
    bi = BamIterator(b.fp(), b.idx(), b.header().h(), region);
    return bi; // BamIterator(b.fp(), b.idx(), b.header().h(), region);
}

int main() {


    const char *fn1 = "../data/range.bam";
    std::string fn2 = "../data/xx_minimal.sam";

    Bam b(fn1, "r");

    BamIterator bi_0;
    BamIterator bi_1(b.fp(), b.idx(), b.header().h());
    BamIterator bi_2(b.fp(), b.idx(), b.header().h(), "CHROMOSOME_IV:1-50000");
    BamIterator bi_3 = BamIterator(b.fp(), b.idx(), b.header().h());

    bi_0 = ret_bi(b, "CHROMOSOME_IV");

//    read_br(bi_3);
    read_br(bi_2);
//    read_br(bi_1);

    return 0;
}