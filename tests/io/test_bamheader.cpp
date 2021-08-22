// Author: Shujia Huang
// Date: 2021-08-22
#include <iostream>
#include <string>

#include <ngslib/_bam_header.h>

int main() {

    using ngslib::BamHeader;
    std::string fn1 = "../data/range.bam";
    std::string fn2 = "../data/range.cram";
    std::string fn3 = "../data/no_hdr_sq_1.expected.sam";

    BamHeader bh0;
    BamHeader bh1(fn1);
    BamHeader bh2(fn2);
    BamHeader bh3(bh2.hdr());

    BamHeader bh4 = bh2;
    bh4 = bh3;
    bh4 = bh2.hdr();

    // Test NULL to NULL
    BamHeader bh5;
    BamHeader bh6 = bh5;
    bh6 = bh5;

    // bh7 and bh8 is equally.
    BamHeader bh7 = fn3;
    BamHeader bh8 = BamHeader(fn3);

    std::cout << ">>> The header of NULL: " << bh0 << "\n";
    std::cout << ">>> 1. The header of file: " + fn1 << "\n" << bh1 << "\n";
    std::cout << ">>> 2. The header of file: " + fn2 << "\n" << bh2 << "\n";
    std::cout << ">>> 3. The header of file: bh3\n"  << bh3 << "\n";
    std::cout << ">>> 4. The header of file: bh4\n"  << bh4 << "\n";
    std::cout << ">>> 5. The header of file: bh5\n"  << bh5 << "\n";
    std::cout << ">>> 6. The header of file: bh6\n"  << bh6 << "\n";
    std::cout << ">>> 7. The header of file: bh7\n"  << bh7 << "\n";
    std::cout << ">>> 8. The header of file: bh8\n"  << bh8 << "\n";

    return 0;
}

