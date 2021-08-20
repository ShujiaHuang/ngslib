#include <iostream>
#include <string>

#include <ngslib/fasta.h>

int main() {
    using ngslib::Fasta;

    std::string fn = "../data/tinyfasta.fa";
    std::string c_fn = "../data/tinyfasta2.fa.gz";

    Fasta fa0;                // default constructor function
    Fasta fa(fn);             // constructor function: `Fasta(const string file_name)`
    Fasta fa1(c_fn.c_str());  // constructor function: `Fasta(const char *file_name)`
    Fasta fa2 = fn;           // assignment function
    Fasta fa3 = fa;
    Fasta fa4;
    fa4 = c_fn;
    fa4 = fa;

    std::cout << "The path of FASTA file0: " << fa0 << std::endl;
    std::cout << "The path of FASTA file1: " << fa  << std::endl;
    std::cout << "The path of FASTA file2: " << fa2 << std::endl;
    std::cout << "The path of FASTA file3: " << fa3 << std::endl;
    std::cout << "The path of FASTA file4: " << fa4 << std::endl;
    std::cout << "fetch(\"ref1\", 0, 10) : " << fa.fetch("ref1", 0, 10)   << std::endl;
    std::cout << "fetch(\"ref1\", 10) :    " << fa.fetch("ref1", 10)      << std::endl;
    std::cout << "fetch(\"ref3\") :        " << fa.fetch("ref3")          << std::endl;
    std::cout << "fetch(\"ref1\", 0, 10):  " << fa1.fetch("ref1", 0, 10)  << std::endl;
    std::cout << "fetch(\"ref1\", 0, 1):   " << fa1.fetch("ref1", 0, 1)   << std::endl;
    std::cout << "fetch(\"ref1\", 1, 12):  " << fa2.fetch("ref1", 1, 12)  << std::endl;
    std::cout << "fetch(\"ref1\", 1, 12):  " << fa3.fetch("ref1", 1, 12)  << std::endl;
    std::cout << "fetch(\"ref1\", 1, 100): " << fa4.fetch("ref1", 1, 100) << std::endl;
    std::cout << "sequence_length(\"ref1\"): " << fa.sequence_length("ref1") << std::endl;
//    std::cout << "The sequence: " << fa.fetch("ref1", 12, 10) << std::endl;

    return 0;
}
