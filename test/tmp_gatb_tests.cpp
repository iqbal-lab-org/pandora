#include "gtest/gtest.h"
#include <gatb/gatb_core.hpp>
#include <iostream>


//TEST(GatbTmp, create) {
//    // We get a handle on a fake bank made of 3 sequences.
//    IBank *bank = new BankStrings(
//            "ATCGTACGACGCTAGCTAGCA",
//            "ACTACGTATCGGTATATATTTTCGATCGATCAG",
//            "TGACGGTAGCATCGATCAGGATCGA",
//            NULL
//    );
//    try {
//        // We create the graph with the bank and other options
//        Graph graph = Graph::create(bank, "-kmer-size 5  -abundance-min 1  -out mygraph");
//        // We dump some information about the graph.
//        std::cout << graph.getInfo() << std::endl;
//        // Note: Graph::create will take care about 'bank' object and will delete it if nobody else needs it.
//        // In other words: there is no need here to call 'delete' on 'bank' here.
//    }
//    catch (Exception &e) {
//        std::cerr << "EXCEPTION: " << e.getMessage() << std::endl;
//    }
//}

TEST(GatbFastq, create) {
    std::string fastqPath = "../../test/test_cases/test.fastq";

    // we define a try/catch block in case some method fails
    try {
        // we decalre an input Bank
        IBank *inputBank = Bank::open(fastqPath);

        // create an iterator over this Bank
        Iterator<Sequence> *it = inputBank->iterator();

        // loop over the sequences
        for (it->first(); !it->isDone(); it->next()) {
            // shortcut
            Sequence &seq = it->item();

            // we dump the sequence size and the comment
            std::cout << "[" << seq.getDataSize() << "]" << seq.getComment() << "\n";

            // we dump the sequence
            std::cout << seq.toString() << "\n";

            // we dump the quality score
            std::cout << seq.getQuality() << "\n";
        }
    }
    catch (Exception &e) {
        std::cerr << "EXCEPTION: " << e.getMessage() << "\n";
    }
}