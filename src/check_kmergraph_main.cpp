#include <cstring>
#include <cassert>
#include <vector>
#include <iostream>
#include <fstream>
#include "utils.h"
#include "localPRG.h"

using namespace std;

int pandora_check_kmergraph(int argc, char *argv[]) // the "pandora check_kmergraph" comand
{
    // can either provide a prgfile with 1 prg and a sequence file (or top/bottom) and return the path through the prg for each sequence
    // OR can provide a prgfile with multiple sequences in it, and a seq file with the same number, with 1-1 correspondance and check
    // seq n in prg n.
    if (argc < 5) {
        fprintf(stderr, "Usage: pandora check_kmergraph <prg.fa> <seq.fa> <k> <w> <flag>\n");
        return 1;
    }

    // load prg graphs and kmergraphs, for now assume there is only one PRG in this file
    vector<LocalPRG *> prgs;
    read_prg_file(prgs, argv[1]);
    load_PRG_kmergraphs(prgs, stoi(argv[4]), stoi(argv[3]), argv[1]);
    assert(!prgs.empty());
    bool flag = false; // output success/fail rather than the node path

    if (strcmp(argv[2], "--top") == 0) {
        for (uint i = 0; i != prgs.size(); ++i) {
            cout << "Top node path along PRG " << prgs[i]->name << ": ";
            vector<LocalNodePtr> npath = prgs[i]->prg.top_path();
            for (uint j = 0; j != npath.size(); ++j) {
                cout << "->" << npath[j]->id;
            }
            cout << endl;
        }
        /*vector<KmerNodePtr> kpath = prgs[0]->find_kmernodes_on_localnode_path(npath);
        cout << kpath[0]->id;
        for (uint j=1; j != kpath.size(); ++j)
            {
            if (find(kpath[j]->inNodes.begin(), kpath[j]->inNodes.end(), kpath[j-1])!=kpath[j]->inNodes.end())
            {
                    cout << "->" << kpath[j]->id;
            } else {
            cout << "  " << kpath[j]->id;
            }
            }
        cout << endl;*/
        return 0;
    } else if (strcmp(argv[2], "--bottom") == 0) {
        for (uint i = 0; i != prgs.size(); ++i) {
            cout << "Bottom node path along PRG " << prgs[i]->name << ": ";
            vector<LocalNodePtr> npath = prgs[i]->prg.bottom_path();
            for (uint j = 0; j != npath.size(); ++j) {
                cout << "->" << npath[j]->id;
            }
            cout << endl;
        }
        /*vector<KmerNodePtr> kpath = prgs[0]->find_kmernodes_on_localnode_path(npath);
        cout << kpath[0]->id;
        for (uint j=1; j != kpath.size(); ++j)
        {
	    if (find(kpath[j]->inNodes.begin(), kpath[j]->inNodes.end(), kpath[j-1])!=kpath[j]->inNodes.end())
            {
                cout << "->" << kpath[j]->id;
            } else {
                cout << "  " << kpath[j]->id;
            }
        }
        cout << endl;*/
        return 0;
    }
    if (argc > 5 and strcmp(argv[5], "--flag") == 0) {
        flag = true;
    }

    // for each read in readfile,  infer node path along sequence
    vector<LocalNodePtr> npath;
    string name, read, line;
    uint read_num = 0;

    ifstream myfile(argv[2]);
    if (myfile.is_open()) {
        while (getline(myfile, line).good()) {
            if (line.empty() || line[0] == '>') {
                if (!name.empty() && !read.empty()) {
                    if (prgs.size() == 1) {
                        cout << "Node path for read " << read_num << " " << name << " along PRG " << prgs[0]->name
                             << ": ";
                        npath = prgs[0]->prg.nodes_along_string(read);
                        if (npath.empty()) {
                            npath = prgs[0]->prg.nodes_along_string(rev_complement(read));
                        }
                    } else if (read_num < prgs.size()) {
                        cout << "Node path for read " << read_num << " " << name << " along PRG "
                             << prgs[read_num]->name << ": ";
                        npath = prgs[read_num]->prg.nodes_along_string(read);
                        if (npath.empty()) {
                            npath = prgs[read_num]->prg.nodes_along_string(rev_complement(read));
                        }
                    } else {
                        cout << "Different numbers of PRGs and reads, exiting" << endl;
                        break;
                    }
                    if (flag) {
                        if (npath.empty() and read.size() < 300) {
                            cout << "short fail!" << endl;
                        } else if (npath.empty() and read.size() >= 300) {
                            cout << "long fail!" << endl;
                        } else {
                            cout << "success!" << endl;
                        }
                    } else {
                        for (uint j = 0; j != npath.size(); ++j) {
                            cout << "->" << *npath[j];
                        }
                        cout << endl;
                    }
                    /*vector<KmerNodePtr> kpath = prgs[0]->find_kmernodes_on_localnode_path(npath);

                    cout << "kmers on path: " << endl;
                            for (uint j=0; j != kpath.size(); ++j)
                            {
                                cout << kpath[j]->id << " " << kpath[j]->path << endl;
                            }

                        cout << "kmer path: " << kpath[0]->id;
                        for (uint j=1; j != kpath.size(); ++j)
                        {
                      if (find(kpath[j]->inNodes.begin(), kpath[j]->inNodes.end(), kpath[j-1])!=kpath[j]->inNodes.end())
                                {
                                    cout << "->" << kpath[j]->id;
                                } else {
                                    cout << endl << "no edge from " << kpath[j-1]->path << " to " << kpath[j]->path << endl;
                        cout << "outnodes are: " << endl;
                        for (uint n=0; n!= kpath[j-1]->outNodes.size(); ++n)
                        {
                        cout << kpath[j-1]->outNodes[n]->path << endl;
                        }

                                }
                        }
                        cout << endl;*/
                    /*for (uint j=0; j != kpath.size(); ++j)
                            {
                    cout << kpath[j]->khash << endl;
                    }*/
                    read_num += 1;
                }
                name.clear();
                read.clear();
                if (!line.empty()) {
                    name = line.substr(1);
                }
            } else {
                read += line;
            }
        }
        // and last entry
        if (!name.empty() && !read.empty()) {
            if (prgs.size() == 1) {
                cout << "Node path for read " << read_num << " " << name << " along PRG " << prgs[0]->name << ": ";
                npath = prgs[0]->prg.nodes_along_string(read);
                if (npath.empty()) {
                    npath = prgs[0]->prg.nodes_along_string(rev_complement(read));
                }
            } else if (read_num < prgs.size()) {
                cout << "Node path for read " << read_num << " " << name << " along PRG " << prgs[read_num]->name
                     << ": ";
                npath = prgs[read_num]->prg.nodes_along_string(read);
                if (npath.empty()) {
                    npath = prgs[read_num]->prg.nodes_along_string(rev_complement(read));
                }
            } else {
                cout << "Different numbers of PRGs and reads, exiting" << endl;
                myfile.close();
                return 1;
            }
            if (flag) {
                if (npath.empty() and read.size() < 300) {
                    cout << "short fail!" << endl;
                } else if (npath.empty() and read.size() >= 300) {
                    cout << "long fail!" << endl;
                } else {
                    cout << "success!" << endl;
                }
            } else {
                for (uint j = 0; j != npath.size(); ++j) {
                    cout << "->" << *npath[j];
                }
                cout << endl;
            }
            /*vector<KmerNodePtr> kpath = prgs[0]->find_kmernodes_on_localnode_path(npath);

                cout << "kmers on path: " << endl;
                for (uint j=0; j != kpath.size(); ++j)
                {
                    cout << kpath[j]->id << " " << kpath[j]->path << endl;
                }

                cout << "kmer path: " << kpath[0]->id;
                for (uint j=1; j != kpath.size(); ++j)
                {
            if (find(kpath[j]->inNodes.begin(), kpath[j]->inNodes.end(), kpath[j-1])!=kpath[j]->inNodes.end())
                    {
                        cout << "->" << kpath[j]->id;
                    } else {
                cout << endl << "no edge from " << kpath[j-1]->path << " to " << kpath[j]->path << endl;
                        cout << "outnodes are: " << endl;
                        for (uint n=0; n!= kpath[j-1]->outNodes.size(); ++n)
                        {
                            cout << kpath[j-1]->outNodes[n]->path << endl;
                        }

                    }
                }
                cout << endl;*/
            /*for (uint j=0; j != kpath.size(); ++j)
            {
                cout << kpath[j]->khash << endl;
            }*/
        }
        myfile.close();
    } else {
        cerr << "Unable to open sequence file " << argv[2] << endl;
        return 1;
    }
    return 0;
}

