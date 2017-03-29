#include <cstring>
#include <cassert>
#include <vector>
#include <iostream>
#include <fstream>
#include "utils.h"
#include "localnode.h"
#include "localPRG.h"

using namespace std;

int pandora_check_kmergraph(int argc, char *argv[]) // the "pandora check_kmergraph" comand
{
    if (argc != 5) {
	fprintf(stderr, "Usage: pandora walk <prg.fa> <seq.fa> <k> <w>\n");
	return 1;
    }

    // load prg graphs and kmergraphs, for now assume there is only one PRG in this file
    vector<LocalPRG*> prgs;
    read_prg_file(prgs, argv[1]);
    load_PRG_kmergraphs(prgs, string(argv[1]) + ".k" + string(argv[3]) + ".w" + string(argv[4]));
    assert(prgs.size() > 0);
	
    if (strcmp(argv[2],"--top") == 0)
    {
	vector<LocalNode*> npath = prgs[0]->prg.top_path();
        for (uint j=0; j != npath.size(); ++j)
        {
            cout << "->" << npath[j]->id;
        }
        cout << endl;
	vector<KmerNode*> kpath = prgs[0]->find_kmernodes_on_localnode_path(npath);
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
	cout << endl;
        return 0;
    } else if (strcmp(argv[2],"--bottom") == 0) {
        vector<LocalNode*> npath = prgs[0]->prg.bottom_path();
        for (uint j=0; j != npath.size(); ++j)
        {
            cout << "->" << npath[j]->id;
        }
        cout << endl;
        vector<KmerNode*> kpath = prgs[0]->find_kmernodes_on_localnode_path(npath);
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
        cout << endl;
        return 0;
    }

    // for each read in readfile,  infer node path along sequence
    vector<LocalNode*> npath;
    string name, read, line;

    ifstream myfile (argv[2]);
    if (myfile.is_open())
    {
        while ( getline (myfile,line).good() )
        {
            if (line.empty() || line[0] == '>' )
            {
                if (!name.empty() && !read.empty())
                {
		    cout << name << endl;
		    cout << "node path: ";
		    npath = prgs[0]->prg.nodes_along_string(read);
		    for (uint j=0; j != npath.size(); ++j)
		    {
			cout << "->" << *npath[j];
                    }
		    cout << endl;
		    vector<KmerNode*> kpath = prgs[0]->find_kmernodes_on_localnode_path(npath);
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
        	    cout << endl;
		    for (uint j=0; j != kpath.size(); ++j)
                    {
			cout << kpath[j]->khash << endl;
		    }
                }
                name.clear();
                read.clear();
                if (!line.empty())
                {
                    name = line.substr(1);
                }
            }
            else
            {
                read += line;
            }
        }
        // and last entry
        if (!name.empty() && !read.empty())
        {
	    cout << name;
            npath = prgs[0]->prg.nodes_along_string(read);
            for (uint j=0; j != npath.size(); ++j)
            {
                cout << "->" << npath[j]->id;
            }
	    cout << endl;
	    vector<KmerNode*> kpath = prgs[0]->find_kmernodes_on_localnode_path(npath);
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
            cout << endl;
            for (uint j=0; j != kpath.size(); ++j)
            {
                cout << kpath[j]->khash << endl;
            }
        }
        myfile.close();
    } else {
        cerr << "Unable to open sequence file " << argv[2] << endl;
        return 1;
    }
    return 0;
}

