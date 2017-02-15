#include <cstring>
#include <vector>
#include <iostream>
#include <fstream>
#include "localgraph.h"
#include "localnode.h"

using namespace std;

int pandora_walk(int argc, char *argv[]) // the "pandora walk" comand
{
    if (argc != 3) {
	fprintf(stderr, "Usage: pandora walk <in.gfa> <seq.fa>\n");
	return 1;
    }

    // create graph from gfa
    LocalGraph lg;
    lg.read_gfa(argv[1]);
	
    if (strcmp(argv[2],"--top") == 0)
    {
	vector<LocalNode*> npath = lg.top_path();
        for (uint j=0; j != npath.size(); ++j)
        {
            cout << "->" << npath[j]->id;
        }
        cout << endl;
        return 0;
    } else if (strcmp(argv[2],"--bottom") == 0) {
        vector<LocalNode*> npath = lg.bottom_path();
        for (uint j=0; j != npath.size(); ++j)
        {
            cout << "->" << npath[j]->id;
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
		    cout << name;
		    npath = lg.nodes_along_string(read);
		    for (uint j=0; j != npath.size(); ++j)
		    {
			cout << "->" << npath[j]->id;
                    }
		    cout << endl;
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
            npath = lg.nodes_along_string(read);
            for (uint j=0; j != npath.size(); ++j)
            {
                cout << "->" << npath[j]->id;
            }
	    cout << endl;
        }
        myfile.close();
    } else {
        cerr << "Unable to open sequence file " << argv[2] << endl;
        return 1;
    }
    return 0;
}

