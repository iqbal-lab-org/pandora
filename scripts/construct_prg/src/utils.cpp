#include <cstring>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cassert>
#include "msarecord.h"

void read_msa_fasta_file(vector<MSARecord*>& msa, const string& filepath)
{
    time_t now;
    uint32_t id = 0;
    string name, read, line, dt, sdt;
    MSARecord *s;

    ifstream myfile (filepath);
    if (myfile.is_open())
    {
        //cout << "Opened msa file: " << filepath << endl;
        uint i = 0;
        while ( getline (myfile,line).good() )
        {
            //cout << "reading line " << i << endl;
            if (line.empty() || line[0] == '>' )
            {
                //cout << "line empty or starts with >" << endl;
                if (!name.empty() && !read.empty())
                {
                    now = time(0);
                    dt = ctime(&now);
                    sdt = dt.substr(0,dt.length()-1);
                    cout << sdt << " Found record " << name << endl;
                    s = new MSARecord(id, name, read);
                    
                    if (s!=nullptr)
                    {   
                        msa.push_back(s);
                        id++;
                    }
                }
                name.clear();
                read.clear();
                if (!line.empty())
                {
                    name = line.substr(1);
                    //cout << "new name " << name << endl;
                }
            }
            else
            {
                read += line;
                //cout << "read starts " << read.substr(0,3) << endl; 
            }
            i++;
        }
        // and last entry
        if (!name.empty() && !read.empty())
        {
            now = time(0);
            dt = ctime(&now);
            sdt = dt.substr(0,dt.length()-1);
            cout << sdt << " Found record " << name << endl;
            s = new MSARecord(id, name, read);
            if (s!=nullptr)
                {msa.push_back(s);}
        }
        now = time(0);
        dt = ctime(&now);
        sdt = dt.substr(0,dt.length()-1);
        cout << sdt <<  " Number of records read from MSA: " << msa.size() << endl;
        myfile.close();
    } else {
        cerr << "Unable to open MSA file " << filepath << " - is it in fasta-like format?" << endl;
        exit (EXIT_FAILURE);
    }
    return;
}
