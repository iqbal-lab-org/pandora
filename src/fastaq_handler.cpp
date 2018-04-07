#include <string>
#include <cstring>
#include <iostream>
//#include <fstream>
#include "fastaq_handler.h"

using namespace std;

FastaqHandler::FastaqHandler(const string& filepath) : num_reads_parsed(0) {
    fastaq_file.open(filepath);
    if (not fastaq_file.is_open()) {
        cerr << "Unable to open fastaq file " << filepath << endl;
        exit(EXIT_FAILURE);
    }
}

void FastaqHandler::get_next(){
    if (!line.empty() and (line[0] == '>' or line[0] == '@')) {
        name = line.substr(1);
        ++num_reads_parsed;
        read.clear();
    }

    while (getline(fastaq_file, line).good()){
        if (line.empty() || line[0] == '>' || line[0] == '@') {
            if (!read.empty()) // ok we'll allow reads with no name, removed
            {
                return;
            }
            name = line.substr(1);
            ++num_reads_parsed;
            read.clear();
        } else if (line[0] == '+') {
            //skip this line and the qual score line
            getline(fastaq_file, line);
        } else {
            read += line;
        }
    }
}

void FastaqHandler::skip_next(){
    if (!line.empty() and (line[0] == '>' or line[0] == '@')) {
        ++num_reads_parsed;
    }

    while (getline(fastaq_file, line).good()){
        if (line[0] == '>' || line[0] == '@') {
	    return;
	}
    }
}

void FastaqHandler::get_id(const uint32_t& id){
    if (id < num_reads_parsed) {
        num_reads_parsed = 0;
        name.clear();
        read.clear();
        line.clear();
        fastaq_file.clear();
        fastaq_file.seekg(0, fastaq_file.beg);
    }

    while (id > 1 and num_reads_parsed < id-2) {
	skip_next();
        int c = fastaq_file.peek();
        if (c == EOF)
            break;
    }

    while (num_reads_parsed <= id) {
        get_next();
        int c = fastaq_file.peek();
        if (c == EOF)
            break;
    }
}

void FastaqHandler::close(){
    fastaq_file.close();
}


