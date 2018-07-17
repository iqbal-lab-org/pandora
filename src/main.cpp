#include <iostream>
#include <cstring>

int pandora_index(int argc, char *argv[]);
int pandora_walk(int argc, char *argv[]);
int pandora_map(int argc, char *argv[]);
int pandora_compare(int argc, char *argv[]);
int pandora_check_kmergraph(int argc, char *argv[]);
int pandora_get_vcf_ref(int argc, char *argv[]);


static int usage()
{
    std::cerr << "\n"
              << "Program: pandora\n"
	      << "Contact: Rachel Colquhoun <rmnorris@well.ox.ac.uk>\n\n"
	      << "Usage:   pandora <command> [options]\n\n"
              << "Command: index         index PRG sequences from FASTA format\n"
	      << "         walk          outputs a path through the nodes in a GFA corresponding\n"
	      << "                       to input sequence, provided it exists\n"
	      << "         map           identify PRG ordering and sequence from reads for a single sample\n"
	      << "	   compare	 identify and compare the PRG ordering and sequences for a set of samples\n"
	      << "\n"
	      << "Note: To map reads against PRG sequences, you need to first index the\n"
	      << "      PRGs with pandora index\n"
              << std::endl;
	return 1;
}

int main(int argc, char *argv[])
{
	int ret;
	if (argc < 2) return usage();
	if (strcmp(argv[1], "index") == 0) ret = pandora_index(argc-1, argv+1);
	else if (strcmp(argv[1], "walk") == 0) ret = pandora_walk(argc-1, argv+1);
	else if (strcmp(argv[1], "map") == 0) ret = pandora_map(argc-1, argv+1);
	else if (strcmp(argv[1], "compare") == 0) ret = pandora_compare(argc-1, argv+1);
	else if (strcmp(argv[1], "check_kmergraph") == 0) ret = pandora_check_kmergraph(argc-1, argv+1);
	else if (strcmp(argv[1], "get_vcf_ref") == 0) ret = pandora_get_vcf_ref(argc-1, argv+1);
	else {
		fprintf(stderr, "[main] unrecognized command '%s'\n", argv[1]);
		return 1;
	}
	return ret;
}
