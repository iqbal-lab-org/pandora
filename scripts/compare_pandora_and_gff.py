import logging
import argparse
from BCBio import GFF
from Bio.Seq import Seq
import glob
import subprocess

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("GFF", action="store", type=str,
                        help='GFF file for reference genome')
    parser.add_argument("READ_FQ", action="store", type=str,
                        help='fastq file of nanopore reads')
    parser.add_argument("DIR", action="store", type=str,
                        help='Directory containing pandora output, no trailing /')
    parser.add_argument("-v", "--verbosity", dest='verbosity', action='store_true', help='If flagged, puts logger in DEBUG mode')
    args = parser.parse_args()

    if args.verbosity == True:
        logging.basicConfig(filename='%s/alignment_results.log' %args.DIR,level=logging.DEBUG, format='%(asctime)s %(message)s', datefmt='%d/%m/%Y %I:%M:%S')
    else:
        logging.basicConfig(filename='%s/alignment_results.log' %args.DIR,level=logging.INFO, format='%(asctime)s %(message)s', datefmt='%d/%m/%Y %I:%M:%S')

    output_file = open('%s/alignment_results.txt' %args.DIR, 'w')
    #output_racon_file = open('%s/alignment_racon_results.txt' %args.DIR, 'w')

    logging.info("Loading reference GFF %s", args.GFF)
    in_handle = open(args.GFF)
    limit_info = dict(gff_type = ["CDS"])
    for rec in GFF.parse(in_handle, limit_info=limit_info):
        for feat in rec.features:
            if 'gene' in feat.qualifiers.keys():

                # get RefSeq id
		ref_id = False
                for item in feat.qualifiers['inference']:
                    if "RefSeq" in item:
                        ref_id = item.split(":")[2]
		if not ref_id:
		    continue

		# search in pandora dir for corresponding output
		pandora_files = glob.glob("%s/*%s*_kmlp.fasta" %(args.DIR, ref_id))
		if len(pandora_files) == 1:
		    logging.debug("Found output file %s for gene with RefSeq id %s", pandora_files[0], ref_id)
		elif len(pandora_files) > 1:
		    logging.debug("Found multiple output files for gene with RefSeq id %s", ref_id)
		    continue
		else:
		    continue

		
		# get ref gene sequence
                seq1 = rec.seq[feat.location.nofuzzy_start:feat.location.nofuzzy_end]
                if feat.strand == -1:
                    seq1 = seq1.reverse_complement()

		# compare pandora sequence to ref
		with open(pandora_files[0], 'r') as f:
		    content = f.read().splitlines()
		    output_file.write(content[0]+"\n")
		    logging.debug("echo -e '%s\n%s' | ~/apps/cortex/scripts/analyse_variants/needleman_wunsch/needleman_wunsch --stdin --zam", seq1, content[1])
                    retvalue1 = subprocess.check_output("echo -e '%s\n%s' | ~/apps/cortex/scripts/analyse_variants/needleman_wunsch/needleman_wunsch --stdin --zam" %(seq1,content[1]), shell=True)
		    sd1 = int(retvalue1.split('\n')[3].split()[0])#+int(retvalue1.split('\n')[3].split()[1])
                    seq = Seq(content[1])
		    content[1] = seq.reverse_complement()
		    retvalue2 = subprocess.check_output("echo -e '%s\n%s' | ~/apps/cortex/scripts/analyse_variants/needleman_wunsch/needleman_wunsch --stdin --zam" %(seq1,content[1]), shell=True)
                    sd2 = int(retvalue2.split('\n')[3].split()[0])#+int(retvalue2.split('\n')[3].split()[1])
		    if sd1 <= sd2:
		        output_file.write(retvalue1)
		    else:
			output_file.write(retvalue2)

		# polish with racon
		#retvalue = subprocess.check_output("/data2/apps/minimap/minimap %s %s > %s.paf" %(args.READ_FQ, pandora_files[0], pandora_files[0]), shell=True)
		#if (int(retvalue) != 0):
		#    continue
		#rfile = '.'.join(pandora_files[0].split('.')[:-1]) + "racon.fasta"
		#retvalue = subprocess.check_output("/apps/well/racon/20170606/bin/racon %s %s.paf %s > %s" %(args.READ_FQ, pandora_files[0], pandora_files[0], rfile), shell=True)
		#if (int(retvalue) != 0):
                #    continue

		# create the results file for this too
		#with open(rfile, 'r') as f:
		#    content = f.read().splitlines()
		#    output_racon_file.write(content[0]+"\n")
                #    logging.debug("echo -e '%s\n%s' | ~/apps/cortex/scripts/analyse_variants/needleman_wunsch/needleman_wunsch --stdin --zam", seq1, content[1])
                #    retvalue1 = subprocess.check_output("echo -e '%s\n%s' | ~/apps/cortex/scripts/analyse_variants/needleman_wunsch/needleman_wunsch --stdin --zam" %(seq1,content[1]), shell=True)
                #    sd1 = int(retvalue1.split('\n')[3].split()[0])#+int(retvalue1.split('\n')[3].split()[1])
                #    seq = Seq(content[1])
                #    content[1] = seq.reverse_complement()
                #    retvalue2 = subprocess.check_output("echo -e '%s\n%s' | ~/apps/cortex/scripts/analyse_variants/needleman_wunsch/needleman_wunsch --stdin --zam" %(seq1,content[1]), shell=True)
                #    sd2 = int(retvalue2.split('\n')[3].split()[0])#+int(retvalue2.split('\n')[3].split()[1])
                #    if sd1 <= sd2:
                #        output_racon_file.write(retvalue1)
                #    else:
                #        output_racon_file.write(retvalue2)
		    
			
	

    in_handle.close()
    output_file.close()
    #output_racon_file.close()

if __name__ == "__main__":
    main()

