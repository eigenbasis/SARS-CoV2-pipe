import os, subprocess, re, time, sys, argparse
from collections import deque
from pathlib import Path

###This is a wrapper script to submit demultiplexing jobs to hpc
full_path = Path(os.path.dirname(os.path.realpath(__file__))).parents[0] #fastq processing folder
script_dir = os.path.dirname(os.path.realpath(__file__))
parser = argparse.ArgumentParser(description='A script to demultiplex raw bcl files to fastq.gz.')
req_arg_grp = parser.add_argument_group('required arguments') #add optional argument group
req_arg_grp.add_argument('-i', '--input', metavar='\b', help = 'Full path to the run folder containing raw data to demultiplex.', required=True)
req_arg_grp.add_argument('-o', '--output', metavar='\b', help = 'The desired name for the output folder, which will be saved in the "demultiplexed" folder.', required=True)
#req_arg_grp.add_argument('-d', '--dest', metavar='\b', help = 'Full path to the folder where demultiplexed files should be copied (e.g. covid_input or bact_data)', required=False)

if len(sys.argv)==1: #if no command-line arguments provided - display help and stop script excecution
    parser.print_help(sys.stderr)   
    sys.exit(1)
args = parser.parse_args() #args list from command-line input

#If command line options contain forward slash at the end
if args.input[-1]=="/":
    input = args.input[:-1]
if args.output[-1]=="/":
    output = args.output[:-1]
#if args.dest[-1]=="/":
#    dest = args.dest[:-1]

#print(['qsub', '-F', f'{input} {output}', 'run_fastq_screen.sh'])
subprocess.check_call(['qsub', '-F', f'{input} {output} {script_dir}', 'run_bcl2fastq.sh'])
