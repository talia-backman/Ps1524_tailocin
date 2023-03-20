# put all combined*.fna VIBRANT output files into one directory and run the following code:
mash sketch -m 2 -s 10000 -o reference *.fna 
mash info reference.msh

# -m 2 improves results by ignoring single copy k-mers which are more likely to be erroneous
# -s 10000 increases the sketch size from the default of 1,000 to 10,000
# sketch size corresponds to the number of (non-redundant) min-hashes that are kept. Larger sketches will better represent the sequences, but at the cost of larger sketch files and longer comparison times.

mash dist -t reference.msh *.fna > final_mash_results.txt
