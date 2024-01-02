This assignment is based on a script I made at my UMCG internship.
It compares two (G)VCF files and creates a csv output file showing the crossover between the two.
It was originaly developed to be ran in a Nextflow pipeline and to run for 22 chromosomes.
The program was changed so that it may be ran from the terminal. 

Before running:
Make sure you are SSH connected to a nuc with atleast two cores.

Running:
To run, make sure you are located in the "Assignment6" directory.
Then simply enter the command below.
sbatch assignment6.sh