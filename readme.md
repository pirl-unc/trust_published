This docker image sets up the environment to run the version of TRUST that was published:  
https://media.nature.com/original/nature-assets/ng/journal/v49/n4/extref/ng.3820-s2.zip  
  
This needs to be run on bam files.  They previously recommended Bowtie2 or MapSlice to run the bams so I would suggest making your bams with those to avoid having to optomize parameters for STAR.  
https://bitbucket.org/liulab/trust  
  
---- Helpful commands ----  
  
Init:  
`cd /datastore/alldata/shiny-server/rstudio-common/dbortone/docker/trust_published/; ls  
my_image=trust_published:1`  
  
Run interactive session as root:  
`srun --pty -c 1 --mem 1g -p dockerbuild docker run -v /datastore:/datastore:shared -it dockerreg.bioinf.unc.edu:5000/${my_image} bash`  

Run interactive session as yoourself:  
`srun --pty -c 1 --mem 1g -p docker docker run -v /datastore:/datastore:shared -v /home/dbortone/scratch:/home/dbortone/scratch -it dockerreg.bioinf.unc.edu:5000/${my_image} bash`  

Build the image:  
`srun --pty -c 2 --mem 1g -w c6145-docker-2-0.local -p docker docker build -t dockerreg.bioinf.unc.edu:5000/${my_image} .`  

Push the image:  
`srun --pty -c 2 --mem 1g -w c6145-docker-2-0.local -p docker docker push dockerreg.bioinf.unc.edu:5000/${my_image}`  

Pull the image to other nodes:  
`srun --pty -c 1 --mem 1g -w r820-docker-2-0.local -p docker docker pull dockerreg.bioinf.unc.edu:5000/${my_image}  
srun --pty -c 1 --mem 1g -w r820-docker-2-1.local -p docker docker pull dockerreg.bioinf.unc.edu:5000/${my_image}  
srun --pty -c 1 --mem 1g -w fc830-docker-2-0.local -p docker docker pull dockerreg.bioinf.unc.edu:5000/${my_image}  
srun --pty -c 1 --mem 1g -w c6100-docker-2-0.local -p dockerbuild docker pull dockerreg.bioinf.unc.edu:5000/${my_image}`  
  
Check build works:  
`srun --pty -c 1 --mem 1g -p docker -w c6145-docker-2-0.local docker run -v /datastore:/datastore:shared dockerreg.bioinf.unc.edu:5000/trust_3.0.1:1 bash -c "export HOME=/tmp; trust --help"`
  
srun --pty -c 8 --mem 16g -p docker bash  
num_threads=8  
cd /home/dbortone/scratch/test_bowtie  
sam_file=test_file.sam  
trust_bam_file=test_file_bowtie2.bam  
samtools view -@ $num_threads -S -b $sam_file > unsorted_$trust_bam_file  
samtools sort -@ $num_threads -o $trust_bam_file unsorted_$trust_bam_file  
samtools index -@ $num_threads $trust_bam_file  
exit  
  
`srun --pty -c 8 --mem 16g -p docker -w c6145-docker-2-0.local docker run -v /datastore:/datastore:shared -v /home/dbortone:/home/dbortone dockerreg.bioinf.unc.edu:5000/trust_3.0.1:1 bash -c "export HOME=/tmp; trust --help"`
  
  
`srun --pty -c 8 --mem 16g -p docker -w c6145-docker-2-0.local docker run -v /datastore:/datastore:shared -v /home/dbortone:/home/dbortone dockerreg.bioinf.unc.edu:5000/trust_published:1 bash -c "\  
  python /opt/TRUST/TRUST.py \  
    -H \  
    -a \  
    -f /home/dbortone/scratch/test_bowtie/test_file_bowtie2.bam \  
    -o /home/dbortone/scratch/test_bowtie/"`  
  
  
Usage: TRUST.py [options]  
  
Options:  
-h, --help  
. . . . . . . .show this help message and exit  
-d DIRECTORY, --directory=DIRECTORY  
. . . . . . . .Input bam directory  
-f FILE, --file=FILE  Input bam file: if given, overwite -d option  
-F FILES, --fileList=FILES  
. . . . . . . .Alternative input: a file containing the full path to  
. . . . . . . .all the files. If given, overwrite -d and -f option  
-m RUNNINGMODE, --mode=RUNNINGMODE  
. . . . . . . .Running mode. Accept Cov and Full. Cov: only report  
. . . . . . . .coverage information on each gene; Full: run full  
. . . . . . . .analysis, slow. Default: Full  
-e ERROR, --error=ERROR  
. . . . . . . .Maximum number of sequencing error per repeating unit  
. . . . . . . .tolerated. Default: 1  
-l LENGTH, --overlaplength=LENGTH  
. . . . . . . .Minimum length of overlap sequence for calling reads  
. . . . . . . .overlap. Default 10  
-a, --fasta           Whether or not output fasta format, only in Full mode.
. . . . . . . .Default False  
-s, --Single          If set True, TRUST will always run in single end mode  
-H, --HeavyChainOnly  To save time, in single-end mode, TRUST only search  
. . . . . . . .for beta and delta chains unless this flag is set.  
-I INSERTTHRESHOLD, --InsertThreshold=INSERTTHRESHOLD  
. . . . . . . .For PE library, when two mates overlap, TRUST cannot  
. . . . . . . .properly detect CDR3s based on mapped mates. Set  
. . . . . . . .larger value to force TRUST to run in PE mode despite  
. . . . . . . .small library insert size. Default 10.  
-o WD, --OutputDirectory=WD  
. . . . . . . .Directory to store intermediate files and final TRUST  
. . . . . . . .report. User must have write privilege. If omitted,  
. . . . . . . .the current directory will be applied.  
-B, --Bcell  
. . . . . . . .B cell receptor inference is currently under active  
. . . . . . . .development.  
  