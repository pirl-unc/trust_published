This docker image sets up the environment to run the version of TRUST that was published:
https://media.nature.com/original/nature-assets/ng/journal/v49/n4/extref/ng.3820-s2.zip

This needs to be run on bam files.  Theyu previously recommended Bowtie2 or MapSlice to run the bams so I would suggest making your bams with those to avoid having to optomize parameters for STAR.
https://bitbucket.org/liulab/trust

Interestingly the program and the readme say all rights reserved while the license is GNU.

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

`srun --pty -c 8 --mem 16g -p docker -w c6145-docker-2-0.local docker run -v /datastore:/datastore:shared -v /home/dbortone:/home/dbortone dockerreg.bioinf.unc.edu:5000/trust_3.0.1:1 bash -c "\
  trust \
    --CoreN=8 \
    --genome=hg38 \
    -H \
    -E \
    -f /home/dbortone/scratch/test_bowtie/test_file_bowtie2.bam \
    -o /home/dbortone/scratch/test_bowtie/"
