cp lamdaChr to /mnt/50tb/publicdata/GRCh38/assembly/fastq

/mnt/50tb/repository/3rdparty/bismark_v0.10.1/bismark_genome_preparation --path_to_bowtie /mnt/50tb/repository/3rdparty/bowtie2/bowtie2-2.1.0/ --bowtie2 --verbose /mnt/50tb/publicdata/GRCh38/assembly/fastq

cd /mnt/50tb/privatedata/non-Adams/Suhaib_Shiels
mkdir data
mv ./sequences/* data/
rm -rf sequences/
cd data
tmux new -s "Suhaib"
bismarkPipeline . /mnt/50tb/publicdata/GRCh38/assembly/fastq auto 60


