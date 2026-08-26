# Download instructions and files for proper formatting

FOR WINDOWS, follow each step carefully as presented

1.	Open windows powershell with admin permissions and type:
wsl --install
2.	Restart your computer when prompted (this will install a linux distribution on your computer.
3.	Download and open the ubuntu app (just search and download ubuntu from google) and create a username and password
a.	When you make your password it’s secret so you can’t see it typed on the screen.
4.	Run these commands to install miniconda for linux.
a.	You need to use Ctrl + shift + v to paste and then enter to run.
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh (if you get an error or failed after this command it means it’s a differed version for your computer)
bash miniconda.sh

During this you would need to type yes a bunch of times. After this is done close and re-open ubuntu.
5.	Copy this into ubuntu (if it doesn’t work because terms of service have not been accepted: simply copy the 2 commands it gives you to accept these terms):
conda create -n rna_env -c bioconda -c conda-forge -y \
cutadapt hisat2 samtools subread fastqc multiqc \
python=3.10 pandas matplotlib seaborn scikit-learn
conda activate rna_env

Then copy each of these one at a time as it will need to download it, so may take a bit.

HISAT2 index (GRCm38 / mm10)
mkdir -p ~/mouse_hisat2_index
cd ~/mouse_hisat2_index
wget https://cloud.biohpc.swmed.edu/index.php/s/grcm38_tran/download -O grcm38_tran.tar.gz
tar -xzf grcm38_tran.tar.gz

Gene annotation GTF
mkdir -p ~/mouse_reference
cd ~/mouse_reference
wget ftp://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz
gunzip Mus_musculus.GRCm38.102.gtf.gz

6.	Create the pipeline script: 
nano ~/run_full_gene_pipeline.py

7.	Copy and paste pipeline script in, (right click to paste), then ctrl + X, Y, and enter to save.
8.	Create a folder with subfolders of each sample name and each file. (Each sample should have 2 files and you will need to extract the fastq.gz files from the core).
~/raw_data_all/                     <-- the main folder you will create
├── sample1/                         <-- subfolder for a sample
│   ├── sample1_R1.fastq.gz
│   └── sample1_R2.fastq.gz
├── sample3/
│   ├── sample3_R1.fastq.gz
│   └── sample3_R2.fastq.gz
├── sample5/
│   ├── sample5_R1.fastq.gz
│   └── sample5_R2.fastq.gz

Copy this command to of the file directory in the ubuntu terminal (I have bolded my file directory from c drive): cp -r /mnt/c/Users/….

(Change from Users onwards to your directory)

Then we will need to rename and compress files in that folder with this code (please rename accordingly the first line to the last folder of the directory:
cd ~/…..
for d in */; do
cd "$d"
for f in *_1.fq; do mv "$f" "${f%_1.fq}_R1.fq"; done
for f in *_2.fq; do mv "$f" "${f%_2.fq}_R2.fq"; done
gzip *.fq
cd ..
done
9.	Now run the pipeline (keep in mind this might take a while):
python ~/run_full_gene_pipeline.py
10.	Answer each question accorindly but please rearrange for your username on ubuntu:
/home/Vasili/…(Vasili is my username: name of the folder of the directory mentioned 2 steps before)

/home/Vasili/final_results

(just press enter, no need to type anything, it is set in the code)

(just press enter, no need to type anything, it is set in the code)

/home/Vasili/mouse_hisat2_index/grcm38_tran/genome_tran

/home/Vasili/mouse_reference/Mus_musculus.GRCm38.102.gtf

(Usually 4 or 8, depends on your CPU, if you do not know just pick 4 but it will run slightly slower)

Now use your real folder names for each sample and then after the colon put the group such as: sample1:A sample2:B ...

sample2:B, sample3:A, sample4:A, sample5:B, sample6:B, sample7:A, sample8:A


11.	Now you can copy your full gene counts.csv to a directory on your pc of your liking by running this command (this will go to my documents foldier but change accordingly):

cp ~/final_results/full_gene_counts.csv "/mnt/c/Users….."

if you want specific files see below but I would save the SUBREAD output from the terminal at the end somewhere to keep 

Differential expression analysis in anaconda prompt on Windows
1.	Open anaconda prompt and create a new environment and install dependencies as below:
conda create -n de_analysis -c conda-forge -y python=3.10 pandas numpy
conda activate de_analysis
pip install pydeseq2
2.	In the same directory where your full_gene_counts.csv file is, create another csv file. TO do this open notepad and simply create a file as per below with your sample number and condition to whatever was defined in the previous pipeline. TO save this, go save as, select All files and add a .csv extension at the end of the file name (make sure it is saved as metadata.csv):
sample,condition
sample1,male
sample2,female
...
3.	You can now run the script but would need to change line 18 to be what you have called your conditions in the metadata file for instance mine are male and female-simply just change male and female to your groups:
contrast = ("condition", "male", "female")  

You will also need to change the User settings from line 13-15 where your full gene counts, metadata file (which was previously created with sample/condition), and output
4.	Now run the script in your python environment:
cd  (your file directory with the script)
python de_analysis.py

For analysis run Analysis.py, for matplotlib figures of volcano plots, and heatmaps.
For analysis, you will need your deseq2_results.csv, full_gene_counts.csv, and metadata.csv.
Volcano plot
-Can specific two-tiers of significance and log2foldchange if needed
-Custom colours for each parameter
-Custom colours for specific genes
-Gene name conversion
-Labelling and figure size

Heatmap
-Single heatmap from specific top N of genes
-Multi-block heatmap (you will need to upload an excel/csv with gene and group for column titles and then ENSEMBL names and group name
-Can use metadata for group labelling
-Colouring, gene renaming, custom colours, and figure size
