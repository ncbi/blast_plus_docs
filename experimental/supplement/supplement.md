This repository is an extension to the Official NCBI BLAST+ Docker Image [repository.](https://github.com/ncbi/blast_plus_docs)  If you are not familiar with Docker or general cloud computing concepts, please complete the tutorials in the aforementioned repository first.


## Appendix D - Proof-of-concept for running BLAST Singularity Image
*Note this is using a modified BLAST+ Docker image.

This section contains documentation for the NCBI BLAST+ command line applications in a [Singularity](https://singularity.lbl.gov/docs-installation) image. We will demonstrate how to use the image to run BLAST analysis on the Google Cloud Platform using a small basic example and a more advanced production-level example. Some basic knowledge of Unix/Linux commands and BLAST+ is useful in completing this tutorial.  

In this section, We will first use the same small example from Section 1 of the official BLAST+ Docker repository mentioned above.

Input data
* Query – 1 sequence, 44 nucleotides, file size 0.2 KB
* Database – 7 sequences, 922 nucleotides, file size 1.7 KB

### System Configurations
Region us-east4 (Northern Virginia)   
Zone us-east4c   
16vCPU 104 GB Memory n1-highmem-16   
New 200 GB standard persistent disk  
Ubuntu 18.04 LTS  

### GCP VM Cost (June 2019 - provided by GCP when VM is created)
$559.96 monthly estimate  
$0.767 hourly   

### Install Singularity
https://www.sylabs.io/guides/3.0/user-guide/installation.html

```
## Section 1. Install Singularity June 2019
## https://www.sylabs.io/guides/3.0/user-guide/installation.html

### Install dependencies
sudo apt-get update && sudo apt-get install -y \
    build-essential \
    libssl-dev \
    uuid-dev \
    libgpgme11-dev \
    squashfs-tools \
    libseccomp-dev \
    pkg-config

### Install GO
export VERSION=1.11 OS=linux ARCH=amd64 && \
    wget https://dl.google.com/go/go$VERSION.$OS-$ARCH.tar.gz && \
    sudo tar -C /usr/local -xzvf go$VERSION.$OS-$ARCH.tar.gz && \
    rm go$VERSION.$OS-$ARCH.tar.gz

echo 'export GOPATH=${HOME}/go' >> ~/.bashrc && \
    echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> ~/.bashrc && \
    source ~/.bashrc

go get -u github.com/golang/dep/cmd/dep
go get -d github.com/sylabs/singularity
### Go will complain that there are no Go files, but it will still download the Singularity source code to the appropriate directory within the $GOPATH.

export VERSION=v3.0.3 # or another tag or branch if you like && \
    cd $GOPATH/src/github.com/sylabs/singularity && \
    git fetch && \
    git checkout $VERSION # omit this command to install the latest bleeding edge code from master

### By default Singularity will be installed in the /usr/local directory hierarchy. 
./mconfig && \
    make -C ./builddir && \
    sudo make -C ./builddir install

### Optional bash completion
. /usr/local/etc/bash_completion.d/singularity

### Confirm installation
singularity --version
### output singularity version 3.0.3
```
### Run BLAST+ Singularity Image Using a Small Example

```
## blastsing image has edirect commands copied from /root/edirect to /blast/bin	
cd ~
singularity pull docker://stevetsa/blastsing	
singularity run blastsing_latest.sif

### inside container # specify full paths
blastn -version
### Output 
### blastn: 2.9.0+
### Package: blast 2.9.0, build Apr  2 2019 17:23:44
efetch -version
### Output 10.4

## Section 2. Run BLAST with Singularity using a small example

cd ; mkdir blastdb queries fasta results blastdb_custom
efetch -db protein -format fasta -id P01349 > ~/queries/P01349.fsa
efetch -db protein -format fasta -id Q90523,P80049,P83981,P83982,P83983,P83977,P83984,P83985,P27950 > ~/fasta/nurse-shark-proteins.fsa
cd blastdb_custom
makeblastdb -in ~/fasta/nurse-shark-proteins.fsa -dbtype prot -parse_seqids -out nurse-shark-proteins -title "Nurse shark proteins" -taxid 7801 -blastdb_version 5
cd ~
blastp -query queries/P01349.fsa -db blastdb_custom/nurse-shark-proteins

## output on screen

## Exit container
exit

```

### Run BLAST+ Singularity Image at Production Scale
  
One of the promises of cloud computing is scalability. In this section, we will demonstrate how to use the BLAST+ Singularity image at production scale on the Google Cloud Platform. We will perform a BLAST analysis similar to the approach described in this [publication](https://www.ncbi.nlm.nih.gov/pubmed/31040829) to compare de novo aligned contigs from bacterial 16S-23S sequencing against the nucleotide collection (nt) database.

To test scalability, we will use inputs of different sizes to estimate the amount of time to download the nucleotide collection database and run BLAST search using the latest version of the BLAST+ Docker image. Expected results are summarized in the following tables.

Input files: 28 samples (multi-FASTA files) containing de novo aligned contigs from the publication.  
(Instructions to [download]((https://figshare.com/s/729b346eda670e9daba4)) and create the input files are described in the [code block](#commands) below.)    
  
Database: Pre-formatted BLAST nucleotide collection database, version 5 (nt_v5): 68.7217 GB  (Database size could vary because of regular updates)
  
```
## Section 3. RUN BLAST with Singularity at production scale
## Install Singularity if not already done
## This section assumes using recommended hardware requirements below
## 16 CPUs, 104 GB memory and 200 GB persistent hard disk

## Modify the number of CPUs (-num_threads) in Step 3 if another type of VM is used.

## Step 1. Prepare for analysis
## Create directories
cd ; mkdir -p blastdb queries fasta results blastdb_custom

## Import and process input sequences
sudo apt install unzip
wget https://ndownloader.figshare.com/articles/6865397?private_link=729b346eda670e9daba4 -O fa.zip
unzip fa.zip -d fa

### Create three input query files
### All 28 samples
cat fa/*.fa > query.fa

### Sample 1
cat fa/'Sample_1 (paired) trimmed (paired) assembly.fa' > query1.fa

### Sample 1 to Sample 5
cat fa/'Sample_1 (paired) trimmed (paired) assembly.fa' \
    fa/'Sample_2 (paired) trimmed (paired) assembly.fa' \
    fa/'Sample_3 (paired) trimmed (paired) assembly.fa' \
    fa/'Sample_4 (paired) trimmed (paired) assembly.fa' \
    fa/'Sample_5 (paired) trimmed (paired) assembly.fa' > query5.fa
    
### Copy query sequences to $HOME/queries folder
cp query* $HOME/queries/.

### Start singularity container
singularity run blastsing_latest.sif

### Inside container

## Step 2. Display BLAST databases on the GCP
update_blastdb.pl --showall pretty --source gcp

## Download nt_v5 (nucleotide collection version 5) database
## This step takes approximately 7 min.  The following command runs in the background.
cd blastdb
update_blastdb.pl --source gcp nt_v5 &
cd ..

## At this point, confirm query/database have been properly provisioned before proceeding

## Check the size of the directory containing the BLAST database
## nt_v5 should be around 68 GB
du -sk $HOME/blastdb

## Check for queries, there should be three files - query.fa, query1.fa and query5.fa
ls -al $HOME/queries

## From this point forward, it may be easier if you run these steps in a script. 
## A sample script is available in the next code block

## Step 3. Run BLAST
## Run BLAST using query1.fa (Sample 1) 
## This command will take approximately 8 minutes to complete.
## Expected output size: 3.1 GB  
blastn -query $HOME/queries/query1.fa -db $HOME/blastdb/nt_v5 -num_threads 16 -out $HOME/results/blastn.query1.denovo16s.out

## Run BLAST using query5.fa (Samples 1-5) 
## This command will take approximately 27  minutes to complete.
## Expected output size: 10.4 GB  
blastn -query $HOME/queries/query5.fa -db $HOME/blastdb/nt_v5 -num_threads 16 -out $HOME/results/blastn.query5.denovo16s.out

## Run BLAST using query.fa (All 28 samples) 
## This command will take approximately 145 minutes to complete.
## Expected output size: 47.8 GB  
blastn -query $HOME/queries/query.fa -db $HOME/blastdb/nt_v5 -num_threads 16 -out $HOME/results/blastn.query.denovo16s.out

## Stdout and stderr will be in script.out
## BLAST output will be in $HOME/results

## Exit Singularity container
exit
```

Alternatively, you can run commands in the singularity image as if it were running directly on the host system.  Unlike Docker, you do not need to have a running container. For example, you can copy and paste all the commands below into a file called `script.sh` then run the script in the background using `nohup bash script.sh > script.out &`.  

Sample script.sh - 
```
#singularity exec blastsing_latest.sif blastn -version
date
echo "start Analysis 1"
singularity exec blastsing_latest.sif \
    blastn -query $HOME/queries/query1.fa -db $HOME/blastdb/nt_v5 -num_threads 16 \
    -out $HOME/results/blastn.query1.denovo16s.out
date
echo "start Analysis 2"
singularity exec blastsing_latest.sif \
    blastn -query $HOME/queries/query5.fa -db $HOME/blastdb/nt_v5 -num_threads 16 \
    -out $HOME/results/blastn.query5.denovo16s.out
date
echo "start Analysis 3"
singularity exec blastsing_latest.sif \
    blastn -query $HOME/queries/query.fa -db $HOME/blastdb/nt_v5 -num_threads 16 \
    -out $HOME/results/blastn.query.denovo16s.out
date
echo "Done"
```


### Stop the GCP instance
Remember to [stop](https://cloud.google.com/compute/docs/instances/stop-start-instance) or [delete](https://cloud.google.com/compute/docs/instances/stop-start-instance) the VM to prevent incurring additional cost. You can do this at the GCP Console.


## Appendix E - Proof-of-concept for a BLAST Jupyter Notebook
*Note this is using a modified BLAST+ Docker image.
  
The Jupyter Notebook is a great way to combine free text description with code in the same space and the notebook can be easily shared and reproduced.  In this section, it is assumed that you are familiar with the instructions to create a VM [in the main documentation](https://github.com/ncbi/blast_plus_docs/blob/master/README.md#google-cloud-platform-setup).  The following describes one way to run a Jupyter Notebook server on the GCP.

*The following Steps 1 and 2 are covered in more detail in the NIH Data Science Cookbook. (Availability and location of the Cookbook to be determined)

### Step 1. Create new firewall rules
1. On the console dashboard, expand the navigation menu by clicking on the blue menu button symbol on the left-hand side (if you hover over the symbol, it will label itself as “Navigation menu.”
2. Under the "Networking" section, hover over "VPC Network" and click "Firewall rules"
3. On the upper toolbar, click "Create Firewall Rule"
4. Create a firewall rule to allow all egress traffic
    Under "Name", enter "open-all-egress"
    Under "Direction of traffic", select "Egress"
    Under "Action on match", select "Allow"
    Under "Target tags, enter "open-all-egress"
    Under "Destination IP ranges", enter "0.0.0.0/0"
    Under "Protocols and ports", select "Allow all"
    Leave all other customization options as the default
    Click "Create"
5. Create a firewall rule to allow all ingress traffic
    Under "Name", enter "open-all-ingress"
    Under "Direction of traffic", select "Ingress"
    Under "Action on match", select "Allow"
    Under "Target tags, enter "open-all-ingress"
    Under "Source IP ranges", enter "0.0.0.0/0"
    Under "Protocols and ports", select "Allow all"
    Leave all other customization options as the default
    Click "Create"

### Step 2. Create and set up a VM with the ingress/egress settings
1. On the console dashboard, expand the navigation menu by clicking on the blue menu button symbol on the left-hand side (if you hover over the symbol, it will label itself as “navigation menu.”
2. Hover over “Compute Engine” and click “VM instances”.
3. Create a VM with the following settings
       * Region: us-east4 (Northern Virginia)
       * Machine Type: 16 vCPU, 104 GB memory, n1-highmem-16
       * Boot Disk: Click "Change" and select Ubuntu 18.04 LTS, change the "Boot disk size" to 200 GB Standard persistent disk
4. Under Firewall, check "Allow HTTP traffic" and "Allow HTTPS traffic"
5. Click double down arrows to expand options, click the "Networking" tab and add each of the new firewall rules by typing their names, "open-all-egress" and "open-all-ingress". A grey bubble should appear around the name of the new firewall rule
6. Click "Create"
7. Access the VM by clicking the down arrow next to the "SSH" button.
8. Install Docker
```
## Run these commands to install Docker and add non-root users to run Docker
sudo snap install docker
sudo apt update
sudo apt install -y docker.io
sudo usermod -aG docker $USER
exit
# exit and SSH back in for changes to take effect
```
9. [Install anaconda](https://docs.anaconda.com/anaconda/install/linux/)
```
sudo apt-get install -y libgl1-mesa-glx libegl1-mesa libxrandr2 libxrandr2 libxss1 libxcursor1 libxcomposite1 libasound2 libxi6 libxtst6

wget https://repo.anaconda.com/archive/Anaconda3-2019.03-Linux-x86_64.sh
bash Anaconda3-2019.03-Linux-x86_64.sh
## Press "Enter" a few times to view the Anaconda End User License Agreement
## Type "yes" to consent and begin installation.  Accept all default values.
exit
## Return to the VM by clicking "SSH" in the browser
```
10. Start the Jupyter Notebook server by running the following commands.
```
## The following Docker image is auto-built from
## https://github.com/stevetsa/jupyter-blast-docker

docker pull stevetsa/jupyter-blast-docker

## The following repository contains the Jupyter Notebook
## Change as necessary

git clone https://github.com/stevetsa/addon.git
chmod 777 addon
cd addon
docker run -it --rm -v `pwd`:`pwd` -w `pwd` -p 8888:8888 stevetsa/jupyter-blast-docker
```
Follow on-screen instructions and copy-and-paste the url with token in a brower's url field.  Then modify the URL so it has the form - http://\<external-ip-address\>:8888/?token=xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

Now you are ready to run BLAST from a Jupyter Notebook!  Click to open the notebook named `BLAST+ Docker Jupyter Notebook.ipynb` and run all command blocks in the notebook.

### Stop the GCP instance
Remember to [stop](https://cloud.google.com/compute/docs/instances/stop-start-instance) or [delete](https://cloud.google.com/compute/docs/instances/stop-start-instance) the VM to prevent incurring additional cost. You can do this at the GCP Console.

As an alternative, you can also run this notebook using [MyBinder.org](https://mybinder.readthedocs.io/en/latest/). You will be running the notebook in MyBinder's compute resources and you do not need a GCP account or access.  

To start using BLAST+ in a Jupyter Notebook, click [here](https://github.com/stevetsa/jupyter-blast-docker-binder) and click on the "launch binder" button. Storage and compute are provided by MyBinder.org; therefore, your analysis may be limited by the availability of computing resources.


