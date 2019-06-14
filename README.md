# Official NCBI BLAST+ Docker Image

This repository contains documentation for the [NCBI BLAST+](https://www.ncbi.nlm.nih.gov/pubmed/2231712) command line applications in a Docker image.  We will demonstrate how to use the Docker image to run BLAST analysis on the Google Cloud Platform using a small basic example and a more advanced production-level example.  Some basic knowledge of Unix/Linux commands and BLAST+ is useful in completing this tutorial.  

## Table of Contents

   * [What Is NCBI BLAST?](#what-is-ncbi-blast)
   * [What Is Cloud Computing?](#what-is-cloud-computing)
   * [What Is Docker?](#what-is-docker)
   * [Section 1 - Getting Started Using the BLAST+ Docker Image with A Small Example](#section-1---getting-started-using-the-blast-docker-image-with-a-small-example)     
   * [Google Cloud Platform Setup](#google-cloud-platform-setup)  
   * [Section 2 - A Step-by-Step Guide Using the BLAST+ Docker Image](#section-2---a-step-by-step-guide-using-the-blast-docker-image)
       * [Step 1. Install Docker](#step-1-install-docker)
           * [Docker `run` command options](#docker-run-command-options)
           * [Docker `run` command structure](#docker-run-command-structure)
           * [Useful Docker commands](#useful-docker-commands)
           * [Using BLAST+ with Docker](#using-blast-with-docker)
           * [Versions of BLAST+ Docker image](#versions-of-blast-docker-image)
           * [Supported tags](#supported-tags)
       * [Step 2. Import sequences and create a BLAST database](#step-2---import-sequences-and-create-a-blast-database)
           * [Show BLAST databases available for download from the Google Cloud bucket](#show-blast-databases-available-for-download-from-the-google-cloud-bucket)
           * [Show BLAST databases available for download from NCBI](#show-blast-databases-available-for-download-from-ncbi)
           * [Show available BLAST databases on local host](#show-available-blast-databases-on-local-host)
       * [Step 3. Run BLAST](#step-3-run-blast)
       * [Stop the GCP instance](#stop-the-gcp-instance)
   * [Section 3 - Using the BLAST+ Docker Image at Production Scale](#section-3---using-the-blast-docker-image-at-production-scale)
       * [Background](#background)
       * [BLAST+ Docker image benchmarks](#blast-docker-image-benchmarks)
       * [Commands to run](#commands-to-run)
   * [Additional Resources](#additional-resources)
   * [Maintainer](#maintainer)
   * [License](#license)
   * [Appendix](#appendix)
       * [Appendix A. Cloud and Docker Concepts](#appendix-a-cloud-and-docker-concepts)
       * [Appendix B. Alternative Ways to Run Docker](#appendix-b-alternative-ways-to-run-docker)
       * [Appendix C. Transfer files to/from a GCP VM](#appendix-c-transfer-files-tofrom-a-gcp-vm)

# What Is NCBI BLAST?

<img align="right" width="150" height="200" src="https://www.nlm.nih.gov/about/logos_nlm_photos/large-White_ncbi_logo_200h.png" alt="ncbi logo">  


The National Center for Biotechnology Information (NCBI) Basic Local Alignment Search Tool [(BLAST)]( https://blast.ncbi.nlm.nih.gov) finds regions of local similarity between sequences. The program compares nucleotide or protein sequences to sequence databases and calculates the statistical significance of matches. BLAST can be used to infer functional and evolutionary relationships between sequences as well as help identify members of gene families.

Introduced in 2009, BLAST+ is an improved version of BLAST command line applications.  For a full description of the features and capabilities of BLAST+, please refer to the [BLAST Command Line Applications User Manual](https://www.ncbi.nlm.nih.gov/books/NBK279690/).

# What Is Cloud Computing?
Cloud computing offers potential cost savings by using on-demand, scalable, and elastic computational resources. While a detailed description of various cloud technologies and benefits is out of the scope for this repository, the following sections contain information needed to get started running the BLAST+ Docker image on the Google Cloud Platform [(GCP)]( https://cloud.google.com/).  

# What Is Docker?
[Docker](https://www.docker.com/) is a tool to perform operating-system level virtualization using software containers. In containerization technology<sup>*</sup>, an image is a snapshot of an analytical environment encapsulating application(s) and dependencies. An image, which is essentially a file built from a list of instructions, can be saved and easily shared for others to recreate the exact analytical environment across platforms and operating systems. A container is a runtime instance of an image. By using containerization, users can bypass the often-complicated steps in compiling, configuring, and installing a Unix-based tool like BLAST+. In addition to portability, containerization is a lightweight approach to make analysis more findable, accessible, interoperable, reusable (F.A.I.R.) and, ultimately, reproducible.  

*There are many containerization tools and standards, such as [Docker](https://www.docker.com/) and [Singularity]( https://www.sylabs.io/singularity/). We will focus solely on Docker, which is considered the de facto standard by many in the field.  
  
# Section 1 - Getting Started Using the BLAST+ Docker Image with a Small Example
   
This section provides a quick run-through of a BLAST analysis in the Docker environment on a Google instance. This is intended as an overview for those who just want an understanding of the principles of the solution.  The Google Cloud Shell, an interactive shell environment, will be used for this example, which makes it possible to run the following small example without having to perform additional setup, such as creating a billing account or compute instance.
More detailed descriptions of analysis steps, alternative commands, and more advanced topics are covered in the later sections of this documentation.  
  
Requirements:  A Google account
  
Flow of the Task:
![Task-Flow](images/task-flow.png)
  
Input data:
*	Query – 1 sequence, 44 nucleotides, file size 0.2 KB
*	Databases 
    *	7 sequences, 922 nucleotides, file size 1.7 KB
    *	PDB protein database (pdb_v5) 0.2945 GB  
    
First, in a separate browser window or tab, sign in at https://console.cloud.google.com/

Click the **Activate Cloud Shell** button at the top right corner of the Google Cloud Platform Console. 
![Activate-Cloud-Shell](images/activate-cloud-shell.png)
   
You now will see your Cloud Shell session window:
![Cloud-Shell-Commandline](images/cloud-shell-commandline.png)

The next step is to copy-and-paste the commands below in your Cloud Shell session.
  
*Please note: In GitHub you can use your mouse to copy; however, in the command shell you must use your keyboard. In Windows or Unix/Linux, use the shortcut `Control+C` to copy and `Control+V` to paste. On macOS, use `Command+C` to copy and `Command+V` to paste.*   

To scroll in the Cloud Shell, enable the scrollbar in `Terminal settings` with the wrench icon.
![Cloud-Shell-wrench](images/cloud-shell-wrench.png)

```
# Time needed to complete this section: <10 minutes

# Step 1. Retrieve sequences
## Create directories for analysis
cd ; mkdir blastdb queries fasta results blastdb_custom

## Retrieve query sequence
docker run --rm ncbi/blast efetch -db protein -format fasta \
    -id P01349 > queries/P01349.fsa
    
## Retrieve database sequences
docker run --rm ncbi/blast efetch -db protein -format fasta \
    -id Q90523,P80049,P83981,P83982,P83983,P83977,P83984,P83985,P27950 \
    > fasta/nurse-shark-proteins.fsa
    
## Step 2. Make BLAST database 
docker run --rm \
    -v $HOME/blastdb_custom:/blast/blastdb_custom:rw \
    -v $HOME/fasta:/blast/fasta:ro \
    -w /blast/blastdb_custom \
    ncbi/blast \
    makeblastdb -in /blast/fasta/nurse-shark-proteins.fsa -dbtype prot \
    -parse_seqids -out nurse-shark-proteins -title "Nurse shark proteins" \
    -taxid 7801 -blastdb_version 5
    
## Step 3. Run BLAST+ 
docker run --rm \
    -v $HOME/blastdb:/blast/blastdb:ro \
    -v $HOME/blastdb_custom:/blast/blastdb_custom:ro \
    -v $HOME/queries:/blast/queries:ro \
    -v $HOME/results:/blast/results:rw \
    ncbi/blast \
    blastp -query /blast/queries/P01349.fsa -db nurse-shark-proteins
    
## Output on screen
## Scroll up to see the entire output
## Type "exit" to leave the Cloud Shell or continue to the next section
```
  
At this point, you should see the output on the screen. With your query, BLAST identified the protein sequence P80049.1 as a match with a score of 14.2 and an E-value of 0.96.  
  
For larger analysis, it is recommended to use the `-out` flag to save the output to a file.  For example, append `-out /blast/results/blastp.out` to the last command in Step 3 above and view the content of this output file using `more $HOME/results/blastp.out`.  

You can also query P01349.fsa against the PDB as shown in the following code block.

```
## Extend the example to query against the Protein Data Bank
## Time needed to complete this section: <10 minutes

## Confirm query
ls queries/P01349.fsa

## Download Protein Data Bank Version 5 database (pdb_v5)
docker run --rm \
     -v $HOME/blastdb:/blast/blastdb:rw \
     -w /blast/blastdb \
     ncbi/blast \
     update_blastdb.pl --source gcp pdb_v5

## Run BLAST+ 
docker run --rm \
     -v $HOME/blastdb:/blast/blastdb:ro \
     -v $HOME/blastdb_custom:/blast/blastdb_custom:ro \
     -v $HOME/queries:/blast/queries:ro \
     -v $HOME/results:/blast/results:rw \
     ncbi/blast \
     blastp -query /blast/queries/P01349.fsa -db pdb_v5

## Output on screen
## Scroll up to see the entire output
## Leave the Cloud Shell

exit
```

You have now completed a simple task and seen how BLAST+ with Docker works. To learn about Docker and BLAST+ at production scale, please proceed to the next section.  
  

# Google Cloud Platform Setup
The following sections include instructions to create a Google virtual machine, install Docker, and run BLAST+ commands using the Docker image.  
  
In [Section 2 - A Step-by-Step Guide Using the BLAST+ Docker Image](#section-2---a-step-by-step-guide-using-the-blast-docker-image), we will use the same small example from the previous section and discuss alternative approaches, additional useful Docker and BLAST+ commands, and Docker command options and structures.  In [Section 3](#section-3---using-the-blast-docker-image-at-production-scale), we will demonstrate how to run the BLAST+ Docker image at production scale.  
  
First, you need to set up a Google Cloud Platform (GCP) virtual machine (VM) for analysis.  

## Requirements
* A GCP account linked to a billing account
* A GCP VM running Ubuntu 18.04LTS

## Set up your GCP account and create a VM for analysis
  
### 1. Creating your GCP account and registering for the free $300 credit program. (If you already have a GCP billing account, you can skip to step 2.)

* First, in a separate browser window or tab, sign in at https://console.cloud.google.com/
    * If you need to create one, go to https://cloud.google.com/ and click “Get started for free” to sign up for a trial account.
    * If you have multiple Google accounts, sign in using an Incognito Window (Chrome) or Private Window (Safari) or any other private browser window.
     
GCP is currently offering a $300 credit, which expires 12 months from activation, to incentivize new cloud users. The following steps will show you how to activate this credit.
You will be asked for billing information, but GCP will not auto-charge you once the trial ends; you must elect to manually upgrade to a paid account.
  
* After signing in, click **Activate** to activate the $300 credit.
![GCP credit](images/gcp-credit.png)

* Enter your country, for example, **United States,** and check the box indicating that you have read and accept the terms of service.
* Under “Account type,” select “Individual.” (This may be pre-selected in your Google account)
* Enter your name and address.
* Under “How you pay," select “Automatic payments.” (This may be pre-selected in your Google account)
This indicates that you will pay costs after you have used the service, either when you have reached your billing threshold or every 30 days, whichever comes first.  
  
* Under “Payment method,” select “add a credit or debit card” and enter your credit card information.
You will not be automatically charged once the trial ends. You must elect to upgrade to a paid account before your payment method will be charged.  
  
* Click “Start my free trial” to finish registration. When this process is completed, you should see a GCP welcome screen. 
  
### 2. Create a Virtual Machine (VM)  
* On the GCP welcome screen from the last step, click "Compute Engine" or navigate to the "Compute Engine" section by clicking on the navigation menu with the "hamburger icon" (three horizontal lines) on the top left corner.  
<img align="left " width="300" src="images/gcp-instance.png" alt="GCP instance">   
  
* Click on the blue “CREATE INSTANCE” button on the top bar.  
* Create an image with the following parameters: (if parameter is not list below, keep the default setting)
    * Name: keep the default or enter a name
    * Region: **us-east4 (Northern Virginia)**   
    * For Section 2, change these settings - 
        * Machine Type: **micro (1 shared vCPU), 0.6 GB memory, f1-micro**
        * Boot Disk: Click "Change," select **Ubuntu 18.04 LTS,** and click "Select" (Boot disc size is default 10 GB).
    * For Section 3, change these settings - 
        * Machine Type: **16 vCPU, 104 GB memory, n1-highmem-16**
        * Boot Disk: Click "Change" and select **Ubuntu 18.04 LTS**, change the "Boot disk size" to **200 GB** Standard persistent disk, and click "Select." 

At this point, you should see a cost estimate for this instance on the right side of your window.  
![GCP VM cost](images/gcp-vm-cost1.png)  
   
* Click the blue “Create” button. This will create and start the VM.  
  
*Please note: Creating a VM in the same region as storage can provide better performance. We recommend creating a VM in the us-east4 region. If you have a job that will take several hours, but less than 24 hours, you can potentially take advantage of [preemptible VMs.](https://cloud.google.com/compute/docs/instances/preemptible)*
  
Detailed instructions for creating a GCP account and launching a VM can be found [here.](https://cloud.google.com/compute/docs/quickstart-linux)
   
### 3. Access a GCP VM from a local machine

Once you have your VM created, you must access it from your local computer. There are many methods to access your VM, depending on the ways in which you would like to use it. On the GCP, the most straightforward way is to SSH from the browser.

* Connect to your new VM instance by clicking the "SSH" button
![GCP SSH](images/gcp-ssh.png)

You now have a command shell running and you are ready to proceed.

Remember to [stop](https://cloud.google.com/compute/docs/instances/stop-start-instance) or [delete](https://cloud.google.com/compute/docs/instances/stop-start-instance) the VM to prevent incurring additional cost.    
    
    
# Section 2 - A Step-by-Step Guide Using the BLAST+ Docker Image
In this section, we will cover Docker installation, discuss various `docker run` command options, and examine the structure of a Docker command.  We will use the same small example from Section 1 and explore alternative approaches in running the BLAST+ Docker image. However, we are using a real VM instance, which provides greater performance and functionality than the Google Cloud Shell.

Input data
* Query – 1 sequence, 44 nucleotides, file size 0.2 KB
* Database – 7 sequences, 922 nucleotides, file size 1.7 KB
   
## Step 1. Install Docker
In a production system, Docker has to be installed as an application.  

```
## Run these commands to install Docker and add non-root users to run Docker
sudo snap install docker
sudo apt update
sudo apt install -y docker.io
sudo usermod -aG docker $USER
exit
# exit and SSH back in for changes to take effect
```
To confirm the correct installation of Docker, run the command `docker run hello-world`. If correctly installed, you should see "Hello from Docker!..."(https://docs.docker.com/samples/library/hello-world/)  
  
### Docker run command options
*This section is optional.*    
  
Below is a list of `docker run` command line [options](https://docs.docker.com/engine/reference/commandline/run/) used in this tutorial.

| Name, short-hand(if available) | Description |
| :----------------------------  | :---------- |
|`--rm`|Automatically remove the container when it exits|
|`--volume` , `-v`|Bind mount a volume|
|`--workdir` , `-w`| Working directory inside the container|

### Docker run command structure
*This section is optional.*  
 
For this tutorial, it would be useful to understand the structure of a Docker command. The following command consists of three parts.

```
docker run --rm ncbi/blast \
```
```
    -v $HOME/blastdb_custom:/blast/blastdb_custom:rw \
    -v $HOME/fasta:/blast/fasta:ro \
    -w /blast/blastdb_custom \
```
```
    makeblastdb -in /blast/fasta/nurse-shark-proteins.fsa -dbtype prot \
    -parse_seqids -out nurse-shark-proteins -title "Nurse shark proteins" \
    -taxid 7801 -blastdb_version 5
```
  
The first part of the command `docker run --rm ncbi/blast` is an instruction to run the docker image `ncbi/blast` and remove the container when the run is completed.  
  
The second part of the command makes the query sequence data accessible in the container. [Docker bind mounts]( https://docs.docker.com/storage/bind-mounts/) uses `-v` to mount the local directories to directories inside the container and provide access permission rw (read and write) or ro (read only). For instance, assuming your subject sequences are stored in the $HOME/fasta directory on the local host, you can use the following parameter to make that directory accessible inside the container in /blast/fasta as a read-only directory `-v $HOME/fasta:/blast/fasta:ro`.  The `-w /blast/blastdb_custom` flag sets the working directory inside the container.
  
The third part of the command is the BLAST+ command. In this case, it is executing makeblastdb to create BLAST database files.  
  
You can start an interactive bash session for this image by using `docker run -it ncbi/blast /bin/bash`. For the BLAST+ Docker image, the executables are in the folder /blast/bin and /root/edirect and added to the variable $PATH.  
  
For additional documentation on the `docker run` command, please refer to [documentation](https://docs.docker.com/engine/reference/commandline/run/).  
  
### Useful Docker commands
*This section is optional.*  
  
| Docker Command | Description |
| :----------------------------  | :---------- |
|`docker ps -a`|Displays a list of containers|
|`docker rm $(docker ps -q -f status=exited)`|Removes all exited containers, if you have at least 1 exited container|
|`docker rm <CONTAINER_ID>`|Removes a container|
|`docker images`|Displays a list of images|
|`docker rmi <REPOSITORY (IMAGE_NAME)>`|Removes an image|
  
### Using BLAST+ with Docker
*This section is optional.*    
  
With this Docker image you can run BLAST+ in an isolated container, facilitating reproducibility of BLAST results. As a user of this Docker image, you are expected to provide BLAST databases and query sequence(s) to run BLAST
as well as a location outside the container to save the results. The following is a list of directories used by BLAST+. You will create them in Step 2.  
  
| Directory | Purpose | Notes |
| --------- | ------  | ----- |
| `$HOME/blastdb` | Stores NCBI-provided BLAST databases | If set to a _single, absolute_ path, the `$BLASTDB` environment variable could be used instead (see [Configuring BLAST via environment variables](https://www.ncbi.nlm.nih.gov/books/NBK279695/#_usermanual_Configuring_BLAST_via_environ_).) |
| `$HOME/queries` | Stores user-provided query sequence(s) | |
| `$HOME/fasta`   | Stores user-provided FASTA sequences to create BLAST database(s) | |
| `$HOME/results` | Stores BLAST results | Mount with `rw` permissions |
| `$HOME/blastdb_custom` | Stores user-provided BLAST databases | |

### Versions of BLAST Docker image
*This section is optional.*  
  
The following command displays the latest BLAST version.  
```docker run --rm ncbi/blast blastn -version```

Appending a tag to the image name (`ncbi/blast`) allows you to use a
different version of BLAST+ (see “Supported Tags and Respective Release Notes” section for supported versions).  

Different versions of BLAST+ exist in different Docker images. The following command will initiate download of the BLAST+ version 2.7.1 Docker image. 
```
docker run --rm ncbi/blast:2.7.1 blastn -version
## Display a list of images
docker images
```

For example, to use the BLAST+ version 2.7.1 Docker image instead of the latest version, replace the first part of the command

```docker run --rm ncbi/blast``` with ```docker run --rm ncbi/blast:2.7.1 ```

### Supported tags
*This section is optional.*   
  
* [2.9.0](https://github.com/ncbi/docker/blob/master/blast/2.9.0/Dockerfile): [release notes](https://www.ncbi.nlm.nih.gov/books/NBK131777/#_Blast_ReleaseNotes_BLAST_2_9_0_April_01)
* [2.8.1](https://github.com/ncbi/docker/blob/master/blast/2.8.1/Dockerfile): [release notes](https://www.ncbi.nlm.nih.gov/books/NBK131777/#_Blast_ReleaseNotes_BLAST_2_8_1_DECEMBER_1_)
* [2.8.0](https://github.com/ncbi/docker/blob/master/blast/2.8.0/Dockerfile): [release notes](https://www.ncbi.nlm.nih.gov/books/NBK131777/#_Blast_ReleaseNotes_BLAST_2_8_0_March_28_)
* [2.7.1](https://github.com/ncbi/docker/blob/master/blast/2.7.1/Dockerfile): [release notes](https://www.ncbi.nlm.nih.gov/books/NBK131777/#_Blast_ReleaseNotes_BLAST_2_7_1_October_2_)

## Step 2. Import sequences and create a BLAST database
In this example, we will start by fetching query and database sequences and then create a custom BLAST database.  

```
# Start in a directory where you want to perform the analysis
## Create directories for analysis
cd ; mkdir blastdb queries fasta results blastdb_custom

## Retrieve query sequences
docker run --rm ncbi/blast efetch -db protein -format fasta \
    -id P01349 > queries/P01349.fsa
    
## Retrieve database sequences
docker run --rm ncbi/blast efetch -db protein -format fasta \
    -id Q90523,P80049,P83981,P83982,P83983,P83977,P83984,P83985,P27950 \
    > fasta/nurse-shark-proteins.fsa
    
## Make BLAST database 
docker run --rm \
    -v $HOME/blastdb_custom:/blast/blastdb_custom:rw \
    -v $HOME/fasta:/blast/fasta:ro \
    -w /blast/blastdb_custom \
    ncbi/blast \
    makeblastdb -in /blast/fasta/nurse-shark-proteins.fsa -dbtype prot \
    -parse_seqids -out nurse-shark-proteins -title "Nurse shark proteins" \
    -taxid 7801 -blastdb_version 5
```

To verify the newly created BLAST database above, you can run the following command to display the accessions, sequence length, and common name of the sequences in the database.

```
docker run --rm \
    -v $HOME/blastdb:/blast/blastdb:ro \
    -v $HOME/blastdb_custom:/blast/blastdb_custom:ro \
    ncbi/blast \
    blastdbcmd -entry all -db nurse-shark-proteins -outfmt "%a %l %T"
```
  
As an alternative, you can also download preformatted BLAST databases from NCBI or the NCBI Google storage bucket.  
   
### Show BLAST databases available for download from the Google Cloud bucket  
  
```docker run --rm ncbi/blast update_blastdb.pl --showall pretty --source gcp```

For a detailed description of `update_blastdb.pl`, please refer to the [documentation.](https://www.ncbi.nlm.nih.gov/books/NBK537770/)
  
### Show BLAST databases available for download from NCBI  
*This section is optional.*  
  
```docker run --rm ncbi/blast update_blastdb.pl --showall --source ncbi```  
  
### Show available BLAST databases on local host 
*This section is optional.*   
  
The command below mounts the `$HOME/blastdb` path on the local machine as
`/blast/blastdb` on the container, and `blastdbcmd` shows the available BLAST
databases at this location.  
  
```
## Download Protein Data Bank Version 5 database (pdb_v5)
docker run --rm \
     -v $HOME/blastdb:/blast/blastdb:rw \
     -w /blast/blastdb \
     ncbi/blast \
     update_blastdb.pl --source gcp pdb_v5

## Display database(s) in $HOME/blastdb
docker run --rm \
    -v $HOME/blastdb:/blast/blastdb:ro \
    ncbi/blast \
    blastdbcmd -list /blast/blastdb -remove_redundant_dbs
```
  
You should see an output `/blast/blastdb/pdb_v5 Protein`.  
  
```
## For the custom BLAST database used in this example -
docker run --rm \
    -v $HOME/blastdb_custom:/blast/blastdb_custom:ro \
    ncbi/blast \
    blastdbcmd -list /blast/blastdb_custom -remove_redundant_dbs
```
You should see an output `/blast/blastdb_custom/nurse-shark-proteins Protein`.  
   
## Step 3. Run BLAST
When running BLAST in a Docker container, note the mounts specified to the `docker run` command to make the input and outputs accessible. In the examples below, the first two mounts provide access to the BLAST databases, the third
mount provides access to the query sequence(s), and the fourth mount provides a directory to save the results. (Note the `:ro` and `:rw` options, which mount the directories as read-only and read-write respectively.)  
  
```
docker run --rm \
    -v $HOME/blastdb:/blast/blastdb:ro \
    -v $HOME/blastdb_custom:/blast/blastdb_custom:ro \
    -v $HOME/queries:/blast/queries:ro \
    -v $HOME/results:/blast/results:rw \
    ncbi/blast \
    blastp -query /blast/queries/P01349.fsa -db nurse-shark-proteins \
    -out /blast/results/blastp.out
```  
At this point, you should see the output file ```$HOME/results/blastp.out```. With your query, BLAST identified the protein sequence P80049.1 as a match with a score of 14.2 and an E-value of 0.96. To view the content of this output file, use the command ```more $HOME/results/blastp.out```.    

## Stop the GCP instance
Remember to [stop](https://cloud.google.com/compute/docs/instances/stop-start-instance) or [delete](https://cloud.google.com/compute/docs/instances/stop-start-instance) the VM to prevent incurring additional cost. You can do this at the GCP Console as shown below.
![GCP instance stop](images/gcp-instance-stop.png)
  
# Section 3 - Using the BLAST+ Docker Image at Production Scale
## Background
One of the promises of cloud computing is scalability. In this section, we will demonstrate how to use the BLAST+ Docker image at production scale on the Google Cloud Platform. We will perform a BLAST analysis similar to the approach described in this [publication](https://www.ncbi.nlm.nih.gov/pubmed/31040829) to compare de novo aligned contigs from bacterial 16S-23S sequencing against the nucleotide collection (nt) database.

To test scalability, we will use inputs of different sizes to estimate the amount of time to download the nucleotide collection database and run BLAST search using the latest version of the BLAST+ Docker image. Expected results are summarized in the following tables.

Input files: 28 samples (multi-FASTA files) containing de novo aligned contigs from the publication.  
(Instructions to [download]((https://figshare.com/s/729b346eda670e9daba4)) and create the input files are described in the [code block](#commands-to-run) below.)    
  
Database: Pre-formatted BLAST nucleotide collection database, version 5 (nt_v5): 68.7217 GB  
  
|       | Input file name | File content | File size | Number of sequences | Number of nucleotides | Expected output size |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| Analysis 1 | query1.fa | only sample 1 | 59 KB | 121 | 51,119 | 3.1 GB |
| Analysis 2 | query5.fa | only samples 1-5 | 422 KB | 717 | 375,154 | 10.4 GB |
| Analysis 3 | query.fa | all 28 samples | 2.322 MB | 3798 | 2,069,892 | 47.8 GB |


## BLAST+ Docker image benchmarks  
| VM Type/Zone | CPU | Memory (GB) | Hourly Cost* | Download nt (min) | Analysis 1 (min) | Analysis 2 (min) | Analysis 3 (min)| Total Cost**
| :-: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :--: |
| n1-standard-8 us-east4c | 8 | 30 | $0.312 | 9 | 22 | - | - | - |
| n1-standard-16 us-east4c | 16 | 60 | $0.611 | 9 | 14 | 53 | 205 | $2.86 |
| n1-highmem-16 us-east4c | 16 | 104 | $0.767 | 9 | 9 | 30 | 143 | $2.44 |
| n1-highmem-16 us-west2a  | 16 | 104 | $0.809 | 11 | 9 | 30 | 147 | $2.60 |
| n1-highmem-16 us-west1b | 16 | 104 | $0.674 | 11 | 9 | 30 | 147 | $2.17 |
| BLAST website (blastn) | - | - | - | - |Searches exceed current restrictions on usage|Searches exceed current restrictions on usage|Searches exceed current restrictions on usage| - |


All GCP instances are configured with 200 GB of persistent standard disk. 
  

*Hourly costs were provided by Google Cloud Platform (May 2019) when VMs were created and are subject to change.  
**Total costs were estimated using the hourly cost and total time to download nt and run Analysis 1, Analysis 2, and Analysis 3. Estimates are used for comparison only; your costs may vary and are your responsibility to monitor and manage.
   
Please refer to GCP for more information on [machine types](https://cloud.google.com/compute/docs/machine-types),
[regions and zones,](https://cloud.google.com/compute/docs/regions-zones/) and [compute cost.](https://cloud.google.com/compute/pricing)


## Commands to run
```
## Install Docker if not already done
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

## Step 2. Display BLAST databases on the GCP
docker run --rm ncbi/blast update_blastdb.pl --showall pretty --source gcp

## Download nt_v5 (nucleotide collection version 5) database
## This step takes approximately 10 min.  The following command runs in the background.
docker run --rm \
  -v $HOME/blastdb:/blast/blastdb:rw \
  -w /blast/blastdb \
  ncbi/blast \
  update_blastdb.pl --source gcp nt_v5 &

## At this point, confirm query/database have been properly provisioned before proceeding

## Check the size of the directory containing the BLAST database
## nt_v5 should be around 68 GB
du -sk $HOME/blastdb

## Check for queries, there should be three files - query.fa, query1.fa and query5.fa
ls -al $HOME/queries

## From this point forward, it may be easier if you run these steps in a script. 
## Simply copy and paste all the commands below into a file named script.sh
## Then run the script in the background `nohup bash script.sh > script.out &`

## Step 3. Run BLAST
## Run BLAST using query1.fa (Sample 1) 
## This command will take approximately 9 minutes to complete.
## Expected output size: 3.1 GB  
docker run --rm \
  -v $HOME/blastdb:/blast/blastdb:ro -v $HOME/blastdb_custom:/blast/blastdb_custom:ro \
  -v $HOME/queries:/blast/queries:ro \
  -v $HOME/results:/blast/results:rw \
  ncbi/blast \
  blastn -query /blast/queries/query1.fa -db nt_v5 -num_threads 16 \
  -out /blast/results/blastn.query1.denovo16s.out

## Run BLAST using query5.fa (Samples 1-5) 
## This command will take approximately 30 minutes to complete.
## Expected output size: 10.4 GB  
docker run --rm \
  -v $HOME/blastdb:/blast/blastdb:ro -v $HOME/blastdb_custom:/blast/blastdb_custom:ro \
  -v $HOME/queries:/blast/queries:ro \
  -v $HOME/results:/blast/results:rw \
  ncbi/blast \
  blastn -query /blast/queries/query5.fa -db nt_v5 -num_threads 16 \
  -out /blast/results/blastn.query5.denovo16s.out

## Run BLAST using query.fa (All 28 samples) 
## This command will take approximately 147 minutes to complete.
## Expected output size: 47.8 GB  
docker run --rm \
  -v $HOME/blastdb:/blast/blastdb:ro -v $HOME/blastdb_custom:/blast/blastdb_custom:ro \
  -v $HOME/queries:/blast/queries:ro \
  -v $HOME/results:/blast/results:rw \
  ncbi/blast \
  blastn -query /blast/queries/query.fa -db nt_v5 -num_threads 16 \
  -out /blast/results/blastn.query.denovo16s.out

## Stdout and stderr will be in script.out
## BLAST output will be in $HOME/results
``` 

You have completed the entire tutorial.  At this point, if you do not need the downloaded data for further analysis, please [delete](https://cloud.google.com/compute/docs/instances/deleting-instance) the VM to prevent incurring additional cost.  
  
To delete an instance, follow instructions in the section [Stop the GCP instance.](#stop-the-gcp-instance)
    
For additional information, please refer to Google Cloud Platform's documentation on [instance life cycle.](https://cloud.google.com/compute/docs/instances/instance-life-cycle)
  
# Additional Resources
* BLAST:
    * [BLAST Command Line Applications User Manual](https://www.ncbi.nlm.nih.gov/books/NBK279696/)  
    * [BLAST Knowledge Base](https://support.nlm.nih.gov/knowledgebase/category/?id=CAT-01239)
* Docker: 
    * [Docker Community Forums](https://forums.docker.com)
    * [Docker Community Slack](https://blog.docker.com/2016/11/introducing-docker-community-directory-docker-community-slack/)
    * [Stack Overflow](https://stackoverflow.com/search?tab=newest&q=docker+blast)
* Other:
    * [Common Workflow Language (CWL)](https://www.commonwl.org/) is a specification to describe tools and workflows.  This [GitHub Repository](https://github.com/ncbi/cwl-demos/tree/master/blast-pipelines) contains sample CWL workflows using containerized BLAST+.
    * [Google Cloud Platform](https://cloud.google.com/)
    * [NIH/STRIDES](https://datascience.nih.gov/strides)
    * [GitHub](https://github.com/ncbi)
    
or [email us.](mailto:blast-help@ncbi.nlm.nih.gov)

# Maintainer

[National Center for Biotechnology Information (NCBI)](https://www.ncbi.nlm.nih.gov/)  
[National Library of Medicine (NLM)](https://www.nlm.nih.gov/)  
[National Institutes of Health (NIH)](https://www.nih.gov/)

# License

View refer to the [license](https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/scripts/projects/blast/LICENSE) and [copyright](http://ncbi.github.io/blast-cloud/dev/copyright.html) information for the software contained in this image.

As with all Docker images, these likely also contain other software which may be under other licenses (such as bash, etc., from the base distribution, along with any direct or indirect dependencies of the primary software being contained).

As with any pre-built image usage, it is the image user's responsibility to ensure that any use of this image complies with any relevant licenses for all software contained within.

# Appendix
## Appendix A. Cloud and Docker Concepts
![Cloud-Docker-Simple](images/cloud-docker-simple.png)
Figure 1. Docker and Cloud Computing Concept. Users can access compute resources provided by cloud service providers (CSPs), such as the Google Cloud Platform, using SSH tunneling (1). When you create a VM (2), a hard disk (also called a boot/persistent disk) (3) is attached to that VM. With the right permissions, VMs can also access other storage buckets (4) or other data repositories in the public domain. Once inside a VM with Docker installed, you can run a Docker image (5), such as NCBI's BLAST image. An image can be used to create multiple running instances or containers (6). Each container is in an isolated environment. In order to make data accessible inside the container, you need to use Docker bind mounts (7) described in this tutorial. 

*A Docker image can be used to create a Singularity image.  Please refer to Singularity's [documentation](https://www.sylabs.io/singularity/) for more detail.*

## Appendix B. Alternative Ways to Run Docker
  
As an alternative to what is described above, you can also run BLAST interactively inside a container.   


### Run BLAST+ Docker image interactively  
__When to use__: This is useful for running a few (e.g., fewer than 5-10) BLAST searches on small BLAST databases where you expect the search to complete in seconds/minutes.  
  
```
docker run --rm -it \
    -v $HOME/blastdb:/blast/blastdb:ro -v $HOME/blastdb_custom:/blast/blastdb_custom:ro \
    -v $HOME/queries:/blast/queries:ro \
    -v $HOME/results:/blast/results:rw \
    ncbi/blast \
    /bin/bash

# Once you are inside the container (note the root prompt), run the following BLAST commands.
blastp -query /blast/queries/P01349.fsa -db nurse-shark-proteins \
    -out /blast/results/blastp.out

# To view output, run the following command
more /blast/results/blastp.out

# Leave container
exit

```

In addition, you can run BLAST in [detached mode](https://docs.docker.com/engine/reference/run/#detached--d) by running a container in the background.  

### Run BLAST+ Docker image in detached mode 

__When to use__: This is a more practical approach if you have many (e.g., 10 or
more) BLAST searches to run or you expect the search to take a long time to execute. In this case it may be better to start the BLAST container in detached mode and execute commands on it. 
  
**NOTE**: Be sure to mount _all_ required directories, as these need to be
specified when the container is started.

```
# Start a container named 'blast' in detached mode
docker run --rm -dit --name blast \
    -v $HOME/blastdb:/blast/blastdb:ro -v $HOME/blastdb_custom:/blast/blastdb_custom:ro \
    -v $HOME/queries:/blast/queries:ro \
    -v $HOME/results:/blast/results:rw \
    ncbi/blast \
    sleep infinity

# Check the container is running in the background
docker ps -a
docker ps --filter "status=running"
```
Once the container is confirmed to be [running in detached mode](https://docs.docker.com/engine/reference/commandline/ps/), run the following BLAST command.
  
```
docker exec blast blastp -query /blast/queries/P01349.fsa \
    -db nurse-shark-proteins -out /blast/results/blastp.out

# View output
more $HOME/results/blastp.out

# stop the container
docker stop blast
```
If you run into issues with `docker stop blast` command, reset the VM from the GCP Console or restart the SSH session.  

## Appendix C. Transfer Files to/from a GCP VM

To copy the file `$HOME/script.out` in the home directory on a local machine to the home directory on a GCP VM named `instance-1` in project `My First Project` using GCP Cloud SDK.

GCP [documentation](https://cloud.google.com/compute/docs/instances/transfer-files)

First install GCP [Cloud SDK]( https://cloud.google.com/sdk/) command line tools for your operating system.

```
# First, set up gcloud tools
# From local machine's terminal

gcloud init

# Enter a configuration name
# Select the sign-in email account
# Select a project, for example “my-first-project”
# Select a compute engine zone, for example, “us-east4-c”

# To copy the file $HOME/script.out to the home directory of GCP instance-1 
# Instance name can be found in your Google Cloud Console -> Compute Engine -> VM instances

gcloud compute scp $HOME/script.out instance-1:~

# Optional - to transfer the file from the GCP instance to a local machine's home directory

gcloud compute scp instance-1:~/script.out $HOME/.
```
