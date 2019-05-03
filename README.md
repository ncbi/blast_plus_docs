# Official NCBI BLAST+ Docker image

This repository contains documentation for the [NCBI BLAST+](http://blast.ncbi.nlm.nih.gov/) command line applications in a Docker image.     

*This is only a preview outline of the documentation.*
*Please continue to check this site for the latest updates!*

# What is NCBI BLAST?
This section includes the NCBI logo, an introduction to BLAST/BLAST+, and a description of the content of this repository.

# Get Started
This section includes a brief introduction of cloud computing and Docker, focusing only on relevant components used in this repository.

## What is Cloud Computing?
This section includes a brief introduction to cloud computing concepts and terms as well as reasons to migrate data and analysis to the cloud.
Reference to NIH/STRIDES material.
 
## What is Docker?
This section includes a brief introduction to containers, images and Docker, including a reference to Singularity.

## Cloud computing and Docker concept
This section uses a figure to explain various technological components, concepts and terms used in this repository.

## When is it appropriate to use the BLAST+ Docker image?
This section describes the appropriate usage for BLAST WebUI vs BLAST+ Docker image in terms of standard/custom databases and query size. 

## Setting up the Cloud Computing Environment
This section provides references (Cookbook, GCP, etc.) to set up the computing environment and provides performance benchmarking information to assist users in selecting the right computational resources.

## Software installation
This section covers Docker installation, examines the structure of a Docker command, and discusses command options used in this repository.
* Docker installation
* Docker command structure
* Docker command options

# How to use the BLAST+ Docker image?

This section covers the required inputs for BLAST+ and detailed step-by-step command line instructions for running BLAST in two steps.
+ Step 1. Setting up BLAST databases and query sequences 
    - Download and Install NCBI-provided BLAST databases
        - Show available BLAST databases on local host
        - Show BLAST databases available for download from GCP
        - Show BLAST databases available for download from NCBI
        - Download NCBI-provided BLAST databases
    - Make and install my own BLAST databases

+ Step 2. Run BLAST+ Docker image
This section covers various versions of BLAST in this Docker image, how to run the image in an interactive bash shell *inside a container* as well as scripted BLAST. 
    - Interactive BLAST
    - Scripted BLAST
 
# Advanced use case(s) of BLAST+ Docker image
This section highlights potential scientific application(s) of using BLAST+ Docker image, referencing examples from publications.

# Supported tags and respective Dockerfile links
This section includes tags, references and release notes for the versions of BLAST+ in this Docker image.

# Additional resources
This section covers additional resources to help user and NCBI/BLAST contact information.
* BLAST
* Docker
* GCP
* NIH/STRIDES
* Contact info
