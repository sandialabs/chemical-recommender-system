<div align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="src/App/static/assets/CRSwhite.png">
    <img alt="Text changing depending on mode. Light: 'So light!' Dark: 'So dark!'" src="src/App/static/assets/CRS.png" width="300px" align="center">
  </picture>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="src/App/static/assets/SNL.png">
    <img alt="Text changing depending on mode. Light: 'So light!' Dark: 'So dark!'" src="src/App/static/assets/SNL_Stacked.png" width="265px" align="center">
  </picture>
</div>

# Chemical Recommender System (CRS)

The Chemical Recommender System (CRS) is an advanced tool designed to assist researchers in identifying and comparing chemical compounds based on various criteria such as structural similarity, thermophysical properties, and toxicity. Utilizing state-of-the-art machine learning models and vector databases, CRS streamlines the process of chemical discovery and evaluation, making it an invaluable resource for scientific research and development. The CRS is elastic to user needs, transparent in its methods, and robust in allowing users to incorporate their own models and comparison metrics within CRS runs. Explore the documentation below to use the CRS.

## Instructions to Run

Warning: The current release of the CRS does not work on ARM64 machines, including M-series Macbooks. To run the Chemical Recommender System, first ensure that your machine has docker installed and configured (https://docs.docker.com/engine/install/). Then make sure the docker engine is running and download the docker-compose.yml file from this repo. Store it in a new directory, enter the directory. In here, create a .env file and create your settings as described below in the documentation. This is a necessary step to ensure all functionality. Finally, and then use the following command:

`docker compose up -d`

This will automatically open the web version of the application in your browser at `localhost:5005`. You will open up on the home page and the functionality of the CRS is found at the search and batch pages. When you are done using the CRS, you can shut down all containers and networks with the command `docker compose down`. If you are running into issues, make sure installation has occured correctly by deleting all volumes once the images are pulled and rerunning the compose up command.

## Command Line Interface (CLI)

The easiest way to use the CRS is described above through the webapp. However, using the CLI, developers can integrate their own machine learning models or comparison metrics they have created. To do this, enter the container shell and use the CLI. Continue from the previous steps and then:

`docker exec -it CRS /bin/bash`

The created terminal in a linux command line for the container and can generally be used as a normal linux machine. A variety of ways to use the CRS begin from here and can be further explore below at [Command Line / Developer Integrations](#command-line--developer-integrations). Each run will produce an report in the form of a and csv in the container.

When finished using the CRS, copy resultant output files from the Docker container to your local machine, use the following command. You will want to copy over the output directory recursively:

`docker cp CRS:/app/output /local/path/to/output`

## Table of Contents

- [Chemical Recommender System (CRS)](#chemical-recommender-system-crs)
  - [Instructions to Run](#instructions-to-run)
  - [Command Line Interface (CLI)](#command-line-interface-cli)
  - [Table of Contents](#table-of-contents)
  - [Command Line / Developer Integrations](#command-line--developer-integrations)
    - [Input File Format](#input-file-format)
  - [Environment Variables:](#environment-variables)
    - [Parameter Descriptions:](#parameter-descriptions)
    - [Example Inputs:](#example-inputs)
    - [Running the CRS](#running-the-crs)
      - [Example Commands:](#example-commands)
  - [Preprocessing](#preprocessing)
    - [PubChem](#pubchem)
  - [Running Fingerprint Similarity](#running-fingerprint-similarity)
  - [Further Sorting](#further-sorting)
    - [Thermophysical Comparison](#thermophysical-comparison)
    - [Toxicity Evaluation](#toxicity-evaluation)
    - [Structural Similarity](#structural-similarity)
    - [Synthetic Accessibility Scoring](#synthetic-accessibility-scoring)
  - [Technical Details](#technical-details)
    - [Vector Database Solution](#vector-database-solution)
    - [OPERA Integration](#opera-integration)
    - [RDKit SA Scoring](#rdkit-sa-scoring)
  - [Integrating Your Own Models](#integrating-your-own-models)
    - [Steps to Integrate Your Model](#steps-to-integrate-your-model)

## Command Line / Developer Integrations

Now, run the CRS by using the command `python src/main.py (args)`. For a list of options, run `python src/main.py -h` to get help. Below is a basic overview of the options.

### Input File Format

Input search parameters for the run must be written into a separate file. You can use text editors like nano or vim. The format of the input should be as follows:

```sh
query, final_number, thermo_array, include_all_elements, include_specific_elements, substructure_search, number_substructure_search
```

## Environment Variables:

Create a .env file in the directory of your `docker-compose.yml`. If you have the need for proxies, these are specified here. Also, if the default CRS port, 5005, is in use, you can change this to whatever may be open on your machine. Finally, ensure you use the correct architecture for your image pull. The options for TARGETARCH are amd64 and arm64. If you are using a Windows machine backed by WSL, you are likely using amd64. A Silicon M-Series Macboook will use arm64. If you add no env file, the default will run with no proxies on port 5005, assuming an amd64 machine. Update the settings in the file as so:

```sh
HTTP_PROXY=http://your-proxy-server:port
HTTPS_PROXY=https://your-proxy-server:port
CRS_PORT=xxxx
TARGETARCH=xxx64
```

#### Parameter Descriptions:

- **query**: The query to be searched for. This can be a PubChem CID, IUPAC Name, or SMILES.
- **final_number**: The number of resultant candidates to be outputted in the final report.
- **thermo_array**: An array of 5 boolean values (True or False) to decide if the user wants to include the following thermophysical properties (in this order, no spaces between the commas):
  - Melting Point
  - Boiling Point
  - Log P
  - Vapor Pressure
  - Henry's Law Constant
- **include_all_elements**: CRS by default only searches these elements: H, C, N, O, F, P, S, Cl, Se, Br, I. Setting this to parameter to True will include all elements instead. (Must be True or False)
- **include_specific_elements**: Add specific elements to search in addition to the default ones. Should be in a comma-separated format with no spaces in between. Must be either this or None.
- **substructure_search**: Provide a SMARTS representation of a substructure to require in all candidates. If not using, must leave as None.
- **number_substructure_search**: Signifies how many occurrences of the substructure must appear in the candidate. If not a number, leave as None. If a substructure is given and this is left as None, the search will look for at least 1 or more occurrences.
- **weights (optional)**: An array of weights signifying how to weigh each comparison value in the total rankings. This is [1,1,1,1,1] by default. These correspond to:
  - Structural Similarity
  - Molecular Weight Similarity
  - Thermophysical Similarity
  - Evaluated Toxicity
  - Synthetic Accesibility Scoring

#### Example Inputs:

```sh
6517, 30, [True,True,False,False,False], False, [Si], CCO, 1
```

This example tells the CRS to search for recommendations for PubChem CID 6517. It would create a report with 30 candidates and use the
thermophysical properties of Melting Point and Boiling Point in
comparison. The elements allowed in candidates are the listed default
plus Silicon. Finally, all returned candidates will have at least one
occurrence of the SMARTS 'CCO'.

```sh
quinolin-8-ol, 30, [False,False,False,False,True], True, None, None, None, [2,1,1,1,0]
```

This example tells the CRS to search for recommendations for the IUPAC Name quinolin-8-ol. It would create a report with 10 candidates and use the thermophysical property of the Henry's Law Constant. Candidates will be allowed to have any elements in it. The final sorting will disregard SA Scoring in the calculation and give increased weightage to structural similarity.

### Running the CRS

Tell the CRS to use this file as input by using the `-i` or `--input` argument followed by the file path. The only other required argument is to specify where the resultant PDF/CSV of the run should be stored. This will go into the output directory by default. Specify a location in this directory by using `-o` or `--output` followed by a file name or path. Another option is to use `-t` or `--time` instead of `-o` to automatically name the PDF/CSV by the date and time when the search began.

#### Example Commands:

1. Using a specified output file name:
   ```sh
   python src/main.py -i input.txt -o report
   ```
2. Automatically naming the output file by date and time:
   ```sh
    python src/main.py -i input.txt -t
   ```
   You can then add your own model and comparison metrics through the command line with the -m or --model option followed by a list of model names. For further explanation on this topic, refer to [Steps to Integrate Your Model](#steps-to-integrate-your-model).

In the case that the webapp of the CRS goes down, it can be restarted with the -w or --webapp option in the following format: `python src/main.py -w`.

## Preprocessing

The content of this section and its corresponding code is already ran by the developers and the data is stored inside the generated Docker volume. This code does not get run again during use of the CRS, but is here for the user's reference in how the CRS was created. The molecular fingerprint is implemented as a 2048-bit long vector. Each bit is set to on/off to represent the existence of a certain structural property. Comparison of these allows for a quick and computationally efficient way of comparing molecules. To do so, we must first create a database of all possible candidates and their respective fingerprints.

### PubChem

This process is done in the `src/Preprocessing` folder in the `SMILEStoFP.py` file. The input file path should be changed to a `.txt` file of PubChem CIDs and SMILES (can be downloaded at [PubChem FTP](https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/)). The program outputs and saves a CSV of each CID and its molecular fingerprint using RDKit. This process can be very lengthy and should be run using HPC. Next, the fingerprints are stored in the vector database which uses an indexing operation for fast searching upon upload. The code to do so with Milvus is shown in `src/Preprocessing/MilvusAdd.py`.

This step is already performed and integrated into the program, so the user will not have to run this on each query.

## Running Fingerprint Similarity

On the user's end, the following steps include going onto the search page of the website and inputting the correct query. This search is done using a vector databasing service, Milvus, which is expanded upon in later sections.

The result of the query is what follows. The following steps are all found in the `Comparison` folder. The process starts by running the program `Comparison.py`. The query's molecular fingerprint is evaluated, and the program uses Tanimoto similarity to go through the entire existing database from the previous step, finding the fingerprints that match the best. A shortlist for similarity is created as all the candidates with the highest Tanimoto similarities move on for further screening.

We utilize a vector database solution, Milvus, to store and manage these fingerprints efficiently. Milvus allows for high-speed retrieval and comparison of molecular fingerprints, significantly speeding up the search process. By partitioning the database, we ensure that searches are both fast and scalable, handling large datasets with ease.

This step ensures that only the most relevant candidates are considered for further analysis, streamlining the process and improving the accuracy of the results.

## Further Sorting

Further sorting is automatically run in the `Comparison.py` program. Sorting is further done by thermophysical properties and predicted toxicity. Most data used for this segment of comparison comes from OPERA, an open-source set of models that provide predictions on physicochemical properties. This is also used for toxicity evaluations. This segment utilizes multiple methods of comparison which are described as follows:

### Thermophysical Comparison

Results are computed by OPERA, and the search allows selection for turning on comparison of the following properties: Melting Point, Boiling Point, LogP, Vapor Pressure, and Henry's Law Constant. Candidates have their value compared against the value of the query, which is factored into the similarity score.

### Toxicity Evaluation

This also occurs during the OPERA evaluation. OPERA computes the following endpoints for toxicity: Log_BCF, CATMoS_EPA, CATMoS_LD50. These values are evaluated and factored into the total similarity score for each candidate. Higher toxicity decreases the similarity score.

### Structural Similarity

Another stage of structural similarity is performed past the fingerprinting method. Here, we also have OPERA compute several predicted structural properties of the candidates/query. These are mainly the Molecular Weight, Number of Rings, and Number of Lipinski Failures. These have little effect on the overall similarity score but are considered in computation. Like the thermophysical properties, these values for the candidates are compared against that of the query and factored into the similarity score.

### Synthetic Accessibility Scoring

The SA score is computed using a function provided by RDKit. The SA score predicts on a scale of 1-10 how difficult it would be to obtain the candidate. This can be used as a method of including cost-effectiveness in the similarity search. These values for each candidate are factored into the similarity score.

## Technical Details

### Vector Database Solution

The preprocessing step now uses a vector database solution for its functionality. This allows for efficient storage and retrieval of molecular fingerprints. The vector database is optimized for handling high-dimensional data, making it ideal for storing the 2048-bit long vectors used in this project. This ensures that the fingerprint comparison process is both fast and scalable.

### OPERA Integration

OPERA is an open-source software that uses machine learning techniques to predict thermophysical properties of compounds. It is used extensively in this project for both thermophysical comparison and toxicity evaluation. OPERA computes various properties such as Melting Point, Boiling Point, LogP, Vapor Pressure, and Henry's Law Constant. These properties are compared against those of the query molecule using a percent error format, where similarity is given on a scale of 1-10.

For toxicity evaluation, OPERA computes endpoints such as Log_BCF, CATMoS_EPA, and CATMoS_LD50. These values are factored into the total similarity score, with higher toxicity decreasing the similarity score.

### RDKit SA Scoring

RDKit provides a function to compute the Synthetic Accessibility (SA) score of a molecule. The SA score predicts on a scale of 1-10 how difficult it would be to synthesize the candidate molecule. This score is used to factor in the cost-effectiveness of the candidate in the similarity search. The SA score is computed for each candidate and factored into the overall similarity score.

## Integrating Your Own Models

CRS allows users to integrate their own machine learning models to enhance the chemical comparison process. By specifying the container names for your models, CRS will run these models during the comparison process and factor their outputs into the similarity score. This feature provides flexibility and customization, enabling users to tailor the system to their specific needs.

### Steps to Integrate Your Model

1. **Create a Docker Image**: Ensure that your model is accessible via a Flask API. The Flask API should be set up to receive SMILES strings and return the computed results.

2. **Create Dockerfile and Requirements**: Create a `Dockerfile` and `requirements.txt` for your model. The `Dockerfile` should specify the base image, copy the necessary files, install dependencies, and expose the required port. An example on how to do so with further explanation can be found in the `src/Model` folder of this repository.

3. **Build Docker Image**: Build your Docker image by running the following command in your terminal:
   ```sh
   docker build -t insert-image-name .
   ```
4. **Update Docker Compose**: Update the `docker-compose.yml` file of the CRS application to include your new Docker image as a service. Here is an exmaple of how to add your service as another one of the services:

   ```sh
   model_service_name:
    container_name: image-name
    image: image-name
    environment:
      - HTTP_PROXY=${HTTP_PROXY:-}
      - HTTPS_PROXY=${HTTPS_PROXY:-}
   ```

5. **Run the Model**: Run the `docker compose up -d` command as usual and enter the interactive shell using `docker exec -it CRS /bin/bash`. Now, when running the CRS by calling `src/main.py`, add an extra argument, `-m` or `--model`, that is followed by the name(s) of models you would like to add. These models should be named the same as the image-name written in the step above.

For further questions, please contact panair@sandia.gov
