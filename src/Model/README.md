_CRS Model Integration through Images_

This code serves as an example for the setup of a Docker Image in order to be integrated with the CRS. This code inputs SMILES and returns the molecular weight. Notice that the added files which a model requires to be compatible with the CRS are in this case, the app.py and the dockerfile/requirements.txt. App.py opens a port for the CRS to access your model during a run. The requirements.txt should be updated with all python packages required by the program to be installed by pip. The other file in the repository, function.py, is an example of the code where the added models are calculated.

Once all files are created and organized, cd into the directory of the program and run use the terminal to run 'docker build -t {image-name} .', where the image-name is whatever you would like to call your image. Now, your docker image exists on your machine for your docker daemon to access.

Return to the CRS application and update the docker-compose file that you run to add in this image as a service for the CRS application. Then, by using the CLI option of the CRS, these models can be added on with the -c option and specifying the name of your image.
