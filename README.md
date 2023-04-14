# pubSight: creating a summary figure of PUBMED publications direction #

**pubSight** , a visualization tools.

---
**Citation:**


Mahdi Baghbanzadeh, Daniel Kerchner, Ranojoy Chatterjee, Bahar Sayoldin, Ali Rahnavard (2023+).
**pubSight: creating a summary figure of publications in PUBMED database for literature direction**,
https://github.com/omicsEye/pubSight/.

---
# pubSight user manual

## Contents ##
* [pubSight](#pubsight)
* [Installation](#installation)
* [Getting Started with pubSight](#getting-started-with-pubsight)
* [Examples](#examples)
* [Support](#support)
------------------------------------------------------------------------------------------------------------------------------
## pubSight ##
pubSight is a Python package that visualizes the number of publications in specific fields using PubMed records.
It allows users to input a list of queries and make API calls to PubMed with Entrez to retrieve data.
The retrieved data is then processed and visualized using Matplotlib.


## Installation ##

To install pubSight, use pip:
```commandline
pip install pubSight
```
or use github:
```commandline
python -m pip install git+https://github.com/omicsEye/pubSight
```
* First install *conda*  
Go to the [Anaconda website](https://www.anaconda.com/) and download the latest version for your operating system.  
* For Windows users: do not forget to add `conda` to your system `path`
* Second is to check for conda availability  
open a terminal (or command line for Windows users) and run:
```
conda --version
```
it should out put something like:
```
conda 4.9.2
```
if not, you must make *conda* available to your system for further steps.
if you have problems adding conda to PATH, you can find instructions
[here](https://docs.anaconda.com/anaconda/user-guide/faq/).  

### Windows Linux Mac ###
If you are using an **Apple M1/M2 MAC** please go to the [Apple M1/M2 MAC](#apple-m1m2-mac) for installation
instructions.  
If you have a working conda on your system, you can safely skip to step three.  
If you are using windows, please make sure you have both git and Microsoft Visual C++ 14.0 or greater installed.
install [git](https://gitforwindows.org/)
[Microsoft C++ build tools](https://visualstudio.microsoft.com/visual-cpp-build-tools/)
In case you face issues with this step, [this link](https://github.com/pycaret/pycaret/issues/1254) may help you.
1) Create a new conda environment (let's call it pubSight_env) with the following command:
```
conda create --name pubSight_env python=3.9
```
2) Activate your conda environment:
```commandline
conda activate pubSight_env 
```
3) Install *pubSight*:
install with pip:
```commandline
pip install pubSight
```
or you can directly install if from GitHub:
```commandline
python -m pip install git+https://github.com/omicsEye/pubSight
```
### Apple M1/M2 MAC ###
1) Update/install Xcode Command Line Tools
  ```commandline
  xcode-select --install
  ```
2) Install [Brew](https://brew.sh/index_fr)
  ```commandline
  /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
  ```
3) Install libraries for brew
  ```commandline
  brew install cmake libomp
  ```
4) Install miniforge
  ```commandline
  brew install miniforge
  ```
5) Close the current terminal and open a new terminal
6) Create a new conda environment (let's call it pubSight_env) with the following command:
  ```commandline
  conda create --name pubSight_env python=3.9
  ```
7) Activate the conda environment
  ```commandline
  conda activate pubSight_env
  ```
8) Install packages from Conda
  ```commandline
  conda install lightgbm
  pip install xgboost
  ```
9) Finally, install *pubSight*:
install with pip:
```commandline
pip install pubSight
```
or you can directly install if from GitHub:
```commandline
python -m pip install git+https://github.com/omicsEye/pubSight
```
-----------------------------------------------------------------------------------------------------------------------

## Getting Started with pubSight ##
You can use pubSight in a Python IDE or through the command line.
The following are the available command line arguments:
```commandline
optional arguments:
  -h, --help            show this help message and exit
  --input INPUT, -i INPUT
                        path to file that has all the query terms
  --email EMAIL, -e EMAIL
                        a valid email address
  --out_dir OUT_DIR, -o OUT_DIR
                        path to report directory
  --color_palette COLOR_PALETTE, -c COLOR_PALETTE
                        Matplotlib color palette
  --group_legend, -g   make a group legend for the whole plot
  --pubmed_data PUBMED_DATA, -p PUBMED_DATA
                        path to a CSV dataframe that has four columns ['year', 'count', 'main_term', 'sub_term']
```
Here's an example usage:
```commandline
pubSight --input queries.txt --email example@example.com --out_dir report --color_palette viridis --group_legend
```
This command retrieves data from PubMed based on the queries listed in the file queries.txt. It uses the email 
address example@example.com for `Entrez`, saves the report to the directory report,
uses the color palette `viridis`, and creates a group legend for the whole plot.

## Examples ##
Here are some examples of visualizations created using pubSight:
![alt text](https://github.com/omicsEye/pubSight/blob/main/img/Figure_2.png)
![alt text](https://github.com/omicsEye/pubSight/blob/main/img/SFig1.png)
![alt text](https://github.com/omicsEye/pubSight/blob/main/img/supp1.png)


### Support ###

* Please submit your questions or issues with the software at [Issues tracker](https://github.com/omicsEye/pubSight/issues).
