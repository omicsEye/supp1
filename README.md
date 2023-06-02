# pubSight: creating a summary figure of PUBMED publications direction #

**pubSight** , a visualization tools.

---
**Citation:**


Mahdi Baghbanzadeh, Bahar Sayoldin, Ali Rahnavard (2023+).
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

install [git](https://gitforwindows.org/)

1) Create a new conda environment (let's call it pubSight_env) with the following command:
```commandline
conda create --name pubSight_env python=3.9
```
2) Activate your conda environment:
```commandline
conda activate pubSight_env 
```
3) Install *pubSight*:
install with pip (NOT AVAILABLE NOW, please use the next option):
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
