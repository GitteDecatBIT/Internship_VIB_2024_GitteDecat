# Internship_VIB_2024_GitteDecat
Have to install poetry ro be able to install biocypher 
pipx install poetry


poetry new [name_of_project]
cd [name_of_project]
poetry add biocypher

Poetry creates a virtual environment for you

install pypath-omnipath--> cause pypath was not recognised 

poetry add pypath-omnipath

same for bioregistry  
poetry add bioregistry 


Do it with the Crossbar github link

git clone link
cd folder that was created 
peotry install
poetry shell 
peotry run python file (knowledge graph)

when running it --> it keeps going for a long time at the update_process part
^CTraceback (most recent call last):
  File "/home/guest/.cache/pypoetry/virtualenvs/bccb-FmSLJcNW-py3.11/lib/python3.11/site-packages/pypath/share/curl.py", line 1548, in update_progress
    def update_progress(self, download_total, downloaded, upload_total,
