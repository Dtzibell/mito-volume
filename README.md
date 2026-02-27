# Installation instructions
Download the repo:
```git clone https://www.github.com/dtzibell/mito-volume```

Open the directory and edit the ```config.ini``` file:
- ```DefaultDirectory``` is the directory where raw datasets are typically being kept, although they do not strictly have to be kept there.
- ```OutputDirectory``` is the directory where the output file ```output.xlsx``` will be created after running the script

Before running the script, ensure you have ```uv``` installed: ```pip install uv```\
Change into the directory ```mito-volume``` and run the script: ```uv run main.py```
