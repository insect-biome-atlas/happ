#!/usr/bin/env python

import requests
import os
import gzip as gz
from io import BytesIO
import yaml
import sys

# Read test config
sys.stderr.write("Loading config from test/configfile.yml\n")
configfile = "test/configfile.yml"
with open(configfile, "r") as fhin:
    config = yaml.safe_load(fhin)

sys.stderr.write("Download sintax reference from Figshare\n")
# Download sintax ref from figshare
BASE_URL = "https://api.figshare.com/v2"
id = 38787078
os.makedirs("resources/sintax", exist_ok=True)
response = requests.get(f"{BASE_URL}/file/download/{id}")
with gz.open(BytesIO(response.content), "rt") as fhin, open(
    "resources/sintax/sintax.ref.fasta", "wt"
) as fhout:
    x = fhout.write(fhin.read())
sys.stderr.write("Sintax reference downloaded to resources/sintax/sintax.ref.fasta\n")

sys.stderr.write("Updating config\n")
# Set sintax ref in config
config["sintax"] = {"ref": "resources/sintax/sintax.ref.fasta"}

sys.stderr.write(
    "Downloading epa-ng reference files from GitHub (https://github.com/insect-biome-atlas/paper-bioinformatic-methods)\n"
)
# Download Chesters tree files
CHESTERS_BASE_URL = "https://raw.githubusercontent.com/insect-biome-atlas/paper-bioinformatic-methods/refs/heads/main/chesters/data/chesters_tree"
chesters_data = {
    "fasta": {
        "name": "chesters_new_outgroups_aligned.trim0.9.fasta",
        "url": f"{CHESTERS_BASE_URL}/chesters_new_outgroups_aligned.trim0.9.fasta",
    },
    "tree": {
        "name": "chesters_new_outgroups.nwk",
        "url": f"{CHESTERS_BASE_URL}/chesters_new_outgroups.nwk",
    },
    "taxonomy": {"name": "taxonomy.tsv", "url": f"{CHESTERS_BASE_URL}/taxonomy.tsv"},
}

os.makedirs("resources/chesters/", exist_ok=True)
for key, value in chesters_data.items():
    response = requests.get(value["url"])
    with open(f"resources/chesters/{value['name']}", "wt") as fhout:
        fhout.write(response.text)
sys.stderr.write("Reference files downloaded to resources/chesters\n")

sys.stderr.write("Updating config\n")
config["epa-ng"] = {
    "tree": "resources/chesters/chesters_new_outgroups.nwk",
    "ref_taxonomy": "resources/chesters/taxonomy.tsv",
    "msa": "resources/chesters/chesters_new_outgroups_aligned.trim0.9.fasta",
    "chunk_size": 100,
}

with open("test/configfile.yml", "w") as fhout:
    yaml.safe_dump(config, fhout)
