# this script takes the repo name as the first cmd line argument and creates the files:
# 'setup.sh' in the code directory of the project
# and 'README.md' in the notes directory of the project

# get the repo name from the first cmd line argument
repo_name=$1

# check if repo_name is equal to "template_speli_project"
if [ "$repo_name" = "template_speli_project" ]; then
  echo "repo_name is equal to template_speli_project, exiting script"
  exit 0
fi

# create the setup.sh file in the code directory using a here document to create a multi-line string and write it to a file.
cat <<EOF > code/setup.sh
#!/bin/bash
# This script is the basic setup script for the project
# run from the code directory with the following cmd:
# $ . ./setup.sh
# it will create the environment variables PROJECT_LOC & DATA_OUT in the current shell
cd ..
PROJECT_LOC=\$(pwd)
export PROJECT_LOC
cd code
mkdir -p /mnt/speliotes-lab/Projects/UK_ATLAS/IndivProj/\${USER}
DATA_OUT=/mnt/speliotes-lab/Projects/UK_ATLAS/IndivProj/\${USER}/${repo_name}/
export DATA_OUT
mkdir -p /mnt/speliotes-lab/Projects/UK_ATLAS/IndivProj/\${USER}/${repo_name}/
EOF

# create the markdown notes file in the notes directory using a here document to create a multi-line string and write it to a file.
cat <<EOF > notes/README.md
# Notes for ${repo_name}
---
## $(date +%Y%m%d)
- started project from template 🤖
EOF

echo "created setup.sh and README.md files in the code and notes directories respectively"
