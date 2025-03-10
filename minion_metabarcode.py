#!/usr/bin/env /opt/miniconda3/envs/snakemake/bin/python

import snakemake
import os
import glob
import pandas as pd
from shutil import copyfile

##### load config and sample sheets #####
input_directory = os.getcwd()
source_con='/home/shared_docs/pipelines/minion_metabarcode/config.yaml'
des_con=os.path.abspath(input_directory) + "/config.yaml"
copyfile(source_con,des_con)


print("##########\nThis pipeline run MinION metabarcoding pipeline. Must be run from a directory with demultiplexed reads in a directory named 'data'.\n##########\n")


usr_db = input("Reference database to use? (phytoplasma or phytophthora): ")
usr_minlength = input("Minimum read length to use?: ")
usr_maxlength = input("Maximum read length to use?: ")
usr_minqual = input("Minimum read quality to use?: ")
usr_cores = input("How many threads to use (default 40)?: ")


snakemake.snakemake("/home/shared_docs/pipelines/minion_metabarcode/Snakefile", conda_prefix="/home/shared_docs/pipeline_envs", cores=int(usr_cores), use_conda=True, keepgoing=True, resources={"load":1}, config={"minlength":usr_minlength, "maxlength":usr_maxlength, "minqual":usr_minqual })

