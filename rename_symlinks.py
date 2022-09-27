import os, sys
import shutil


directory_initial = "/g/korbel/tsapalou/anaconda3/mosaicatcher-pipeline/POOLING_PSEUDOPOOL/pseudopool/all/"

directory = "/g/korbel2/weber/MosaiCatcher_files/POOLING/PSEUDOPOOL/all/"


[os.remove(directory + e) for e in os.listdir(directory)] 

[os.symlink(directory_initial + e, directory + e) for e in os.listdir(directory_initial) if "bam" in e and "NA12878" not in e]

print(os.listdir(directory))
[os.rename(directory + e, directory + "_".join(e.split("_")[0:1] + e.split("_")[2:]))  for e in os.listdir(directory) if "bam" in e]
[os.rename(directory + e, directory + e.replace('x02PE20', "_").replace("x01PE20", "_").replace("102PE20", "_"))  for e in os.listdir(directory) if "bam" in e]
print(os.listdir(directory))