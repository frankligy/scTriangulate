## Instructions

So in each python script I listed here, you can see I have a shebang line at the very top, indicating which python environment I used when executing the code. When reproducing the results or inspecting the intermediate files I have on Synapse, it may be necessary to make sure the dependencies are the same. For example, if you'd like to utilize my pickle file and the dependencies are different, you may not be able to read in my pickle file.

Due to the long time range of this project, I alternatively used two environments, one I called `old_env`, which based on python 3.6, another I called `new_env`, which based on python 3.7, both are on Linux system. I have two yml files in this folder for you to exactly reproduce my environments.

```
# using old
conda env -f sctri_old_env_py36_linux.yml -p ./py36_old_linux_env

# using new
conda env -f sctri_new_env_py37_linux.yml -p ./py37_new_linux_env
```