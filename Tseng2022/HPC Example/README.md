# Running a job using R script on an HPC cluster
The following Unix command lines (Terminal in Mac) provides a walk-through of how to submit an R job on Kamiak, Washington State University's high performance computing cluster. You can adapt the script to your HPC cluster.

## Description
To follow along, download the folder ~/Tseng2022/[HPC Example](https://github.com/viralemergence/PoxHost/tree/main/Tseng2022/HPC%20Example) from the [PoxHost repository](https://github.com/viralemergence/PoxHost).  In this example, we execute part three (“3. Boosted regression trees (BRT)”) of the markdown file [HostPrediction_Code.Rmd](https://github.com/viralemergence/PoxHost/blob/main/Tseng2022/Host%20Prediction%20Model/HostPrediction_Code.Rmd). This selection of code builds boosted regression tree models trained on infection and competence data to predict mammal hosts of Orthopoxviruses. The following files required for this job are listed below: 

| Files                     | Description                                                                          |
| ------------------------- |------------------------------------------------------------------------------------- |
| Kamiak_BRT_07Sep2022.R    | Model code in R script |
| HostData_clean.RData      | Input data obtained from part one "1. Data preparation" of [HostPrediction_Code.Rmd](https://github.com/viralemergence/PoxHost/blob/main/Tseng2022/Host%20Prediction%20Model/HostPrediction_Code.Rmd) |
| "Output" folder           | Empty folder location for saving model output                                        |
| KatieJob_07Sep2022.sh     | Bash shell script (aka Slurm script) for submitting job                              |

## Instructions 
1. Login to Kamiak to view your allocated compute partitions:
```R
$ ssh [username]@kamiak.wsu.edu
$ sinfo -p fernandez							#view partition
$ squeue -p fernandez							#view current jobs
```

2. Request a compute node (aka idev session)
```R
$ idev -p fernandez 							#receive access to one core for one hour (60 minutes default)
$ idev -p fernandez -t 180						#receive access to one core for three hours (note @cn##)
### $ exit									#returns you to login node (which is shared by all users)
```

3. Pre-install necessary R packages to your local node - EXAMPLE:
```R
### You only need to do this install step once for each package.
### Do not include install.packages() in your R script: attempts to install packages on Kamiak’s R will fail b/c of permissions. 
$ mkdir -pv~/lib/R_libs							#create local directory for R libraries
$ module load r/4.1.0							#load specific R version
$ R											#open R program
$ library(gbm)								#check if package already exists
$ install.packages(“gbm”, dependencies = TRUE, repos = "http://cran.rstudio.com")																		#installs package & dependencies and auto selects CRAN mirror
$ q()									#quit R program; to kill a command, Control+c
$ library(gbm)								#shows that gbm is now loaded
$ ls lib/R_libs/							#view local R libraries 
```

4. Prepare R script 	
```R
### In this example, we title the R script: Kamiak_BRT_07Sep2022
### Once packages are installed from the previous step, remove/comment out any install.packages() commands in your   R script
### Wherever you send your job to (.sh), it will set the file path of your job script as your home directory. Any file paths referenced IN your R script (e.g., reading in data or saving data) needs to match your home directory (file path of your job script). 
```

5. Create your batch job script (.sh) - sample code as follows
```R
#!/bin/bash							#specifies Unix shell to be used
#SBATCH --partition=fernandez					#specify partition to be used
#SBATCH --job-name=[job name]					#job name
#SBATCH --output=[job name]_%j.out				#standard output file
#SBATCH --error=[job name]_%j.err				#standard error file
#SBATCH --mail-type=ALL						#send email on job start, job end and job fault
#SBATCH --mail-user=[username]@wsu.edu				#address where job status emails will be used
#SBATCH --nodes=1						#node count
#SBATCH --ntasks=1						#total number of tasks across all nodes
#SBATCH --cpus-per-task=40					#cpu-cores per task (Fernandez node has 40 cores)
#SBATCH --mem-per-cpu=4G					#memory per cpu-core (4G is default) 
module load r/4.1.0						#load specified R version (newest is default)
Rscript --vanilla [R script title].R  				#’vanilla’ runs R script from clean environment
echo "Completed job $SLURM_JOBID on nodes $SLURM_JOB_NODELIST" 
```

6. Transfer files to Kamiak - EXAMPLE: R script (.R), data file (.RData), Output folder, shell script (.sh), etc.
```R
$ logout									#logout of Kamiak
$ scp -r /~/HPC Example/Kamiak_BRT_07Sep2022.sh [username]@kamiak.wsu.edu:~	#login password will be requested
$ scp -r /~/HPC Example/HostData_clean.RData [username]@kamiak.wsu.edu:~	#login password will be requested
$ scp -r /~/HPC Example/Output/ [username]@kamiak.wsu.edu:~			#login password will be requested
$ scp -r /~/HPC Example/KatieJob_07Sep2022.sh [username]@kamiak.wsu.edu:~	#login password will be requested
$ ssh [username]@kamiak.wsu.edu     #log back into Kamiak
$ ls										#check that contents were saved to home directory
```

7. Before submitting job, let’s test R script 
```R
$ head Kamiak_BRT_07Sep2022.R 			#view your R script
$ vim Kamiak_BRT_07Sep2022.R 			#with vim, you can now edit the file
$ R						#open R
$ load(“HostData_clean.RData”)			#let’s try running a couple lines of code
$ data <- poxdata				#after testing, to exit vim w/o saving: press Esc key, type :q, and hit Enter key
```

8. Submit job script to job queue/scheduler
```R
$ sbatch KatieJob_07Sep2022.sh			#you will get an email from SLURM notifying you the job is running and a 2nd email when it’s finished
$ squeue -u [username]				 
$ squeue -j [job number]									
$ exit
```

9. To check job progress, read contents of the output file; to check for errors, read contents of the error file
```R
$ cat [job name_number].out			#read the contents of the file
$ cat [job name_number].err			#e.g., [job name]_43310144.err
$ tail -4 [job name_number].out		#read the last four lines of the file (default is 10)
```

10. Once the job is complete:
```R
$ ssh [username]@kamiak.wsu.edu	    			#log back into Kamiak
$ ls							#check that data files were saved
$ for f in *\ *; do mv "$f" "${f// /_}"; done     #rename files replacing whitespaces with underscore
$ logout						#log out to proceed with file transfer
$ cd Downloads/Results/Kamiak				#cd [filepath] of where you want to save file
$ scp -r [username]@kamiak.wsu.edu:~/pcr_brts.RData .	#Copy from Kamiak – DO NOT FORGET “ .” at the end
```

11. To edit a file, use vim or nano (text editors for Unix/Linux OS) 
```R
$ vim [job name]_15Jun2022.sh				#view and edit a file
$ i							#insert text in file
$ q:!							#quit without saving changes
$ q:							#quit with saving changes
```

12. FAQ
How do I cancel a job?
```R
$ scancel [job number]    #to kill a job
```
How do I find more information about a specific partition and its allocated node(s)?
```R
sinfo -p fernandez       #view summary information about the fernandez partition (e.g., nodelist/id)
scontrol show node cn128  #view node information (e.g., total CPUs and memory, allocated CPUs, etc.)  
```
```R
$ R_LIBS_USER=~/lib/R_lbs
``` 

13. Additional resources:
- https://www.hpc-carpentry.org/hpc-shell/00-hpc-intro/index.html
- https://researchcomputing.princeton.edu/support/knowledge-base/slurm#arrays
- https://hpc.wsu.edu/users-guide/
- https://s3.wp.wsu.edu/uploads/sites/1122/2021/10/kamia¬¬¬k_cheat_sheet_vm.pdf	
