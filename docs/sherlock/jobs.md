# Jobs on Sherlock

[Running jobs](https://www.sherlock.stanford.edu/docs/user-guide/running-jobs/) on Sherlock is a little different. In order to fairly distribute cluster resources, Sherlock uses a scheduler called SLURM. Users define what work they want done in job files and submit them to the scheduler, which allocates the requested compute resources when it can. When your job is allocated resources, compute nodes with the requested resources will run your job automatically.

## Batch Jobs
> A job is simply an instance of your program, for example your R, Python or Matlab script that is submitted to and executed by the scheduler (Slurm). When you submit a job with the `sbatch` command it's called a batch job and it will either run immediately or will pend (wait) in the queue.

The [Sherlock jobs docs](https://www.sherlock.stanford.edu/docs/user-guide/running-jobs/#batch-jobs) are fairly comprehensive and will provide more detail than we can here. We will provide a few example batch scripts to demonstrate some standard uses. 

## Interactive jobs
It can be helpful when initially creating a job to work on it interactively on a compute node like the ones that will run your job. You can request a compute node via the command `sh_dev`. By default, sh_dev allocates one core and 4 GB of memory on one node for one hour. See [the docs](https://www.sherlock.stanford.edu/docs/user-guide/running-jobs/#interactive-jobs) for more details about requesting an interactive node.

## Interactive applications
Sherlock provides several [interactive applications](https://www.sherlock.stanford.edu/docs/user-guide/ondemand/?h=jupyter#interactive-applications) if you need to run GUI based interactive software. When working with arnie you will most likely use [JupyterNotebooks](https://www.sherlock.stanford.edu/docs/user-guide/ondemand/?h=jupyter#jupyter-notebooks) or [JupyterLab](https://www.sherlock.stanford.edu/docs/user-guide/ondemand/?h=jupyter#jupyterlab) to interactively explore your research questions. 

## Important Notes
We have run into some common issues using Sherlock over the years. Here's a non-comprehensive list of things to watch out for while using Sherlock.

- **Permissions errors**:

  If you plan on working on projects with group members and share files in `$GROUPHOME`, remember to set the permissions of files you create to be group accessible. By default, files will be read-only for group members. `chmod -R 770 /path/to/file` will allow group members to read, write, and execute shared files in `$GROUPHOME`.

- **Large array jobs impacting $GROUPHOME**:

  Be careful about accessing resources in `$GROUPHOME` when running large array jobs. Thousands of the same job accessing the same files on `$GROUPHOME` can slow down file access for other lab members. In many cases, your code may not be the one accessing files in `$GROUPHOME`, but a predictor you're using might (`spotrna` causes this issue often). The best solution is to copy the files you're accessing to your `$SCRATCH` folder and access them there. If you have large array jobs requesting thousands of nodes, you may want to copy the files to the node's `$LSCRATCH` instead. See the array job example sbatch file for more details. 


