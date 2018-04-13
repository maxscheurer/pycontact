from PyContact.core.Scripting import PyContactJob, JobConfig

# define input files and parameters
job = PyContactJob("/path/to/topology", "/path/to/trajectory", "title", JobConfig(5.0, 2.5, 120, [0,0,1,1,0], [0,0,1,1,0], "segid A", "segid B"))
# running the job on 4 cores
job.runJob(2)
# writing the session to file
job.writeSessionToFile()
