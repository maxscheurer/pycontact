from PyContact.core.Scripting import PyContactJob, JobConfig
job = PyContactJob("PyContact/exampleData/rpn11_ubq.psf", "PyContact/exampleData/rpn11_ubq.dcd", "test", JobConfig(5.0, 2.5, 120, [0,0,1,1,0], [0,0,1,1,0], "segid RN11", "segid UBQ"))
job.runJob(1)
job.writeContactDataToFile()
