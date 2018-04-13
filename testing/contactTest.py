from PyContact.core.Scripting import PyContactJob, JobConfig

job = PyContactJob("../PyContact/exampleData/rpn11_ubq.psf",
                   "../PyContact/exampleData/rpn11_ubq.dcd",
                   "title",
                   JobConfig(5.0, 2.5, 120, [0,0,1,1,0], [0,0,1,1,0], "segid RN11", "segid UBQ"))
job.runJob(2)

for contact in job.analyzer.finalAccumulatedContacts:
    print(contact.scoreArray[0])
