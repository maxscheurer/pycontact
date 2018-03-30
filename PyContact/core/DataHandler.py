import pickle
import os

class DataHandler:
    """Handles the import or export of a session."""

    @staticmethod
    def importSessionFromFile(fileName):
        """Imports a saved session from 'filename'."""
        importDict = pickle.load(open(fileName, "rb"))
        contacts = importDict["contacts"]
        arguments = importDict["analyzer"][0:-1]
        trajArgs = importDict["trajectory"]
        maps = importDict["maps"]
        contactResults = importDict["analyzer"][-1]
        del importDict
        return [contacts, arguments, trajArgs, maps, contactResults]

    @staticmethod
    def writeSessionToFile(fileName, analysis):
        """Saves the current Session (analysis) at 'filename'."""
        analyzerArgs = [analysis.psf, analysis.dcd, analysis.cutoff, analysis.hbondcutoff,
                        analysis.hbondcutangle, analysis.sel1text, analysis.sel2text,
                        analysis.contactResults]
        trajArgs = analysis.getTrajectoryData()
        exportDict = {"contacts": analysis.finalAccumulatedContacts, "analyzer": analyzerArgs, "trajectory": trajArgs,
                      "maps": [analysis.lastMap1, analysis.lastMap2]}
        pickle.dump(exportDict, open(fileName, "wb"))


    @staticmethod
    def writeContactsToFile(filename, contacts):
        path, file_extension = os.path.splitext(filename[0])

        requestedParameters = ["contactTypeAsShortcut", "mean_score", "median_score", "hbond_percentage"]
        checkboxdict = {"mean_score": "Mean Score",
                             "hbond_percentage": "HBond Percentage", "median_score": "Median Score",
                             "contactTypeAsShortcut": "Contact Type",
                             "getScoreArray": "Score List", "hbondFramesScan": "Hydrogen Bond Frames"}

        tableHeadings = []
        for par in requestedParameters:
            tableHeadings.append(checkboxdict[par])

        f = open(path + ".txt", "w")
        row_format = " {:>20} " * (len(requestedParameters) + 1)
        f.write(row_format.format("", *tableHeadings))
        f.write("\n")
        for c in contacts:
            currentContactProperties = []
            for p in requestedParameters:
                # Python2:
                # exec('propertyToAdd = c.' + p + '()')
                code = compile('propertyToAdd = c.' + p + '()', '<string>', 'exec')
                ns = {}
                ns['c'] = c
                exec(code, ns)
                propertyToAdd = ns['propertyToAdd']
                if isinstance(propertyToAdd, float):
                    propertyToAdd = "{0:.3f}".format(propertyToAdd)
                currentContactProperties.append(propertyToAdd)
            f.write(row_format.format(c.human_readable_title(), *currentContactProperties))
            f.write("\n")
        f.close()
