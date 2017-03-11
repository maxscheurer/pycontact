import pickle


class DataHandler:
    """docstring"""

    @staticmethod
    def importSessionFromFile(fileName):
        """docstring"""
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
        analyzerArgs = [analysis.psf, analysis.dcd, analysis.cutoff, analysis.hbondcutoff,
                        analysis.hbondcutangle, analysis.sel1text, analysis.sel2text,
                        analysis.contactResults]
        trajArgs = analysis.getTrajectoryData()
        exportDict = {"contacts": analysis.finalAccumulatedContacts, "analyzer": analyzerArgs, "trajectory": trajArgs,
                      "maps": [analysis.lastMap1, analysis.lastMap2]}
        pickle.dump(exportDict, open(fileName, "wb"))
