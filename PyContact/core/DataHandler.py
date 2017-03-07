import pickle


class DataHandler():
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
