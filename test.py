from ExperimentsLogReader.experimentsLogReader import LogReaderFactory, LogTypes

'''
logFile = "exampleLogfiles/fs_example.log"
outputPath = "exampleOutputfiles/fs_example"
logs = LogReaderFactory.getLogReader(LogTypes.DBBC, logFile, outputPath, ['225617.90', '620149.7', '2000.0'], "cepa")
logs.printLogs()
'''

logFile = "exampleLogfiles/SDRtest.log"
outputPath = "exampleOutputfiles/SDRtest"
logs = LogReaderFactory.getLogReader(LogTypes.SDR, logFile, outputPath)
logs.printLogs()