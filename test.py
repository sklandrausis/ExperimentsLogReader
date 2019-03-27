from ExperimentsLogReader.experimentsLogReader import LogReaderFactory, LogTypes


logFile = "exampleLogfiles/fs_example.log"
outputPath = "exampleOutputfiles/fs_example"
logs = LogReaderFactory.getLgReader(LogTypes.DBBC,logFile, outputPath, [], True)
logs.printLogs()


logFile = "exampleLogfiles/SDRtest.log"
outputPath = "exampleOutputfiles/SDRtest"
logs = LogReaderFactory.getLgReader(LogTypes.SDR,logFile, outputPath)
logs.printLogs()