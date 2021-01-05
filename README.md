# ExperimentsLogReader
Experiment log file reader
Parse experiment log file created by field system

# Table of Contents
- [Dependencies](#dependencies)
- [Usage](#usage)
- [Installation](#installation)
- [Changelog](#changelog)
- [Getting Help](#getting-help)
- [Acknowledgements](#acknowledgements)

## Dependencies
- python 3
  - tkinter
  - configparser

## Usage
```python
# For DBBC
logFile = "exampleLogfiles/fs_example.log"
outputPath = "exampleOutputfiles/fs_example"
logs = LogReaderFactory.getLogReader(LogTypes.DBBC, logFile, outputPath, ['225617.90', '620149.7', '2000.0'], "cepa")
logs.printLogs()

# For SDR
logFile = "exampleLogfiles/SDRtest.log"
outputPath = "exampleOutputfiles/SDRtest"
logs = LogReaderFactory.getLogReader(LogTypes.SDR, logFile, outputPath)
logs.printLogs()
```

## Installation
The quickest way to install is using PyPI:
```bash
 sudo pip3 install ExperimentsLogReader
 ```

## Getting Help

Bug reports, feature requests and make contributions (e.g. code patches) can be reported by opening a &quot;new issue&quot; ticket on GitHub. Please give as much information (e.g. the software component, version) as you can in the ticket. For bugs, it is extremely useful if a small self-contained code snippet that reproduces the problem is provided.

## Acknowledgements
This software was written by Jānis Šteinbergs under the supervision of Artis Aberfelds. If you make use of this software to get results that appear in a publication or presentation please include this acknowledgement: &quot;We have made use of ExperimentsLogReader, a tool developed by Jānis Šteinbergs.&quot;

This work is the result of project implementation: «Physical and chemical processes in the interstellar medium», No 1.1.1.1/16/A/213 supported by ERDF​.

This work was supported by Latvian Council of Science Project “Research of Galactic Masers” Nr.: lzp-2018/1-0291.




