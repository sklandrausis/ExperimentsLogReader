#! /usr/bin/python3
# -*- coding: utf-8 -*-

import sys
import numpy as np
import matplotlib
from astropy.stats import LombScargle
from gatspy.periodic import LombScargleFast
import astropy.units as u
import datetime
import json
import argparse
from operator import itemgetter
from astropy.time import Time
import time

from PyQt5.QtWidgets import (QWidget, QGridLayout, QApplication, QPushButton, QLabel, QLineEdit, QDesktopWidget, QComboBox)
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import Qt

from ploting_qt5 import  Plot
from result import  Result
from parsers._configparser import ConfigParser
from help import *
from monitor.months import Months
from monitor.monitoringViewHelper import MonitoringViewHelper
from jsonschema._validators import maximum

matplotlib.use('Qt5Agg')

def parseArguments():
    parser = argparse.ArgumentParser(description='''Monitoring velocity amplitudes in time. ''', epilog="""Monitor.""")
    parser.add_argument("-c", "--config", help="Configuration cfg file", type=str, default="config/config.cfg")
    parser.add_argument("-v","--version", action="version", version='%(prog)s - Version 2.0')
    args = parser.parse_args()

    return args

class Monitoring_View(QWidget):
        def __init__(self, iteration_list, location_list, source, output_path):
            super().__init__()
            self.grid = QGridLayout()
            self.setLayout(self.grid)
            self.grid.setSpacing(10)
            
            comboBox = QComboBox(self)
            comboBox.addItem("polarization AVG")
            comboBox.addItem("polarization U1")
            comboBox.addItem("polarization U9")
            comboBox.addItem("ALL")
            self.grid.addWidget(comboBox, 1, 2)
            comboBox.activated[str].connect(self.getPolarization)
            
            self.polarization = "polarization AVG"
            self.labels = list()
            self.lines = list()
            self.iteration_list = iteration_list
            self.location_list = location_list
            self.source = source
            self.output_path = output_path
            self.new_spectre = True
            self.spectrumSet = set()
            self.plotList = list()
            
            self.setWindowTitle("Monitoring")
            
            self.months = Months()
                        
        def _addWidget(self, widget, row, collon):
            self.grid.addWidget(widget, row, collon)
            
        def getGrid(self):
            return self.grid
        
        def setPolarization(self, polarization):
            self.polarization = polarization
            
        def setLabels(self, labels):
            self.labels = labels
             
        def setLines(self, lines):
            self.lines = lines
            
        def __setVisible(self, label, labels):
            index = labels.index(label)
            self.lines[index].set_visible(True)
            self.lines[index].set_picker(5)
            
        def __unSetVisibel(self, label, labels):
            index = labels.index(label)
            self.lines[index].set_visible(False)
            self.lines[index].set_picker(False)
        
        def getIndexiesOfPolarization(self, labels):
            if self.polarization == "ALL":
                all(i.set_visible(True) for i in self.lines)
                all(i.set_picker(5) for i in self.lines)
                
            else:
                for label in labels:
                    if self.polarization in label:
                        self.__setVisible(label, labels)
                    elif self.polarization not in label:
                        self.__unSetVisibel(label, labels)
                        
        def getPolarization(self, polarization):
            self.setPolarization(polarization)
            self.getIndexiesOfPolarization(self.labels)
            
        def keyPressEvent(self, e):
            if e.key() == Qt.Key_Shift:
                self.new_spectre = True
                
        def chooseSpectrum(self, event):
            thisline = event.artist
            xdata = thisline.get_xdata()
            ind = event.ind
            index = [ind][0]
            polarization = thisline.get_label().split()[1]
            spectraFileName = self.source + "_" + MonitoringViewHelper.formatDate(xdata, index) + "_" + MonitoringViewHelper.getLocation(self.location_list, int(index[0])) + "_"  + MonitoringViewHelper.getIteration(self.iteration_list, int(index[0])) + ".dat"
            self.plotSpecter(spectraFileName, polarization)
        
        def plotSpecter(self, spectraFileName, polarization):
            if polarization == "U1":
                amplitude_colon = 1
            elif polarization == "U9":
                amplitude_colon = 2
            elif polarization == "AVG":
                amplitude_colon = 3
                
            spectraFileName = self.output_path + spectraFileName
            
            data = np.fromfile(spectraFileName, dtype="float64", count=-1, sep=" ") .reshape((file_len(spectraFileName),4))
            tmpDate = spectraFileName.split("/")[-1].split("_")
            tmpDate[-4] = self.months.getMonthNumber([tmpDate[-4]][0])
            plot_name = datetime.datetime.strptime( " ".join(tmpDate[1:-2]), "%H %M %S %d %m %Y") 
            
            x = data[:, [0]]
            y = data[:, [amplitude_colon]]
            
            if self.new_spectre:
                self.Spectre_View = Spectre_View()
                self.spectrumSet.add(self.Spectre_View)
                self.spectrPlot = Plot()
                self.spectrPlot.creatPlot(self.Spectre_View.getGrid(), "Velocity (km sec$^{-1}$)", "Flux density (Jy)", spectraFileName.split("/")[-1].split("_")[0], (1,0))
                self.Spectre_View._addWidget(self.spectrPlot, 0, 0)
                self.Spectre_View.show()
                self.new_spectre = False
                self.plotList.clear()
            
            if plot_name not in self.plotList:
                self.plotList.append(plot_name)
                self.spectrPlot.plot(x,y, "-", label=plot_name)
                self.spectrPlot.canvasShow()
            
class Spectre_View(QWidget):
        def __init__(self):
            super().__init__()
            self.grid = QGridLayout()
            self.setLayout(self.grid)
            self.grid.setSpacing(10)
            
            self.setWindowTitle(" ")
            
        def _addWidget(self, widget, row, collon):
            self.grid.addWidget(widget, row, collon)
            
        def getGrid(self):
            return self.grid         
            
class Monitoring(QWidget):
    def __init__(self, configFilePath):
        super().__init__()
        self.setWindowIcon(QIcon('viraclogo.png'))
        self.center()
        self.configFilePath = configFilePath
        self.grid = QGridLayout()
        self.setLayout(self.grid)
        self.grid.setSpacing(10)
        self.chooseSource()
        self.new_spectre = True
        self.months = Months()
        
    def center(self):
        qr = self.frameGeometry()
        cp = QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())
        
    def chooseSource(self):
        self.setWindowTitle("Choose Source")
        self.chooseLabel = QLabel("Choose source")
        self.grid.addWidget(self.chooseLabel, 0, 0)
        self.sourceInput = QLineEdit()
        self.grid.addWidget(self.sourceInput, 1, 0)
        self.chooseSource = QPushButton("Ok", self)
        self.chooseSource.clicked.connect(self.plot)
        self.chooseSource.setStyleSheet("background-color: green")
        self.grid.addWidget(self.chooseSource, 1, 1)
        
    def keyPressEvent(self, e):
        if e.key() == Qt.Key_Return:
            if len (self.sourceInput.text()) > 1:
                self.plot()
                                                                                                                                            
    def plotMonitoring(self, resultDir, source_velocities, source):
        result_list = list()
        velocitie_dict = {"u1":dict(), "u9":dict(), "avg":dict()}
        iteration_list = list()
        date_list = list()
        
        for velocitie in velocitie_dict:
            for vel in source_velocities:
                velocitie_dict[velocitie][vel] = list()
                                
        resultFileName = source + ".json"
        
        with open(resultDir + resultFileName) as result_data:    
                results = json.load(result_data)
        
        self.location_list = list()
        labels2 = list()
        for experiment in results:
                scanData = results[experiment]
                date = scanData["Date"]
                location = scanData["location"]
                amplitudes_for_u1 = scanData["polarizationU1"] # Got poitns for all experiments for polarization u1
                amplitudes_for_u9 = scanData["polarizationU9"] # Got poitns for all experiments for polarization u9
                amplitudes_for_uAVG = scanData["polarizationAVG"] # Got poitns for all experiments for polarization uAVG
                iter_number = scanData["Iteration_number"]
                specie = scanData["specie"]
                dates = date.split("_")
                monthsNumber = dates[1]
                dates[1] = self.months.getMonthNumber([monthsNumber][0])
                date = scanData["time"].replace(":", " ") + " " +  " ".join(dates) 
            
                result = Result(location, datetime.datetime.strptime(date, '%H %M %S %d %m %Y'), amplitudes_for_u1, amplitudes_for_u9, amplitudes_for_uAVG, iter_number, specie)
                
                result_list.append(dict(result))
              
        result_list = sorted(result_list, key=itemgetter('date'), reverse=False)
         
        for experiment in result_list:
            u1 = experiment["polarizationU1"]
            u9 = experiment["polarizationU9"]
            avg = experiment["polarizationUAVG"]
            iteration_list.append(experiment["iteration_number"])
            date_list.append(experiment["date"])
            location = experiment["location"]
            specie = experiment["specie"]
            
            self.location_list.append(location)
            
            for i in u1:
                for vel in source_velocities:
                    if float(vel) == float(i[0]):
                        velocitie_dict["u1"][vel].append(i[1]) 
            
            for j in u9:
                for vel in source_velocities:
                    if float(vel) == float(j[0]):
                        velocitie_dict["u9"][vel].append(j[1]) 
                        
            for k in avg:
                for vel in source_velocities:
                    if float(vel) == float(k[0]):
                        velocitie_dict["avg"][vel].append(k[1]) 
                         
            label = "Station is " + location + "\n" + "Date is " + experiment["date"].strftime('%d %m %Y') + "\n " + "Iteration number " + str(experiment["iteration_number"]) + "\n " + "Specie " + str(specie)
            labels2.append(label)
            
        self.iteration_list = iteration_list
       
        x = list()
        for a in range(0, len(date_list)):
            x.append(date_list[a])
        
        Symbols =  ["*", "o", "v", "^", "<", ">", "1", "2", "3", "4"]
        
        lines = list()
        
        self.Monitoring_View = Monitoring_View(self.iteration_list, self.location_list, self.source, self.output_path)
        self.monitoringPlot = Plot()
        self.monitoringPlot.creatPlot(self.Monitoring_View.getGrid(), "Time", "Flux density (Jy)", None, (1,0))
        
        def convertDatetimeObjectToMJD(time):
            time=time.isoformat()
            t=Time(time, format='isot')
            return t.mjd 
        
        for i in range(0, len(source_velocities)):
            l1, = self.monitoringPlot.plot(x, velocitie_dict["u1"][source_velocities[i]], Symbols[i]+"r", label="polarization U1 " + "Velocity " + source_velocities[i], visible=False, picker=False)
            l2, = self.monitoringPlot.plot(x, velocitie_dict["u9"][source_velocities[i]], Symbols[i]+"g", label="polarization U9 " + "Velocity " + source_velocities[i], visible=False, picker=False)
            l3, = self.monitoringPlot.plot(x, velocitie_dict["avg"][source_velocities[i]], Symbols[i]+"b", label="polarization AVG " + "Velocity " + source_velocities[i], visible=True, picker=5)
            
            lines.append(l1)
            lines.append(l2)
            lines.append(l3)
            
        t = np.array([convertDatetimeObjectToMJD(i) for i in x], dtype="float64")
        
        self.Monitoring_View._addWidget(self.monitoringPlot, 0, 0)
        self.monitoringPlot.setXtics(x, [date.strftime("%H %M %d %m %Y") for date in  date_list], '30')
        
        self.monitoringPlot.addCursor(labels2)
        labels = [str(line.get_label()) for line in lines]
        
        self.Monitoring_View.setLabels(labels)
        self.Monitoring_View.setLines(lines)
        self.monitoringPlot.addPickEvent(self.Monitoring_View.chooseSpectrum)
        
        y =  velocitie_dict["avg"][source_velocities[0]]

        error = np.array(y) * 0.1
        ls  = LombScargle(t,  y, error, fit_mean=True)
        
        def dateDelta(d1, d2):
            return abs(d1 -d2)
        
        def getMaxDateDelta():
            
            maxDelta = list()
            for x in range(0, len(t) -1):
                maxDelta.append(dateDelta(t[x], t[x + 1]))

            return np.max(maxDelta)

        nyquist_factor = 2 * getMaxDateDelta()
        print ("nyquist_factor", nyquist_factor)
        frequency, power = ls.autopower(method='fastchi2', normalization='model', nyquist_factor=nyquist_factor, minimum_frequency=0.002, maximum_frequency=dateDelta(t[0], t[-1]), samples_per_peak=20)
        #frequency, power = ls.autopower()
        
        period_days = 1. / frequency
        best_period = period_days[np.argmax(power)]
        print("Best period: {0:.2f} hours".format(24 * best_period))
                 
        self.periodPlot = Plot()
        self.periodPlot.creatPlot(self.Monitoring_View.getGrid(), "Period (days)", "Power", None, (1,1))
        self.periodPlot.plot(period_days, power, "r*", label="polarization AVG " + "Velocity " + source_velocities[0], rasterized=True)
        self.Monitoring_View._addWidget(self.periodPlot, 0, 1)
        
        self.Monitoring_View.showMaximized()
        self.Monitoring_View.show()
        
    def plot(self):
        config = ConfigParser.getInstance()
        config.CreateConfig(self.configFilePath)
        resultDir = config.getConfig('paths', "resultFilePath")
        self.output_path = config.getConfig('paths', "outputFilePath")
        self.source = self.sourceInput.text()
        source_velocities = config.getConfig('velocities', self.source).split(",")
            
        self.plotMonitoring(resultDir, source_velocities, self.source)
        
def main():
    # Parse the arguments
    args = parseArguments()
    configFilePath = str(args.__dict__["config"])
     
    #Create App
    qApp = QApplication(sys.argv)
    aw = Monitoring(configFilePath)
    aw.show()
    sys.exit(qApp.exec_())
    
    sys.exit(0)
        
if __name__=="__main__":
    main()
    