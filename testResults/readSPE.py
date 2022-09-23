import matplotlib
#matplotlib.use("Agg")

import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from math import log
import argparse, sys
import re

import matplotlib.colors as mcolors
colors=['blue', 'orange', 'red', mcolors.CSS4_COLORS['limegreen'], mcolors.CSS4_COLORS['darkgreen'], 'purple', 'red', 'gray']


from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('font',**{'family':'serif','serif':['Times'],'size':5})
rc('text', usetex=True)

def timesF(arr1, arr2):
    ret=[]
    for a in arr1:
        for b in arr2:
            ret.append((a,b))
    return ret

class Run:
    def __init__(self,name, test_name, discr, time, length, err, refinements, threads, functionType, jump, guessInitialAngles, initTime=-1.0, endTime=-1.0, power_file=""):
        self.name=name
        self.devName=""
        self.test_name=test_name
        self.discr=discr
        self.time=time
        self.length=length
        self.err=err
        self.refinements=refinements
        self.threads=threads
        self.functionType=int(functionType)
        self.jump=int(jump)
        self.power_file=power_file
        self.power_cons=[]
        self.n=1
        if guessInitialAngles=="true":
            self.guessInitialAngles=True
        else:
            self.guessInitialAngles=False
        self.initTime=initTime
        self.endTime=endTime

    def __eq__(self, other):
        return (self.devName==other.devName and
                self.test_name==other.test_name and 
                self.discr==other.discr and 
                self.refinements==other.refinements and 
                self.threads==other.threads and 
                self.guessInitialAngles==other.guessInitialAngles and 
                self.functionType==other.functionType and 
                self.jump==other.jump)

    def __str__(self):
        return ("dev_name: "+str(self.name)+
                " test_name: "+str(self.test_name)+
                " discr: "+str(self.discr)+
                " refinements: "+str(self.refinements)+
                " threads: "+str(self.threads)+
                " guessInitialAngles: "+str(self.guessInitialAngles)+
                " functionType: "+str(self.functionType)+
                " jump: "+str(self.jump))
        #return (self.name+" "+self.test_name+" "+str(self.discr)+" "+str(self.time)+" "+self.power_file)


def main():
    speRuns=[]
    with open("spe.json") as f:
        data=json.load(f)['SPE']
        for p in data:
            r=Run(p['name'], p['test_name'], p['discr'], p['time'], p['length'], p['err'], p['refinements'], p['threads'], p['functionType'], p['jump'], p['guessInitialAngles'], p['initTime'], p['endTime'])
            if r.refinements==1:
                continue
            found=False
            for i in range(len(speRuns)):
                if speRuns[i]==r:
                    speRuns[i].time+=r.time 
                    speRuns[i].n+=1
                    found=True
            if not found:
                speRuns.append(r)

    testNames=[]
    discrs=[]
    refs=[]
    deviceNames=[]
    for r in speRuns:
        if r.name not in deviceNames:
            deviceNames.append(r.name)
        if r.test_name not in testNames:
            testNames.append(r.test_name)
        if r.discr not in discrs:
            discrs.append(r.discr)
        if r.refinements not in refs:
            refs.append(r.refinements)
        if r.n>1:
            r.time/=r.n

    for tn in testNames:
        for discr in discrs:
            for ref in refs:
                for name in deviceNames:
                    for r in speRuns:
                        if r.name==name and r.discr==discr and r.refinements==ref and r.test_name==tn:
                            print({"name": name, "ref": ref, "discr": discr, "tn": tn, "time": r.time})
                print()

    nGraphs=0
    for tn in range(len(testNames)):
        break
        print(testNames[tn])

        width=0.1
        fig, ax=plt.subplots()

        for n in range(len(deviceNames)):
            values=[]
            labels=timesF(discrs, refs)
            x=np.arange(len(labels))
            for l in range(len(labels)):
                (dis,ref)=labels[l]
                #print("looking for: "+names[n]+" with discr: "+str(dis)+", ref: "+str(ref))
                for r in speRuns:
                    #print(r)
                    if r.name==deviceNames[n]:
                        #print("r.devName")
                        #print("r.test_name <"+r.test_name+"> <"+testNames[tn]+"> "+str(len(testNames)))
                        if r.test_name==testNames[tn]: 
                            if r.discr==dis: 
                                #print("r.discr")
                                if r.refinements==ref:  
                                    values.append(log(r.time))
      
            print(values)
            print(labels)
            print(x)
            a=1
            if len(deviceNames)%2!=1:
                a=0.5
            if n%2==0:
                ax.bar(x-width*(a+int(n/2)), values, width, label=deviceNames[n], color=colors[n])
            else:
                ax.bar(x+width*(a+int(n/2)), values, width, label=deviceNames[n], color=colors[n])
        
        ax.set_title(testNames[tn])
        ax.set_xlabel("Disc $k$, ref $m$")
        ax.set_xticklabels([0]+[str((d,r)) for (d,r) in labels ])
        ax.set_ylabel("Logarithmic times")
        ax.legend()
        plt.savefig("SPEtimesGraphsAllInOneLog"+str(nGraphs)+".pdf", transparent=True)
        nGraphs+=1

    print(deviceNames)
    deviceNames.sort(reverse=True)
    print(deviceNames)
    print("\\subsection{Simulation Data}")
    print("\\begin{center}")
    print("\\begin{tabular}{c|c|c", end="")
    for n in deviceNames:
        print("|c", end="")
    print("}")
    print("\t&&&\\multicolumn{"+str(len(deviceNames))+"}{c}{Times(ms)}\\\\")
    print("\t\\#Points&Disc&Ref", end="")
    for n in deviceNames:
        print("&"+n, end="")
    print("\\\\\n\\hline")
    for tnn in range(len(testNames)):
        tn=testNames[tnn]
        print("\t\\multirow{"+str(len(discrs)*len(refs))+"}{*}{"+tn+"}", end="")
        for d in range(len(discrs)):
            discr=discrs[d]
            if d!=0:
                print("\t\t", end="")
            print("&\\multirow{"+str(len(refs))+"}{*}{"+str(discr)+"}", end="")
            for r in range(len(refs)):
                ref=refs[r]
                if r==0:
                    print("&"+str(ref)+"&", end="")
                else:
                    print("\t\t\t&&"+str(ref)+"&", end="")
                for dv in range(len(deviceNames)): 
                    for r in speRuns:
                        if r.name==deviceNames[dv] and r.test_name==tn and r.discr==discr and r.refinements==ref:
                            if dv!=len(deviceNames)-1:
                                print("{:.2f}".format(r.time/1000.0), end="&")
                            else:
                                print("{:.2f}".format(r.time/1000.0), end="\\\\\n")

        print("\n\\hline")
    print("\\end{tabular}")
    print("\\end{center}\n\n\n\n\n\n\n")





if __name__ == '__main__':
    main()