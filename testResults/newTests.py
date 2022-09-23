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

from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('font',**{'family':'serif','serif':['Times'],'size':14})
rc('text', usetex=True)
LEGENDSIZE=13
AXESSIZE=15
TITLESIZE=AXESSIZE

def myCmp(a, b, c=1):
    d=a<b
    return (d or c) and (not(d) or not(c)) #XOR

def timesF(arr1, arr2):
    ret=[]
    for a in arr1:
        for b in arr2:
            ret.append((a,b))
    return ret

def giveLabels(old, auto=True):
    newLabels=[n for n in old]
    ret=[]
    if auto:
        for l in newLabels:
            if l=="XavierFunc2":
                ret.append("Nvidia® Xavier GPU")
            elif l=="XavierCPUOpenMP":
                ret.append("Nvidia® Xavier CPU OpenMP")
            elif l=="XavierCPUNOOpenMP":
                ret.append("Nvidia® Xavier CPU no OpenMP")
            
            elif l=="NUCi7NOOpenMP":
                ret.append("Intel® NUC7i7BNH no OpenMP")
            elif l=="NUCi7OpenMP":
                ret.append("Intel® NUC7i7BNH OpenMP")

            elif l=="NUCAtomOpenMP":
                ret.append("Bqeel Z83-II OpenMP")

            elif l=="TX2Func1":
                ret.append("Nvidia® Jetson TX2 GPU")
            elif l=="TX2CPUOpenMP":
                ret.append("Nvidia® Jetson TX2 CPU OpenMP")
            elif l=="TX2Func2":
                ret.append("Nvidia® Jetson TX2 GPU function 2")
    else:
        for i in range(len(newLabels)):
            l=input(newLabels[i]+" -> ")
            ret.append(l)
    return ret

class Run:
    def __init__(self,name, test_name, discr, time, length, err, refinements, threads, functionType, jump, guessInitialAngles, initTime=-1.0, endTime=-1.0, power_file=""):
        self.name=name
        self.devName=name
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
        return ("dev_name: "+str(self.devName)+
                " test_name: "+str(self.test_name)+
                " discr: "+str(self.discr)+
                " refinements: "+str(self.refinements)+
                " threads: "+str(self.threads)+
                " guessInitialAngles: "+str(self.guessInitialAngles)+
                " functionType: "+str(self.functionType)+
                " jump: "+str(self.jump))
        #return (self.name+" "+self.test_name+" "+str(self.discr)+" "+str(self.time)+" "+self.power_file)

import re
def scan(r): #Consider the value taken each 50ms
    power=[]
    mean_power=0.0
    frequencies=[]
    t=0
    file_name=r.power_file
    with open(file_name) as file:
        for l in file:
            #https://docs.nvidia.com/jetson/archives/l4t-archived/l4t-3231/index.html#page/Tegra%20Linux%20Driver%20Package%20Development%20Guide/AppendixTegraStats.html
            if re.search("GR3D_FREQ", l)!=None:
                #Find values for tegrastats
                try: #Find values for power consumption
                    p_cpu=re.search("CPU (\d+)/(\d+)", l).groups()
                    p_gpu=re.search("GPU (\d+)/(\d+)", l).groups()
                    inst_cpu=float(int(p_cpu[0])/1000.0)
                    inst_gpu=float(int(p_gpu[0])/1000.0)
                    mean_power+=float(int(p_cpu[1])/1000.0)
                    mean_power+=float(int(p_gpu[1])/1000.0)
                    power.append({"W" : (inst_cpu+inst_gpu), "time" : t})
                except:
                    power.append({"W" : 0, "time" : t})    
                try: #Find values for frequencies
                    frequencies.append({"MHz" : int(re.search("GR3D_FREQ (\d+)%", l).groups()[0]), "time" : t})
                except: 
                    frequencies.append({"MHz" : int(re.search("GR3D_FREQ (\d+)%", l).groups()[0]), "time" : t})
            else:
                #Find values for NUC
                try:
                    # app=re.search("Power: (.+)", l).groups()[0].split(', ')
                    app=[float(w.split(': ')[1]) for w in re.search("Power: (.+), ", l).groups()[0].split(', ')]
                    p=0.0
                    for a in app:
                        p+=a
                    power.append({"W" : p, "time" : t})
                except Exception as e:
                    print(e)
                    power.append({"W" : 0, "time" : t})
            t+=50

    return (power, frequencies)


def getOptions(args=sys.argv[1:]):
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, default="tests.json", help="Input file. Default is tests.json")
    parser.add_argument("-n", "--names", type=str, default="", help="List of devices to print separated by comas without spaces.")
    parser.add_argument("-N", "--test_names", type=str, default="", help="List of the names of the tests to print separated by comas without spaces.")
    parser.add_argument("-r", "--ref", type=str, default="", help="List of refinements to consider. The list is coma-separated without spaces.")
    parser.add_argument("-f", "--func", type=str, default="2", help="List of functions to consider. The list is coma-separated without spaces. The default is 2.")
    parser.add_argument("-j", "--jump", type=str, default="3", help="List of jumps to consider. Taken into account only if 2 is present in the list of functions. The list is coma-separated without spaces. The default is 2.")
    parser.add_argument("-g", "--guess", type=int, default=1, help="Considered only if --func is not 0. It says whether to look for specific angles (1) or not (0). Default is 1.")
    parser.add_argument("-l", "--log", action="store_true", default=False, help="Whether to apply logarithm to data or not.")
    parser.add_argument("-M", "--max", type=int, default=sys.maxsize, help="Only times below this value are kept into consideration. Default is sys.maxsize.")
    parser.add_argument("-d", "--min_discr", type=int, default=0, help="The minimum number of discretizations to consider. Default is 0.")
    parser.add_argument("-D", "--max_discr", type=int, default=sys.maxsize, help="The maximum number of discretizations to consider. Default is sys.maxsize.")
    parser.add_argument("-t", "--min_thread", type=int, default=128, help="The minimum number of threads to consider. Default is 128.")
    parser.add_argument("-T", "--max_thread", type=int, default=1024, help="The maximum number of discretizations to consider. Default is 1024.")
    options=parser.parse_args(args)
    return options

def main():
    arguments=getOptions(sys.argv[1:])
    print(arguments)
    input_file="tests.json"
    acceptAllDevices=False
    acceptAllTests=False
    acceptAllRef=False
    names=[]
    testNames=[]
    defRefinements=[]
    log_b=arguments.log
    timeLimit=sys.maxsize
    defDiscr=[0, sys.maxsize]
    defThreads=[128, 1024]
    defFuncs=[]
    defJumps=[]
    defGuessed=True

    if arguments.input!="" and arguments.input.endswith(".json"):
        input_file=arguments.input
    try:
        names=arguments.names.split(',')
    except:
        pass
    if names==[] or names[0]=="":
        acceptAllDevices=True
        names=[]
    try:
        testNames=arguments.test_names.split(',')
    except:
        pass
    if testNames==[] or testNames[0]=="":
        acceptAllTests=True
        testNames=[]
    try:
        defRefinements.append(int(ref) for ref in arguments.refinements.split(','))
    except:
        pass
    if defRefinements==[] or defRefinements[0]=="":
        acceptAllRef=True
        defRefinements=[]
    try:
        timeLimit=int(arguments.max)
    except:
        pass
    try:
        defDiscr[0]=int(arguments.min_discr)
        defDiscr[1]=int(arguments.max_discr)
    except:
        pass
    try:
        defThreads[0]=int(arguments.min_thread)
        defThreads[1]=int(arguments.max_thread)
    except:
        pass
    
    try:
        print("arguments.func: "+str(arguments.func))
        for func in arguments.func.split(','):
            defFuncs.append(int(func))
    except:
        pass
    if defFuncs==[] or defFuncs[0]=="":
        defFuncs=[2]
    try:
        if 2 in defFuncs:
            try:
                for jump in arguments.jump.split(','):
                    defJumps.append(int(jump))
            except:
                pass
            if defJumps==[] or defJumps[0]=="":
                defJumps=[2]
    except:
        pass
    try:
        if defFunc!=0:
            defGuessed=bool(arguments.guess)
    except:
        pass

    print("input_file: "+str(input_file))
    print("acceptAllDevices: "+str(acceptAllDevices))
    print("acceptAllTests: "+str(acceptAllTests))
    print("acceptAllRef: "+str(acceptAllRef))
    print("names: "+str(names))
    print("testNames: "+str(testNames))
    print("defRefinements: "+str(defRefinements))
    print("log_b: "+str(log_b))
    print("timeLimit: "+str(timeLimit))
    print("defDiscr: "+str(defDiscr))
    print("defThreads: "+str(defThreads))
    print("defFuncs: "+str(defFuncs))
    print("defJumps: "+str(defJumps))
    print("defGuessed: "+str(defGuessed))

    #Read data from json file and add every possible thing
    runs=[]
    lines=0
    with open(input_file) as json_file:
        for l in json_file:
            lines+=1
    n_lines=0
    with open(input_file) as json_file:
        data=json.load(json_file)
        for p in data['run']:
            print("Reading file... {:.2f}%".format(n_lines/lines*100.0), end="\r")
            n_lines+=1
            r=None
            initialTime=-1
            finalTime=-1
            try:
                initialTime=p['initTime']
            except:
                pass
            try:
                finalTime=p['endTime']
            except:
                pass
            
            try:
                r=Run(p['name'], p['test_name'], p['discr'], p['time'], p['length'], p['err'], p['refinements'], p['threads'], p['functionType'], p['jump'], p['guessInitialAngles'], initialTime, finalTime, p['power_file'])
            except:
                r=Run(p['name'], p['test_name'], p['discr'], p['time'], p['length'], p['err'], p['refinements'], p['threads'], p['functionType'], p['jump'], p['guessInitialAngles'], initialTime, finalTime)
            
            realName=r.name.split('_')[0]

            #Check if CPU: threads==0, functype==0, jump==0 guessInitialAngles==False
            if (acceptAllDevices or (not acceptAllDevices and (realName!="" or realName!=" "))):
                #print("acceptAllDevices")
                if (acceptAllTests or (not acceptAllTests and r.test_name in testNames)):
                    #print("acceptAllTests")
                    if (acceptAllRef or (not acceptAllRef and r.refinements in defRefinements)):
                        #print("acceptAllRef")
                        if r.time<timeLimit:
                            #print("time")
                            if  r.discr>=defDiscr[0] and r.discr<=defDiscr[1]:
                                if r.threads==0 and r.jump==0 and r.functionType==0 and r.guessInitialAngles==False:
                                    found=False
                                    r.devName=realName
                                    for i in range(len(runs)):
                                        if runs[i]==r:
                                            runs[i].time+=r.time 
                                            runs[i].n+=1
                                            found=True
                                    if not found:
                                        runs.append(r)


            if (acceptAllDevices or (not acceptAllDevices and (realName!="" or realName!=" "))):
                #print("acceptAllDevices")
                if (acceptAllTests or (not acceptAllTests and r.test_name in testNames)):
                    #print("acceptAllTests")
                    if (acceptAllRef or (not acceptAllRef and r.refinements in defRefinements)):
                        #print("acceptAllRef")
                        if r.time<timeLimit:
                            #print("time")
                            if  r.discr>=defDiscr[0] and r.discr<=defDiscr[1]:
                                #print("discr")
                                if  r.threads>=defThreads[0] and r.threads<=defThreads[1]:
                                    #print("threads")
                                    #print("functionType: "+str(r.functionType)+" "+str(type(r.functionType))+" "+str(defFuncs))
                                    if  r.functionType in defFuncs:
                                        #print("jump "+str(r.jump)+" "+str(type(r.jump))+" "+str(defJumps))
                                        if (r.functionType!=2 or (r.functionType==2 and r.jump in defJumps)):
                                            if  r.guessInitialAngles==defGuessed:
                                                #print("guessInitialAngles")
                                                #if r.power_file!="":
                                                    #(p, f)=scan(r)
                                                    #r.power_cons.append(p)
                                                #Run through all runs and check if there is an equal one
                                                found=False
                                                r.devName=realName
                                                for i in range(len(runs)):
                                                    if runs[i]==r:
                                                        runs[i].time+=r.time 
                                                        runs[i].n+=1
                                                        #runs[i].power_cons.append(r.power_cons[0])
                                                        found=True
                                    
                                                if not found:
                                                    runs.append(r)
    print("                                                 \rReading file... 100%")
    #Run through all runs and check which had multiple times.
    i=0
    for r in runs:
        if r.n>1:
            print(r.n, r.name)
            r.time=(r.time/r.n)
            maxL=0
            for i in range(len(r.power_cons)): #Find the longest log
                if len(r.power_cons[i])>maxL:
                    maxL=len(r.power_cons[i])
            mean_ps=[]
            for i in range(maxL):
                mean_p=0
                nArrOk=0 #Divide just by the number of arrays large at least i
                for arr in range(len(r.power_cons)):
                    if i<len(r.power_cons[arr]):
                        mean_p+=r.power_cons[arr][i]["W"]
                        nArrOk+=1
                if i==0 and len(r.power_cons)!=nArrOk:
                    print(i, len(r.power_cons), nArrOk)
                mean_p/=nArrOk
                mean_ps.append({'W' : mean_p, 'time' : i*50})
            r.power_cons=mean_ps
        i+=1

    discrs=[]
    threads=[]
    funcs=[]
    jumps=[]
    for r in runs:
        if acceptAllRef:
            if r.refinements not in defRefinements:
                defRefinements.append(r.refinements)
        if r.discr not in discrs:
            discrs.append(r.discr)
        if r.threads not in threads:
            threads.append(r.threads)
        if r.functionType not in funcs:
            funcs.append(r.functionType)
        if r.functionType==2 and r.jump not in jumps:
            jumps.append(r.jump)
        if r.devName not in names and acceptAllDevices:
            names.append(r.devName)
        if r.test_name not in testNames and acceptAllTests:
            testNames.append(r.test_name)
    refs=defRefinements
    
    print("discrs: "+str(discrs))
    print("threads: "+str(threads))
    print("funcs: "+str(funcs))
    print("jumps: "+str(jumps))
    print("names: "+str(names))
    print("testNames: "+str(testNames))

    import matplotlib.colors as mcolors
    colors=['blue', 'orange', 'red', mcolors.CSS4_COLORS['limegreen'], mcolors.CSS4_COLORS['darkgreen'], 'purple', 'red', 'gray']

    choice=input("Do you want to show the diagram over threads and different functions? [Y/n] ")
    if choice=="Y" or choice=="y":
        for tn in range(len(testNames)): 
            #break
            width=0.4
            for th in threads:
                nGraphs=0
                fig=plt.figure(figsize=(6.4*3, 4.8*3), constrained_layout=True)
                spec=fig.add_gridspec(int(((len(funcs)+len(jumps)))/2), 2)
    #           fig, ax=plt.subplots()
                ax=[]
                valuess=[]
                for j in jumps:
                    print("j: "+str(j))
                    for fn in funcs: #Create a graph for each function
                        if fn!=2 and j!=jumps[-1]:
                            continue
                        print("fn: "+str(fn))
                        ax.append(fig.add_subplot(spec[int(nGraphs/2), nGraphs%2]))
                        for n in range(len(names)):
                            values=[]
                            labels=timesF(discrs, defRefinements)
                            x=np.arange(len(labels))
                            for l in range(len(labels)):
                                (dis,ref)=labels[l]
                                #print("looking for: "+names[n]+" with discr: "+str(dis)+", ref: "+str(ref))
                                before=len(values)
                                for r in runs:
                                    #print(r)
                                    if r.devName==names[n]:
                                        #print("r.devName")
                                        #print("r.test_name <"+r.test_name+"> <"+testNames[tn]+"> "+str(len(testNames)))
                                        if  r.test_name==testNames[tn]: 
                                            if r.discr==dis: 
                                                #print("r.discr")
                                                if  r.refinements==ref: 
                                                    #print(" r.refinements")
                                                    if r.threads==th: 
                                                        #print("r.threads")
                                                        if  r.functionType==fn: 
                                                            #print(" r.functionType")
                                                            if fn!=2 or (fn==2 and r.jump==j):
                                                                #print("r.jump")
                                                                values.append(r)
                                #if len(values)==before:
                                #    print("not found")
                                #else:
                                #    print("found: "+str(values[-1]))
                            if values!=[]:
                                valuess.append(values)
                                #print(values)
                                ax[nGraphs].barh(x+n*width-width/len(names), [v.time for v in values], width, color=colors[n], label=names[n]) #color=colors[n]
                                ax[nGraphs].set_xlabel('Times in ms')
                                ax[nGraphs].set_yticks(x)
                                ax[nGraphs].set_yticklabels(labels)
                                if func!=2:
                                    ax[nGraphs].set_title(testNames[tn]+", threads: "+str(th)+" func: "+str(fn)+", jump: 0", fontsize=TITLESIZE)
                                else:
                                    ax[nGraphs].set_title(testNames[tn]+", threads: "+str(th)+" func: "+str(fn)+", jump: "+str(j), fontsize=TITLESIZE)
                                ax[nGraphs].legend(fontsize=LEGENDSIZE)
                        nGraphs+=1
                for vv in valuess[0:2]:
                    print("<", end="")
                    for v in vv:
                        print(v, end=", ")
                    print(">")
                #ax.legend(fontsize=LEGENDSIZE)
                plt.show()

    choice=input("Do you want to show the diagrams of times for each example and discretization without OpenMP? [Y/n] ")
    if choice=="Y" or choice=="y" or choice=="":
        print(names)
        Discrs=[d for d in discrs if d in [90,360,720]]
        Refs=[r for r in defRefinements if r in [1,4,16]]
        #Names=[n for n in names if ("NOOpenMP" in n or "XavierFunc2" in n or "TX2Func1" in n)]
        Names=[]
        for n in names:
            if "NOOpenMP" in n or "XavierFunc2" in n or "TX2Func1" in n:
                print("ok: "+n)
                Names.append(n)
            else:
                print("no: "+n+" "+str("NOOpenMP" in n))
        newLabels=[]
        print(Names)
        nGraphs=0
        for tn in range(len(testNames)):
            if testNames[tn]!="Kaya Example 4" and testNames[tn]!="Omega":
                continue
            else: 
                for dis in Discrs:
                    print(testNames[tn])
                    jumpD=[3,0]
                    funcD=[2,1,0]
                    threadD=[128,0]

                    if newLabels==[]:
                        print(Names)
                        choiceC=input("These are the labels in list, do you want to choose new labels? [N/y] ")
                        if choiceC=="Y" or choiceC=="y":
                            newLabels=giveLabels([n for n in Names])
                        print(newLabels)

                    width=0.1
                    fig, ax=plt.subplots()
                    for n in range(len(Names)):
                        values=[]
                        x=np.arange(len(Refs))
                        for ref in Refs:
                            for r in runs:
                                if  r.devName==Names[n] and r.test_name==testNames[tn] and r.discr==dis and \
                                    r.refinements==ref and r.threads in threadD and r.functionType in funcD and \
                                    (r.functionType!=2 or (r.functionType==2 and r.jump in jumpD)):
                                    #print(r)
                                    values.append(r.time/1000.0)
                        a=1
                        if len(newLabels)%2!=1:
                            a=0.5
                        if n%2==0:
                            ax.bar(x-width*(a+int(n/2)), values, width, label=newLabels[n], color=colors[n])
                        else:
                            ax.bar(x+width*(a+int(n/2)), values, width, label=newLabels[n], color=colors[n])

                    ax.set_title(testNames[tn]+", discretization: "+str(dis), fontsize=TITLESIZE)
                    ax.set_xlabel("Refinements $m$", fontsize=AXESSIZE)
                    ax.set_xticklabels([str(r) for r in Refs])
                    ax.set_ylabel("Time (s)", fontsize=AXESSIZE)
                    ax.set_xticks(x)
                    ax.legend(fontsize=LEGENDSIZE)
                    plt.savefig("images/timesGraphsSeparatedNOOpenMP"+str(nGraphs)+"s", format="pgf", transparent=True, bbox_inches='tight', pad_inches=0.1)
                    nGraphs+=1

    choice=input("Do you want to show the diagrams of times for each example and discretization? [Y/n] ")
    if choice=="Y" or choice=="y" or choice=="":
        print(names)
        Discrs=[d for d in discrs if d in [90,360,720]]
        Refs=[r for r in defRefinements if r in [1,4,16]]
        Names=[n for n in names if "NOOpenMP" not in n and "TX2Func2" not in n]
        print(Names)
        nGraphs=0
        newLabels=[]
        for tn in range(len(testNames)):
            if testNames[tn]=="Kaya Example 1" or testNames[tn]=="Kaya Example 2":
                continue
            else: 
                for dis in Discrs:
                    print(testNames[tn])
                    jumpD=[3,0]
                    funcD=[2,1,0]
                    threadD=[128,0]

                    width=0.1
                    fig, ax=plt.subplots()
                    valuess=[]
                    if newLabels==[]:
                        print(Names)
                        choiceC=input("These are the labels in list, do you want to choose new labels? [N/y] ")
                        if choiceC=="Y" or choiceC=="y":
                            newLabels=giveLabels([n for n in Names])
                            
                    for n in range(len(Names)):
                        values=[]
                        x=np.arange(len(Refs))
                        for ref in Refs:
                            #print("looking for: "+names[n]+" with dis: "+str(dis)+", ref: "+str(ref))
                            for r in runs:
                                #print(r)
                                if r.devName==Names[n]:
                                    #print("r.devName")
                                    #print("r.test_name <"+r.test_name+"> <"+testNames[tn]+"> "+str(len(testNames)))
                                    if  r.test_name==testNames[tn]: 
                                        if r.discr==dis: 
                                            #print("r.dis")
                                            if r.refinements==ref: 
                                                #print(" r.refinements")
                                                if r.threads in threadD: 
                                                    #print("r.threads")
                                                    if  r.functionType in funcD: 
                                                        #print(" r.functionType")
                                                        if r.functionType!=2 or (r.functionType==2 and r.jump in jumpD):
                                                            #print("r.jump")
                                                            values.append(r.time/1000.0)
                  
                        print(values)
                        print(x)
                        a=1
                        if len(newLabels)%2!=1:
                            a=0.5
                        print(a)
                        if n%2==0:
                            ax.bar(x-width*(a+int(n/2)), values, width, label=newLabels[n], color=colors[n])
                        else:
                            ax.bar(x+width*(a+int(n/2)), values, width, label=newLabels[n], color=colors[n])
                                        
                    ax.set_title(testNames[tn]+", discretization: "+str(dis), fontsize=TITLESIZE)
                    ax.set_xlabel("Refinements $m$", fontsize=AXESSIZE)
                    ax.set_xticklabels([str(r) for r in Refs])
                    ax.set_ylabel("Time (s)", fontsize=AXESSIZE)
                    ax.set_xticks(x)
                    ax.tick_params(labelsize=AXESSIZE)
                    ax.legend(fontsize=LEGENDSIZE)
                    plt.savefig("images/timesGraphsSeparated"+str(nGraphs)+"s.pdf", transparent=True, bbox_inches='tight', pad_inches=0.1)
                    nGraphs+=1



    choice=input("Do you want to show the diagrams of times for each example? [Y/n] ")
    if choice=="Y" or choice=="y" or choice=="":
        print(names)
        Discrs=[d for d in discrs if d in [90,360,720]]
        Refs=[r for r in defRefinements if r in [1,4,16]]
        Names=[n for n in names if "NOOpenMP" not in n]
        print(Names)
        nGraphs=0
        for tn in range(len(testNames)):
            if testNames[tn]=="Kaya Example 1" or testNames[tn]=="Kaya Example 2":
                continue
            else: 
                print(testNames[tn])
                jumpD=[3,0]
                funcD=[2,0]
                threadD=[128,0]

                width=0.1
                fig, ax=plt.subplots()
                valuess=[]
                for n in range(len(Names)):
                    values=[]
                    labels=timesF(Discrs, Refs)
                    x=np.arange(len(labels))
                    for l in range(len(labels)):
                        (dis,ref)=labels[l]
                        #print("looking for: "+names[n]+" with discr: "+str(dis)+", ref: "+str(ref))
                        for r in runs:
                            #print(r)
                            if r.devName==Names[n]:
                                #print("r.devName")
                                #print("r.test_name <"+r.test_name+"> <"+testNames[tn]+"> "+str(len(testNames)))
                                if  r.test_name==testNames[tn]: 
                                    if r.discr==dis: 
                                        #print("r.discr")
                                        if r.refinements==ref: 
                                            #print(" r.refinements")
                                            if r.threads in threadD: 
                                                #print("r.threads")
                                                if  r.functionType in funcD: 
                                                    #print(" r.functionType")
                                                    if r.functionType!=2 or (r.functionType==2 and r.jump in jumpD):
                                                        #print("r.jump")
                                                        values.append(log(r.time))
              
                    print(values)
                    print(labels)
                    print(x)
                    a=1
                    if len(names)%2!=1:
                        a=0.5
                    if n%2==0:
                        ax.bar(x-width*(a+int(n/2)), values, width, label=Names[n], color=colors[n])
                    else:
                        ax.bar(x+width*(a+int(n/2)), values, width, label=Names[n], color=colors[n])
                
                ax.set_title(testNames[tn], fontsize=TITLESIZE)
                ax.set_xlabel("Disc $k$, ref $m$", fontsize=AXESSIZE)
                ax.set_xticklabels([0]+[str((d,r)) for (d,r) in labels ])
                ax.set_ylabel("Logarithmic times", fontsize=AXESSIZE)
                ax.legend(fontsize=LEGENDSIZE)
                plt.savefig("images/timesGraphsAllInOneLog"+str(nGraphs)+".pdf", transparent=True, bbox_inches='tight', pad_inches=0.1)
                nGraphs+=1



    choice=input("Do you want to print the tables? [Y/n] ")
    if choice=="Y" or choice=="y" or choice=="":
        lthreads=input("For which number of threads? ")
        lfunc=2
        ljump=3
        lguess=True

        fileName=input("Do you wish to save the data to file? [N/<name_of_file>] ")
        origin_stdout=sys.stdout
        f=None
        if input!="N" and input!="n" and input!="":
            f=open(fileName, "a+")
            sys.stdout=f

        #for dv1 in names:
        for tn in testNames:
            print("\\subsection{"+tn+"}")
            print("\\begin{center}")
            print("\\begin{tabular}{c|c|c", end="")
            for n in names:
                print("|c", end="")
            print("}")
            print("\t&&&\\multicolumn{"+str(len(names))+"}{c}{Times(ms)}\\\\")
            print("\tDisc&Ref&Err", end="")
            for n in names:
                print("&"+n, end="")
            print("\\\\\n\\hline")
            for discr in discrs:
                for ref in defRefinements:
                    for dv in range(len(names)): 
                        for r in runs:
                            if r.devName==names[dv] and r.test_name==tn and r.discr==discr and r.refinements==ref:
                                if  str(r.threads) in lthreads and r.functionType==lfunc and \
                                    r.jump==ljump and r.guessInitialAngles==lguess:
                                    if dv==0:
                                        print("{:d}&{:d}&{:.1e}&{:.1e}".format(discr, ref, r.err, r.time), end="")
                                    elif dv==len(names)-1:
                                        print("&{:.1e}\\\\".format(r.time))
                                    else:
                                        print("&{:.1e}".format(r.time), end="")


                                elif r.threads==0 and r.functionType==0 and r.jump==0 and \
                                    r.guessInitialAngles==False:
                                    if dv==0:
                                        print("{:d}&{:d}&{:.1e}&{:.1e}".format(discr, ref, r.err, r.time), end="")
                                    elif dv==len(names)-1:
                                        print("&{:.1e}\\\\".format(r.time))
                                    else:
                                        print("&{:.1e}".format(r.time), end="")
                print("\\hline")
            print("\\end{tabular}")
            print("\\end{center}\n\n\n\n\n\n\n")

        if input!="N" and input!="n" and input!="":
            sys.stdout=origin_stdout
            f.close()
    

    choice=input("Do you want to show the power-consumption graphs? [Y/n] ")
    if choice=="Y" or choice=="y" or choice=="":
        print(names)
        lNamesC=input("These are the devices in list, do you want to show all of them? Otherwise write a coma-separated list [Y/list] ")
        lNames=[]
        if lNamesC=="Y" or lNamesC=="y" or lNamesC=="":
            lNames=names
        else:
            lNames=lNamesC.split(",")
        print("lNames: "+str(lNames))
        print(discrs)
        ldiscrsC=input("These are the discriminations in list, do you want to show all of them? Otherwise write a coma-separated list [Y/list] ")
        ldiscrs=[]
        if ldiscrsC=="Y" or ldiscrsC=="y" or ldiscrsC=="":
        	ldiscrs=discrs
        else:
        	ldiscrs=[int(d) for d in ldiscrsC.split(",")]
        print("ldiscrs "+str(ldiscrs))

        print(refs)
        lrefsC=input("These are the refinements in list, do you want to show all of them? Otherwise write a coma-separated list [Y/list] ")
        lrefs=[]
        if lrefsC=="Y" or lrefsC=="y" or lrefsC=="":
            lrefs=refs
        else:
            lrefs=[int(d) for d in lrefsC.split(",")]
        print("lrefs "+str(lrefs))

        print(threads)
        lthreadsC=input("These are the threads in list, do you want to show all of them? Otherwise write a coma-separated list [Y/list] ")
        lthreads=[]
        if lthreadsC=="Y" or lthreadsC=="y" or lthreadsC=="":
            lthreads=threads
        else:
            lthreads=lthreadsC.split(",")

        print(funcs)
        lfuncsC=input("These are the function types in list, do you want to show all of them? Otherwise write a coma-separated list [Y/list] ")
        lfuncs=[]
        if lfuncsC=="Y" or lfuncsC=="y" or lfuncsC=="":
            lfuncs=funcs
        else:
            lfuncs=lfuncsC.split(",")

        print(jumps)
        ljumpsC=input("These are the jumps in list, do you want to show all of them? Otherwise write a coma-separated list [Y/list] ")
        ljumps=[]
        if ljumpsC=="Y" or ljumpsC=="y" or ljumpsC=="":
            ljumps=jumps
        else:
            ljumps=ljumpsC.split(",")

        print(testNames)
        ltestNameC=input("These are the test names in list. Default is [\"Circuit\"], do you want to choose another? Otherwise write a coma-separated list [N/list] ")
        ltestName=["Circuit"]
        if ltestNameC!="N" and ltestNameC!="n" and ltestNameC!="":
            ltestName=ltestNameC.split(",")

        choice=input("Do you want to show the figures or to save it? [Y/<file_name>.pdf] ")
        filename=""
        if choice!="Y" and choice!="y" and choice!="":
            filename=choice

        newLabels=[]
        for tn in range(len(ltestName)):
            values=[]
            times=[]
            labels=[]
            maxV=0
            
            for r in runs:
                if r.devName in lNames and r.refinements in lrefs and r.discr in ldiscrs:
                    goOn=False
                    #If it's CPU then only the name
                    if r.threads==0 and r.functionType==0 and r.guessInitialAngles==False and r.jump==0:
                        goOn=True
                    #Else also threads, jumps and functiontype
                    else:
                        if r.threads in lthreads and r.functionType in lfuncs and (r.functionType==1 and r.jump==0 or r.jump in ljumps):
                            goOn=True
                    if goOn:
                        if r.test_name==ltestName[tn]:
                            print(r)
                            print(r.initTime, r.endTime)
                            powerFile=r.power_file
                            watts=[]
                            ts=[]

                            with open(powerFile, "r") as powerJson:
                                data=json.load(powerJson)
                                for l in data['power']:
                                    if l['time']>r.initTime and l['time']<r.endTime:
                                        watts.append(l['power'])
                                        ts.append((l['time']-r.initTime)/1000.0)
                                        if l['power']>maxV:
                                            maxV=l['power']
                            values.append(watts)
                            times.append(ts)
                            labels.append(r.devName)
            if newLabels==[]:
                print(labels)
                choiceC=input("These are the labels in list, do you want to choose new labels? [N/y] ")
                if choiceC=="Y" or choiceC=="y":
                    newLabels=giveLabels([l for l in labels])
                else:
                    newLabels=labels
                    
            fig, ax=plt.subplots()
            print(len(times))
            for i in range(len(times)):
                ax.plot(times[i], values[i], color=colors[i], label=newLabels[i])

            ax.set_title(ltestName[tn]+", k=360, m=16", fontsize=TITLESIZE)
            ax.set_xlabel("Time (s)", fontsize=AXESSIZE)
            ax.set_ylabel("Power (W)", fontsize=AXESSIZE)
            #ax.ylim(0, 13/9*maxV)
            ax.legend(fontsize=LEGENDSIZE, loc="upper right")
            if filename!="":
                fig.savefig("images/"+filename+str(tn)+".pdf", transparent=True, bbox_inches='tight', pad_inches=0.1)
            else:
                fig.show()

    choice=input("Do you want to show the graph for normalized power consumption? [Y/n] ")
    if choice=="Y" or choice=="y" or choice=="":
        print(names)
        lNamesC=input("These are the devices in list, do you want to show all of them? Otherwise write a coma-separated list [Y/list] ")
        lNames=[]
        if lNamesC=="Y" or lNamesC=="y" or lNamesC=="":
            lNames=names
        else:
            lNames=lNamesC.split(",")
        print("lNames: "+str(lNames))
        print(discrs)
        ldiscrsC=input("These are the discriminations in list, do you want to show all of them? Otherwise write a coma-separated list [Y/list] ")
        ldiscrs=[]
        if ldiscrsC=="Y" or ldiscrsC=="y" or ldiscrsC=="":
            ldiscrs=discrs
        else:
            ldiscrs=[int(d) for d in ldiscrsC.split(",")]
        print("ldiscrs "+str(ldiscrs))

        print(refs)
        lrefsC=input("These are the refinements in list, do you want to show all of them? Otherwise write a coma-separated list [Y/list] ")
        lrefs=[]
        if lrefsC=="Y" or lrefsC=="y" or lrefsC=="":
            lrefs=refs
        else:
            lrefs=[int(d) for d in lrefsC.split(",")]
        print("lrefs "+str(lrefs))

        print(threads)
        lthreadsC=input("These are the threads in list, do you want to show all of them? Otherwise write a coma-separated list [Y/list] ")
        lthreads=[]
        if lthreadsC=="Y" or lthreadsC=="y" or lthreadsC=="":
            lthreads=threads
        else:
            lthreads=lthreadsC.split(",")

        print(funcs)
        lfuncsC=input("These are the function types in list, do you want to show all of them? Otherwise write a coma-separated list [Y/list] ")
        lfuncs=[]
        if lfuncsC=="Y" or lfuncsC=="y" or lfuncsC=="":
            lfuncs=funcs
        else:
            lfuncs=lfuncsC.split(",")

        print(jumps)
        ljumpsC=input("These are the jumps in list, do you want to show all of them? Otherwise write a coma-separated list [Y/list] ")
        ljumps=[]
        if ljumpsC=="Y" or ljumpsC=="y" or ljumpsC=="":
            ljumps=jumps
        else:
            ljumps=ljumpsC.split(",")

        print(testNames)
        ltestNameC=input("These are the test names in list. Default is [\"Circuit\"], do you want to choose another? Otherwise write a coma-separated list [N/list] ")
        ltestName=["Circuit"]
        if ltestNameC!="N" and ltestNameC!="n" and ltestNameC!="":
            ltestName=ltestNameC.split(",")

        choice=input("Do you want to show the figures or to save it? [Y/<file_name>.pdf] ")
        filename=""
        if choice!="Y" and choice!="y" and choice!="":
            filename=choice

        newLabels=[]
        results=[]
        labels=[]
        maxV=0
        for dni in range(len(lNames)):
            dn=lNames[dni]
            results.append([])
            for tn in range(len(ltestName)):
                watts=[]
                timeM=0
                
                for r in runs:
                    if r.devName==dn and r.refinements in lrefs and r.discr in ldiscrs:
                        goOn=False
                        #If it's CPU then check only the name
                        if r.threads==0 and r.functionType==0 and r.guessInitialAngles==False and r.jump==0:
                            goOn=True
                        #Else check also threads, jumps and functiontype
                        elif r.threads in lthreads and r.functionType in lfuncs and (r.functionType==1 and r.jump==0 or r.jump in ljumps):
                                goOn=True
                        if goOn:
                            if r.test_name==ltestName[tn]:
                                powerFile=r.power_file
                                times=[]

                                with open(powerFile, "r") as powerJson:
                                    data=json.load(powerJson)
                                    for l in data['power']:
                                        if l['time']>r.initTime and l['time']<r.endTime:
                                            times.append((l['time']-r.initTime)/1000.0)
                                            if len(times)>1:
                                                watts.append(l['power']*(times[-1]-times[-2]))
                                            else:
                                                watts.append(l['power']*times[-1])

                                            if l['power']>maxV:
                                                maxV=l['power']

                                timeM=times[-1]
                                if r.devName not in labels:
                                    labels.append(r.devName)

                sum=0.0
                for v in watts:
                    sum+=v
                results[-1].append(sum/timeM)



        if newLabels==[]:
            print(labels)
            choiceC=input("These are the labels in list, do you want to choose new labels? [N/y] ")
            if choiceC=="Y" or choiceC=="y":
                newLabels=giveLabels([l for l in labels])
            else:
                newLabels=labels
        
        final=[]
        for resi in range(len(results)):
            res=results[resi]
            print(newLabels[resi],res)
            sum=0.0
            for r in res:
                sum+=r
            final.append(sum/len(res))

        fig, ax=plt.subplots()
        x=np.arange(len(newLabels))

        dics=[]
        for i in range(len(final)):
            dics.append({"final" : final[i], "color" : colors[i], "label" : newLabels[i]})
        dics.sort(key=lambda x : x['final'], reverse=True)

        for i in range(len(dics)):
            print(str(dics[i]))
            ax.bar(i, dics[i]['final'], width=0.5, color=dics[i]['color'], label=dics[i]['label'])
        ax.set_title("Normalised power consumption over seconds", fontsize=TITLESIZE)
        ax.set_xlabel("Devices", fontsize=AXESSIZE)
        ax.set_ylabel("$\\left({}^W/_s\\right)$", fontsize=AXESSIZE)
        ax.set_xticklabels([])
        ax.set_xticks([])
        ax.legend(fontsize=LEGENDSIZE, loc="upper right")
        #ax.ylim(0, 13/9*maxV)
        if filename!="":
            fig.savefig("images/"+filename+".pdf", transparent=True, bbox_inches='tight', pad_inches=0.1)
        else:
            fig.show()



    choice=input("Do you want to print the table of errors? [Y/n] ")
    if choice=="Y" or choice=="y" or choice=="":
        fileName=input("Do you wish to save the data to file? [N/<name_of_file>] ")
        origin_stdout=sys.stdout
        f=None
        if input!="N" and input!="n" and input!="":
            f=open(fileName, "a+")
            sys.stdout=f
        
        lfunc=2
        ljump=3
        lguess=True
        for tn in testNames:
            print("\\subsection{"+tn+"}")
            print("\\begin{center}")
            print("\\begin{tabular}{c|c", end="")
            for n in names:
                print("|c", end="")
            print("}")
            print("\t&&\\multicolumn{"+str(len(names))+"}{c}{Errors}\\\\")
            print("\tDisc&Ref", end="")
            for n in names:
                print("&"+n, end="")
            print("\\\\\n\\hline")
            for discr in discrs:
                for ref in defRefinements:
                    for dv in range(len(names)): 
                        for r in runs:
                            if r.devName==names[dv] and r.test_name==tn and r.discr==discr and r.refinements==ref:
                                if  r.threads in threads and r.functionType==lfunc and \
                                    r.jump==ljump and r.guessInitialAngles==lguess:
                                    if dv==0:
                                        print("{:d}&{:d}&{:.1e}".format(discr, ref, r.err), end="")
                                    elif dv==len(names)-1:
                                        print("&{:.1e}\\\\".format(r.err))
                                    else:
                                        print("&{:.1e}".format(r.err), end="")


                                elif r.threads==0 and r.functionType==0 and r.jump==0 and \
                                    r.guessInitialAngles==False:
                                    if dv==0:
                                        print("{:d}&{:d}&{:.1e}".format(discr, ref, r.err), end="")
                                    elif dv==len(names)-1:
                                        print("&{:.1e}\\\\".format(r.err))
                                    else:
                                        print("&{:.1e}".format(r.err), end="")
                print("\\hline")
            print("\\end{tabular}")
            print("\\end{center}\n\n\n\n\n\n\n")

        if input!="N" and input!="n" and input!="":
            sys.stdout=origin_stdout
            f.close()
   

    choice=input("\n\nDo you want to enter the SPE part? [N/y] ")
    if choice=="N" or choice=="n" or choice=="":
        sys.exit(0)
    else:
        import readSPE
        readSPE.main()



if __name__=="__main__":
    main()


