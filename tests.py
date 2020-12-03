import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from math import log
import argparse, sys

class Run:
    def __init__(self, name, discr, time, test_name, power_file=""):
        self.name=name
        self.discr=discr
        self.time=time
        self.test_name=test_name
        self.power_file=power_file

    def __str__(self):
        return (self.name+" "+self.test_name+" "+str(self.discr)+" "+str(self.time)+" "+self.power_file)

def getOptions(args=sys.argv[1:]):
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, default="tests.json", help="Input file. Default is tests.json")
    parser.add_argument("-n", "--names", type=str, default="", help="List of devices to print separated by comas without spaces.")
    parser.add_argument("-l", "--log", action="store_true", default=False, help="Whether to apply logarithm to data or not.")
    parser.add_argument("-m", "--max", type=int, default=sys.maxsize, help="Only times below this value are kept into consideration. Default is sys.maxsize.")
    parser.add_argument("-d", "--discr", type=int, default=0, help="The minimum number of discretizations to consider. Default is 0.")
    options=parser.parse_args(args)
    return options

import re
def scan(file_name):
    power=[]
    frequencies=[]
    times=[]
    t=0
    with open(file_name) as file:
        for l in file:
            power.append(float(int(re.search("GPU (\d+)/(\d+)", l).groups()[0])/1000))
            frequencies.append(int(re.search("GR3D_FREQ (\d+)%", l).groups()[0]))
            times.append(t)
            t+=50
    return (power, frequencies, times)


def main():
    arguments=getOptions(sys.argv[1:])
    input_file="tests.json"
    names=[]
    log_b=arguments.log
    limit=sys.maxsize
    discr=0
    print(arguments)
    if arguments.input!="" and arguments.input.endswith(".json"):
        input_file=arguments.input
    try:
        names=arguments.names.split(',')
    except:
        pass
    try:
        limit=int(arguments.max)
    except:
        pass
    try:
        discr=int(arguments.discr)
    except:
        pass


    #Read data from json file
    runs=[]
    with open(input_file) as json_file:
        data=json.load(json_file)
        for p in data['run']:
            try:
                runs.append(Run(p['name'], p['discr'], p['time'], p['test_name'], p['power_file']))
            except:
                runs.append(Run(p['name'], p['discr'], p['time'], p['test_name']))

    #Add data to structures
    tests=[]
    devices=[]
    discrs=[]
    for r in runs:
        if not r.test_name in tests:
            tests.append(r.test_name)
        if (not r.name in devices):
            if len(names)==0 or (len(names)!=0 and (r.name in names)):
                devices.append(r.name)
        if not r.discr in discrs and r.discr>discr:
            discrs.append(r.discr)
    discrs.sort()
    print(tests)
    print(devices)
    print(discrs)

    #Create matrix    
    matrix=[[[] for i in range(len(devices))] for j in range(len(tests))]
    for t in range(len(tests)):
        for i in range(len(devices)):
            for d in range(len(discrs)):
                for r in range(len(runs)):
                    if runs[r].discr==discrs[d] and runs[r].name==devices[i] and runs[r].test_name==tests[t]:
                        if runs[r].time<limit:
                            if log_b:
                                matrix[t][i].append(log(8*(runs[r].time),2))
                            else:
                                matrix[t][i].append(runs[r].time)
                        else: 
                            matrix[t][i].append(0)
        inv=[[0 for j in range(len(devices))] for i in range(len(discrs))]
        for i in range(len(matrix[t])):
            for j in range(len(matrix[t][i])):
                inv[j][i]=matrix[t][i][j]
        # for row in inv:
        #     print(row)
        # print()

    import matplotlib.colors as mcolors
    colors=['blue', 'orange', 'red', mcolors.CSS4_COLORS['limegreen'], mcolors.CSS4_COLORS['darkgreen'], 'purple', 'red', 'gray']

    for i in range(len(tests)):
        labels = discrs
        x = np.arange(len(labels))  # the label locations
        width = 0.1  # the width of the bars

        rects=[]

        fig=plt.figure(constrained_layout=True)
        spec=fig.add_gridspec((int(len(discrs)/2)+1), 2)
        ax=[]
        for d in range(len(discrs)+1):
            if d==0 and (len(discrs))%2==0:
                ax.append(fig.add_subplot(spec[0, :]))
            else:
                x_=int((d+1)/2)
                y_=1
                if d%2==1:
                    y_=0
                ax.append(fig.add_subplot(spec[x_, y_]))
        l=len(devices)
        for j in range(l):
            values=matrix[i][j]
            rect=ax[0].bar(x+(j-(l-1)/2)*width, values, width, label=devices[j], color=colors[j])
            rects.append(rect)

        ax[0].set_ylabel('Times in ms')
        ax[0].set_xlabel('#Discretizations')
        ax[0].set_title(tests[i])
        ax[0].set_xticks(x)
        ax[0].set_xticklabels(labels)
        ax[0].legend()

        longest_x=0.0
        longest_y=0.0
        for d in range(1, len(discrs)+1):
            ax[d].set_ylabel("Power (W)")
            ax[d].set_xlabel("Time in ms")
            ax[d].set_title(str(discrs[d-1])+" discretizations")

            logs=[]
            dev_names=[]
            for r in runs:
                if r.power_file!="" and (r.discr==discrs[d-1] and r.test_name==tests[i]):
                    logs.append(r.power_file)
                    dev_names.append(r.name)
            for l in range(len(logs)):
                p, f, t=scan(logs[l])
                ax[d].step(t, p, label=dev_names[l])
                if t[-1]>longest_x:
                    longest_x=t[-1]
                if max(p)>longest_y:
                    longest_y=max(p)
        for d in range(1, len(discrs)+1):
            ax[d].set_xticks(list(np.arange(0, longest_x, int(longest_x/10)))+[longest_x])
            ax[d].set_yticks(list(np.arange(0, longest_y, longest_y/5))+[longest_y])
            ax[d].legend()
        plt.show()



if __name__=="__main__":
    main()


