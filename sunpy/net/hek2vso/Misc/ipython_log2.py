# IPython log file

get_ipython().magic(u'run mod-load.ipy')
get_ipython().magic(u'store ')
get_ipython().magic(u'logstate ')
import json
h2v = hek2vso.HEK2VSOTool()
temp = h2v.return_hek_from_file()
get_ipython().magic(u'who ')
get_ipython().magic(u'whos ')
print temp
print temp[0]
print json.dumps(temp[0], sort_keys=True, indent=4, separators=(',', ': "))
test = json.loads(temp[0])
test = json.loads(temp)
get_ipython().magic(u'whos ')
print json.dumps(temp[0], sort_keys=True, indent=4, separators=(',', ': "))
print json.dumps(temp[0], sort_keys=True, indent=4, separators=(',', ': '))
get_ipython().magic(u'whos ')
print temp.keys
print temp
print temp[0]
get_ipython().magic(u'whos ')
print temp[0].keys()
print json.dumps(temp[0].keys(), sort_keys=True, indent=4)
print json.dumps(temp[0].keys().sort(), indent=4)
x = temp[0].keys()
x = x.sort()
print x
x = temp[0].keys()
type(x)
#[Out]# list
x = x.sort()
type(x)
#[Out]# NoneType
x = temp[0].keys()
print x
print sorted(x)
x = sorted(x)
print json.dumps(x, indent=4, separators=(','))
print json.dumps(x, indent=4)
def quick_clean(x):
    for i in x:
        if i == "'" or i == ",":
            i = ''
    return x
print quick_clean("hopef,ully thi's will', work")
def quick_clean(x):
    for i in x:
        if i == "'" or i == ",":
            print i
    return x
print quick_clean("hopef,ully thi's will', work")
def quick_clean(x):
    for i in x:
        if i == "'" or i == ",":
            print i
    return x
def quick_clean(x):
    clean = ''
    for i in x:
    	if i is not "'" or i is not ",":
         clean += i
    return clean
print quick_clean("this, i's a te',t")
a = "'"
a == "'"
#[Out]# True
a == '
def quick_clean(x):
    clean = ''
    for i in x:
    	if i is not "'" and i is not ",":
         clean += i
    return clean
print quick_clean("this, i's a te',t")
get_ipython().system(u'ls')
with open(HEK\ attrs\ list, "r") as f:
    temp = f.read()
    temp = quick_clean(temp)
    print temp
    
with open(HEK\ attrs\ list, 'r') as lf:
    temp = lf.read()
    temp = quick_clean(temp)
    print temp
    
with open('HEK\ attrs\ list', 'r') as lf:
    temp = lf.read()
    temp = quick_clean(temp)
    print temp
    
with open('HEK attrs list', 'r') as lf:
    temp = lf.read()
    temp = quick_clean(temp)
    print temp
    
def quick_clean(x):
    clean = ''
    for i in x:
    	if i is not '"' and i is not ",":
         clean += i
    return clean
with open('HEK attrs list', 'r') as lf:
    temp = lf.read()
    temp = quick_clean(temp)
    print temp
    
import pygame
get_ipython().magic(u'logstate ')
get_ipython().magic(u'logstop ')
