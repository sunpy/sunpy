# IPython log file

get_ipython().magic(u'logstate ')
import os
import json
f = open('test.json', 'r')
r = f.read()
print r
print json.dumps(json.loads(r), indent = 4, separators = (',
print json.dumps(json.loads(r), indent = 4, separators = (',
print json.dumps(json.loads(r), indent = 4, separators = (',', ': '))
json.dumps(json.loads(r), indent = 4, separators = (',', ': '))
#[Out]# '{\n    "a": 15,\n    "c": 45,\n    "b": 23\n}'
f.close()
get_ipython().magic(u'logstate ')
get_ipython().magic(u'who ')
get_ipython().magic(u'whos ')
f.close()
with open('test.json', 'r') as f:
    r = f.read()
    
print r
import sunpy.hek
import sunpy
import sunpyt.net as sn
import sunpy.net as sn
tstart = '2011/08/09 07:23:56'
tend = '2011/08/09 12:40:29'
event_type = 'FL'
from sn import hek
from sunpy.net import hek
client = hek.HEKClient()
result = client.query(hek.attrs.Time(tstart, tend), hek.attrs.EventType(event_type))
print result
with open('hek_query-result', 'w') as f:
    f.write(result)
    
get_ipython().system(u'ls')
get_ipython().system(u'less hek_query-result')
get_ipython().magic(u'whos ')
with open('hek_query-result', 'w') as f:
    f.write(json.dump(result))
    
get_ipython().magic(u'pinfo json.dump()')
get_ipython().magic(u'pinfo json.dump')
with open('hek_query-result', 'w') as f:
    f.write(json.dumps(result))
    
with open('hek_query-result', 'r') as f:
    temp = f.read()
    print temp
    
print temp
get_ipython().magic(u'who ')
print json.loads(json.dumps(temp))
print json.dumps(json.loads(temp), sort_keys = True, indent = 4, separators = (',', ': '))
json_buffer = json.dumps(json.loads(temp), sort_keys = True, indent = 4, separators = (',', ': '))
get_ipython().magic(u'whos ')
json_buffer = json.loads(temp)
get_ipython().magic(u'whos ')
_oh
#[Out]# {11: '{\n    "a": 15,\n    "c": 45,\n    "b": 23\n}'}
_ih
#[Out]# ['',
#[Out]#  u"get_ipython().magic(u'logstate ')",
#[Out]#  u"get_ipython().magic(u'logstart -o')",
#[Out]#  u'import os',
#[Out]#  u'import json',
#[Out]#  u"f = open('test.json', 'r')",
#[Out]#  u'r = f.read()',
#[Out]#  u'print r',
#[Out]#  u"print json.dumps(json.loads(r), indent = 4, separators = (',",
#[Out]#  u"print json.dumps(json.loads(r), indent = 4, separators = (',",
#[Out]#  u"print json.dumps(json.loads(r), indent = 4, separators = (',', ': '))",
#[Out]#  u"json.dumps(json.loads(r), indent = 4, separators = (',', ': '))",
#[Out]#  u'f.close()',
#[Out]#  u"get_ipython().magic(u'logstate ')",
#[Out]#  u"get_ipython().magic(u'who ')",
#[Out]#  u"get_ipython().magic(u'whos ')",
#[Out]#  u'f.close()',
#[Out]#  u"with open('test.json', 'r') as f:\n    r = f.read()\n    ",
#[Out]#  u'print r',
#[Out]#  u'import sunpy.hek',
#[Out]#  u'import sunpy',
#[Out]#  u'import sunpyt.net as sn',
#[Out]#  u'import sunpy.net as sn',
#[Out]#  u"tstart = '2011/08/09 07:23:56'",
#[Out]#  u"tend = '2011/08/09 12:40:29'",
#[Out]#  u"event_type = 'FL'",
#[Out]#  u'from sn import hek',
#[Out]#  u'from sunpy.net import hek',
#[Out]#  u'client = hek.HEKClient()',
#[Out]#  u'result = client.query(hek.attrs.Time(tstart, tend), hek.attrs.EventType(event_type))',
#[Out]#  u'print result',
#[Out]#  u"with open('hek_query-result', 'w') as f:\n    f.write(result)\n    ",
#[Out]#  u"get_ipython().system(u'ls')",
#[Out]#  u"get_ipython().system(u'less hek_query-result')",
#[Out]#  u"get_ipython().magic(u'whos ')",
#[Out]#  u"with open('hek_query-result', 'w') as f:\n    f.write(json.dump(result))\n    ",
#[Out]#  u"get_ipython().magic(u'pinfo json.dump()')",
#[Out]#  u"get_ipython().magic(u'pinfo json.dump')",
#[Out]#  u"with open('hek_query-result', 'w') as f:\n    f.write(json.dumps(result))\n    ",
#[Out]#  u"with open('hek_query-result', 'r') as f:\n    temp = f.read()\n    print temp\n    ",
#[Out]#  u'print temp',
#[Out]#  u"get_ipython().magic(u'who ')",
#[Out]#  u'print json.loads(json.dumps(temp))',
#[Out]#  u"print json.dumps(json.loads(temp), sort_keys = True, indent = 4, separators = (',', ': '))",
#[Out]#  u"json_buffer = json.dumps(json.loads(temp), sort_keys = True, indent = 4, separators = (',', ': '))",
#[Out]#  u"get_ipython().magic(u'whos ')",
#[Out]#  u'json_buffer = json.loads(temp)',
#[Out]#  u"get_ipython().magic(u'whos ')",
#[Out]#  u'_oh',
#[Out]#  u'_ih']
get_ipython().magic(u'pinfo hek.query')
