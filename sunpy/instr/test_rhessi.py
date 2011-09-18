import sunpy

from sunpy.instr import rhessi
result = rhessi.backprojection(sunpy.RHESSI_EVENT_LIST)
print(result)

result = rhessi.backprojection(sunpy.RHESSI_EVENT_LIST, detector=7)
print(result)
