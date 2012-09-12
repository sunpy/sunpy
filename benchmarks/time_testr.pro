PRO time_testr, NREP = nrep
;
; Purpose: To test IDL, the Richard Schwartz way.
;   I'll be testing matrix multiplication,
;   sorting, direct addressing, shifting, and indirect addressing.
; 
; Written: Richard Schwartz 

default, nrep, 10
t= fltarr(20)
ntst = 1

;Test 1
a = findgen(500,500L)
b = a
t1 = systime(1)
for i=0, nrep-1 do c = a # b
t[0] = systime(1)-t1
Print, '1', T[0], '  Matrix Multiplication Large Arrays (500,500) ' + num2str(nrep) + ' times'

;Test 2
a = findgen(50,50L)
b = a
t1 = systime(1)
for i=0,nrep*1000-1 do c = a # b
t[1] = systime(1)-t1
Print, '2', T[1], '  Matrix Multiplication Small Array (50,50) ' +num2str(1000*nrep) + ' times'

;Test 3
a = randomu(s,1e6)
t1 = systime(1)
for i=0,nrep-1 do c = sort(a)
t[2] = systime(1)-t1
Print, '3', T[2], '  Sorting 1 million elements ' + num2str(nrep) + ' times.'

;Test 4
a = randomu(s,1e6)
c = fltarr(1e6)
t1 = systime(1)
for i=0,nrep*100-1 do c[0] = a
t[3] = systime(1)-t1
Print, '4', T[3], '  Moving 1 million elements ' + num2str(nrep*100) + ' times.'

;Test 5
a = (randomu(s,1e6)*1e6)
b= long(a)
c = fltarr(1e6)
t1 = systime(1)
for i=0,nrep*10-1 do c[0] = a[b]
t[4] = systime(1)-t1
Print, '5', T[4], '  indirect addressing 1 million elements ' + num2str(nrep*10) + ' times.'

;Test 6
a = ulong(randomu(s,1e6)*1e6)
c = ulonarr(1e6)
t1 = systime(1)
for i=0,nrep*100-1 do c[0] = ishft(a,12)
t[5] = systime(1)-t1
Print, '6', T[5], '  shifting 1 million elements ' + num2str(nrep*100) + ' times.'

;Test 7
a = (randomu(s,1e6)*1e6)

c = fltarr(1e6)
t1 = systime(1)
for i=0,nrep*10-1 do c[0] = cos(a)
t[6] = systime(1)-t1
Print, '7', T[6], '  cosine 1 million elements ' + num2str(nrep*10) + ' times.'

;Test 8
a = (randomu(s,1e6)*1e6)

c = fltarr(1e6)
t1 = systime(1)
for i=0,nrep*10-1 do c[0] = alog(a)
t[7] = systime(1)-t1
Print, '7', T[7], '  alog 1 million elements ' + num2str(nrep*10) + ' times.'
;Result = ISHFT(P1, P2)

;Test 9
a = bindgen(1e6)
b = bytarr(1e6)
openw, lu, /get,'timetest.dat'

t1 = systime(1)
for i=0,nrep*100-1 do begin
	writeu, lu, a
	point_lun,lu, 0
	readu, lu, b
	endfor
t[8] = systime(1)-t1
free_lun, lu
Print, '8', T[8], '  writing and reading bytarr(1e6) ' + num2str(nrep*100) + ' times.'
print, total(t), '  Richard Schwartz Time Test Result'

end

