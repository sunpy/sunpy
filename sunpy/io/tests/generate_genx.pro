; docformat = 'rst'
;+
; NAME:
;  GENERATE_GENX
; PURPOSE:
;  Provide a test file for being tested by SunPy
;  It needs solarsoft to be run and to create the file
;-
pro generate_genx

  ;; variables to create:
  ;; Single text element
  mytext = "Waiter to Mr Creosote: 'Finally, monsieur – a wafer-thin mint.'"
  ;; Text array (one dimension, 3 elements)
  mytext_array = replicate(mytext, 3)
  ;; multi-dimension string
  mytext_array_dimension = replicate(mytext, 2, 3)

  ;; single numbers (unsigned)
  mynumber = 1
  ;; array of one dimension
  mynumber_array = indgen(3)
  ;; multiple dimension:
  mynumber_array_dimension = indgen(2, 3, 4, 5)

  ;; uint
  myUint = uint(2)

  ;; long numbers
  myLnumber = 3L
  ;; array of one dimension
  myLnumber_array = lindgen(3)
  ;; multiple dimension:
  myLnumber_array_dimension = lindgen(2, 3, 4)

  ;; ulong
  myUlong = ulong(4)

  ;; long 64
  myL64number = long64(5)

  ;; ulong64
  myUL64number = ulong64(6)

  ;; float numbers
  myfloatnumber = 7.
  ;; array of one dimension
  myfloatnumber_array = findgen(3)
  ;; multiple dimension:
  myfloatnumber_array_dimension = findgen(2, 3)

  ;; double numbers
  myDnumber = 8D
  ;; array of one dimension
  myDnumber_array = dindgen(3)
  ;; multiple dimension:
  myDnumber_array_dimension = dindgen(3, 2)

  ;; TODO: complex
  ;; myComplex = complex(9, 10)
  ;; myComplex_array = complex(findgen(3), findgen(3)-10)

  ;; TODO: dcomplex
  ;; myDComplex = dcomplex(9, 10)
  ;; myDComplex_array = dcomplex(findgen(3), findgen(3)-10)

  ;; Structure
 myStructure = {dummy, name: 'floats and doubles', $
                 myfnumber: myfloatnumber, myfarray: myfloatnumber_array, myfarrayd: myfloatnumber_array_dimension, $
                 myDnumber: myDnumber, myDarray: myDnumber_array, myDarrayd: myDnumber_array_dimension, $
                 nestedstruct: {subdummy, spam: 'brian, longs and ulongs', $
                                myLnumber: myLnumber, myLarray: myLnumber_array, myLarrayd: myLnumber_array_dimension, $
                                myULnumber: myUlong}, $ ;, myL64number: myL64number, myUL64number: myUL64number}, $
                 mytext: 'adding also complex', $
                 randomNumbers: [87, 105, 116, 104, 32, 108, 111, 118, 101, 32, 102, 114, 111, 109, 32, 83, 97, 110, $
                                 32, 70, 114, 97, 110, 99, 105, 115, 99, 111, 44, 32, 50, 48, 49, 54] }

                 ;; myCnumber: myComplex, myCarray: myComplex_array, $
                 ;; myDCnumber: myDComplex, myDCarray: myDComplex_array, $


  ;; array of structures
  myStructure_array = replicate(myStructure, 2, 3)

  ;; Let's save them into a genx

  savegen, mytext, mytext_array, mytext_array_dimension, $
           mynumber, mynumber_array, mynumber_array_dimension, $
           myUint, $
           myStructure, $
           myStructure_array, $
           names=['mytext', 'mytext_array', 'mytext_array_dimension', $
                  'mynumber', 'mynumber_array', 'mynumber_array_dimension', $
                  'myUint', 'myStructure', 'myStructure_array'], $
           file='testing', /xdr, $
           text='Sample crazy file to test the sunpy genx reader'
end
