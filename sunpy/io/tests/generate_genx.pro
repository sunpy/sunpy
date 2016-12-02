; docformat = 'rst'
;+
; NAME:
;  GENERATE_GENX
;
; PURPOSE:
;  Provide a test file for being tested by SunPy
;  It needs solarsoft to be run and to create the file
;
; EXAMPLE:
;   idl> .r generate_genx.pro
;   idl> generate_genx.pro
;   idl> restgen,file='testing.genx', STRUCT=TESTING, HEAD=HEAD, TEXT=TEXT
;   idl> help, TESTING,/str
;   ** Structure MS_216285419006, 9 tags, length=3720, data length=3658:
;      MYTEXT          STRING    'Waiter to Mr Creosote: 'Finally, monsieur – a wafer-thin mint.''
;      MYTEXT_ARRAY    STRING    Array[3]
;      MYTEXT_ARRAY_DIMENSION
;                      STRING    Array[2, 3]
;      MYNUMBER        INT              1
;      MYNUMBER_ARRAY  INT       Array[3]
;      MYNUMBER_ARRAY_DIMENSION
;                      INT       Array[2, 3, 4, 5]
;      MYUINT          UINT             2
;      MYSTRUCTURE     STRUCT    -> MS_216285419005 Array[1]
;      MYSTRUCTURE_ARRAY
;                      STRUCT    -> MS_216285419005 Array[6]
;   idl> print, TESTING.MYNUMBER_ARRAY_DIMENSION[0,0,0,0]
;   0
;   idl> print, TESTING.MYNUMBER_ARRAY_DIMENSION[0,0,0,4]
;   96
;   idl> print, TESTING.MYNUMBER_ARRAY_DIMENSION[0,2,2,0]
;   16
;   idl> print, TESTING.MYNUMBER_ARRAY_DIMENSION[0,2,3,0]
;   22
;   idl> print, TESTING.MYNUMBER_ARRAY_DIMENSION[0,2,3,4]
;   118
;   idl> print, size(TESTING.MYSTRUCTURE_ARRAY)
;           1           6           8           6
;   idl> print, TESTING.MYSTRUCTURE.MYFARRAY
;      0.00000      1.00000      2.00000
;   idl> print, TESTING.MYSTRUCTURE.MYFARRAYD[0,*]
;      0.00000
;      2.00000
;      4.00000
;   idl> print, TESTING.MYSTRUCTURE.MYDARRAYD[2,1]
;       5.0000000
;   idl> print, TESTING.MYSTRUCTURE.NESTEDSTRUCT.MYLARRAYD[1,0,3]
;       19
;   idl> print, TESTING.MYSTRUCTURE.NESTEDSTRUCT.MYLARRAYD[0,*,2]
;       12
;       14
;       16
;   idl> print, TESTING.MYSTRUCTURE.MYCARRAY[1]
;       (      1.00000,     -9.00000)
;   idl> print, TESTING.MYSTRUCTURE.MYDCARRAY[2]
;       (     12.00000,      1.00000)
;   idl> print, TESTING.MYSTRUCTURE.NESTEDSTRUCT.MYL64NUMBER
;       9223372036854775807
;   idl> print, TESTING.MYSTRUCTURE.NESTEDSTRUCT.MYUL64NUMBER
;       18446744073709551615
;   idl> print, size(TESTING.MYSTRUCTURE.NESTEDSTRUCT.MYLARRAYD[0,*,2])
;       2           1           3           3           3
;   idl> help, HEAD, TEXT,/str
;   ** Structure MS_216109884001, 4 tags, length=72, data length=72:
;     VERSION         LONG                 2
;     XDR             LONG                 1
;     CREATION        STRING    'Sat Oct 29 08:15:08 2016'
;     IDL_VERSION     STRUCT    -> IDLV3_VERSION Array[1]
;     TEXT            STRING    = 'Sample crazy file to test the sunpy genx reader'
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
  myL64number = long64(9223372036854775807)

  ;; ulong64
  myUL64number = ulong64(18446744073709551615)

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

  ;; complex
  myComplex = complex(9, 10)
  myComplex_array = complex(findgen(3), findgen(3)-10)

  ;; dcomplex
  myDComplex = dcomplex(15, 20)
  myDComplex_array = dcomplex(findgen(3)+10, findgen(3)-1)

  ;; Structure
 myStructure = {dummy, name: 'floats and doubles', $
                 myfnumber: myfloatnumber, myfarray: myfloatnumber_array, myfarrayd: myfloatnumber_array_dimension, $
                 myDnumber: myDnumber, myDarray: myDnumber_array, myDarrayd: myDnumber_array_dimension, $
                 nestedstruct: {subdummy, spam: 'brian, longs and ulongs', $
                                myLnumber: myLnumber, myLarray: myLnumber_array, myLarrayd: myLnumber_array_dimension, $
                                myULnumber: myUlong, myL64number: myL64number, myUL64number: myUL64number}, $
                 mytext: 'adding also complex', $
                 randomNumbers: [87, 105, 116, 104, 32, 108, 111, 118, 101, 32, 102, 114, 111, 109, 32, 83, 97, 110, $
                                 32, 70, 114, 97, 110, 99, 105, 115, 99, 111, 44, 32, 50, 48, 49, 54], $
                 myCnumber: myComplex, myCarray: myComplex_array, $
                 myDCnumber: myDComplex, myDCarray: myDComplex_array}


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
