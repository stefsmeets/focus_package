$! KRIBER COMMAND FILE
$!
$ type sys$input


   KRIBER (generates a DLS or LOADAT input file)


$ inquire/p  input     "  DLS76 input file name (default = DLS76.INP) "     
$ if input .eqs. "" then input = "dls76.inp"
$ if f$locate(".",input) .eqs. f$length(input) then input = ''input'+".inp"
$ kriberdir = "sys$products:[xray.kriber]"
$ dosdat = "''kriberdir'DISTDAT.DAT"
$ inquire    check    " Do you need your own DISTDAT.DAT ?"
$ if check .ne. "y" then goto torf
$ fima=f$pars("sys$login",,,"device")+f$pars("sys$login",,,"directory")
$ ouff=f$extract(0,f$locate("]",fima),fima)
$ gont:
$ filese=f$search(''ouff'+"...]DISTDAT.DAT;")
$ if filese .eqs. "" then goto torf
$ inquire    check    " Do you need ''filese'?"
$ if check .ne. "y" then goto gont
$ dosdat=filese
$ torf:
$ assign  'input'  dlsinp
$ assign  'dosdat' distdat
$ assign  'kriberdir'SYMDAT.DAT symdat
$ assign  'kriberdir'COSEQ.DAT coseq 
$ assign  loadat.inp  loadainp
$ assign  structure.dat  strudat
$ assign  kriber.list  list
$ assign  structure_new.dat  strudatn
$ define/user_mode sys$input sys$command
$ run 'kriberdir'kriber
$ deassign/all
$ exit       
