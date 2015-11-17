objects = atominfo.obj,cputime.obj,evalfw.obj,findpnnb.obj,-
frame.obj,io.obj,lattice.obj,lib.obj,matrix.obj,posnmx.obj,-
psearch.obj,ranmar.obj,recycle.obj,scatfact.obj,sgextend.obj,-
setaos.obj,signal.obj,sitefcal.obj,tetrahed.obj,trialoop.obj,-
xtal.obj,main.obj,fourier.obj

cflags = $(cflags)/optimize=inline/prefix_library_entries=all_entries
!linkflags = $(linkflags)/nosysshr

!cflags = $(cflags)/debug /prefix_library_entries=all_entries
!linkflags = $(linkflags)/debug

focus.exe : $(objects)
        $(link)$(linkflags)/notrace $(objects),libsginfo/lib,libfftpack/lib

fourier.obj : fourier.c
        $(cc) $(cflags) /define=("FFTfftpack") fourier.c

atominfo.exe : atominfo.c
        $(cc) $(cflags) atominfo.c /define=("STAND_ALONE")
        $(link)$(linkflags)/notrace atominfo.obj

cputime.exe : cputime.c
        $(cc) $(cflags) cputime.c /define=("STAND_ALONE")
        $(link)$(linkflags)/notrace cputime.obj

ranmar.exe : ranmar.c
        $(cc) $(cflags) ranmar.c /define=("STAND_ALONE")
        $(link)$(linkflags)/notrace ranmar.obj

.c.obj :
        $(cc) $(cflags) $(mms$source)
