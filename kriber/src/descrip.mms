.ignore

objects = kriber.obj, p2clib.obj

cflags = $(cflags)/optimize=inline/prefix_library_entries=all_entries

kriber.exe : $(objects)
        $(link)$(linkflags)/notrace $(objects)

.c.obj :
        $(cc) $(cflags) $(mms$source)
