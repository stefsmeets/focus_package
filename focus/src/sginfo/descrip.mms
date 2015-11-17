cflags = $(cflags)/optimize=inline/prefix_library_entries=all_entries

sginfo.exe : sginfo.obj sgvers.obj sgclib.obj sgio.obj sgfind.obj sghkl.obj sgss.obj sgharker.obj sgutil.obj sgssvfy.obj xtal32ss.obj
        $(link)$(linkflags)/notrace sginfo.obj, sgvers.obj, sgclib.obj, sgio.obj, sgfind.obj, sghkl.obj, sgss.obj, sgharker.obj, sgutil.obj, sgssvfy.obj, xtal32ss.obj

sgtest.exe : sgtest.obj sgvers.obj sgclib.obj sgio.obj sgfind.obj sgutil.obj
        $(link)$(linkflags)/notrace sgtest.obj, sgvers.obj, sgclib.obj, sgio.obj, sgfind.obj, sgutil.obj

sgquick.exe : sgquick.obj sgvers.obj sgclib.obj sgio.obj sgutil.obj
        $(link)$(linkflags)/notrace sgquick.obj, sgvers.obj, sgclib.obj, sgio.obj, sgutil.obj

libsginfo : sgvers.obj sgclib.obj sgio.obj sgfind.obj sghkl.obj sgss.obj sgharker.obj sgutil.obj
        library /create libsginfo sgvers, sgclib, sgio, sgfind, sghkl, sgss, sgharker, sgutil

clean :
        delete *.obj;*
        delete *.olb;*
        delete *.exe;*

.c.obj :
        $(cc) $(cflags) $(mms$source)
