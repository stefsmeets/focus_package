 CC = cc
 CCFLAGS = $(CFLAGS)

sginfo-query.x: sginfo-query.o cgitbx.o
	$(CC) $(CCFLAGS) sginfo-query.o cgitbx.o -o sginfo-query.x

encodeurl: encodeurl.o cgitbx.o
	$(CC) $(CCFLAGS) encodeurl.o cgitbx.o -o encodeurl

.c.o: 
	$(CC) $(CCFLAGS) -c $<

clean:
	rm -f sginfo-query.x sginfo-query.o cgitbx.o
	rm -f encodeurl encodeurl.o
