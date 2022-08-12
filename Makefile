CFLAGS	:= -c
BUILD_PATH	:= build
VPATH	:= eph : lunar : draw : test
OBJ	:= test1.o eph0.o eph.o eph_msc.o eph_rsgs.o eph_rspl.o eph_yspl.o \
eph_szj.o eph_show.o lunar.o lunar_ob.o lunar_ssq.o tool.o lat_lon_data.o
OBJ0	:= test.o eph0.o eph.o eph_msc.o eph_rsgs.o eph_rspl.o eph_yspl.o \
eph_szj.o eph_show.o lunar.o lunar_ob.o lunar_ssq.o tool.o lat_lon_data.o
LIBS	:= -lm

test1: $(OBJ)
	$(CXX) $(LIBS) $^ -o $@
test0: $(OBJ0)
	$(CXX) $(LIBS) $^ -o $@
.PHONY: test1 test0 clean cleanw
clean:
	rm -rf *.o
cleanw:
	del *.o