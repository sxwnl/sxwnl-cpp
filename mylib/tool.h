#ifndef TOOL_H
#define TOOL_H

#include <cstdlib>
#include "../mylib/mystl/my_string.h"
#include "../mylib/math_patch.h"

#define A2R(x) ((x) / 180 * M_PI)

inline int int2(double v) {return (int)floor(v);};


struct Date
{//儒略历结构，包含: 年 月 日 时 分 秒 (浮点型)
	int Y, M, D, h, m;
	double s;
};

//实现数字转字符
// mystl::my_string to_str(int in, int n, int len);
// mystl::my_string to_str(double in, int n, int len);

//mystl::vector<mystl::my_string> split(const mystl::my_string& src, const mystl::my_string& separator);
void string_replace( mystl::my_string &strBig, const mystl::my_string &strsrc, const mystl::my_string &strdst);
mystl::my_string timeStr(double jd);
mystl::my_string rad2strE(double d, bool tim, int ext);
mystl::my_string rad2str(double d, bool tim);
mystl::my_string rad2str2(double d);
mystl::my_string m2fm(double v, int fx, int fs);
double toJD(Date date);
Date setFromJD(double jd);
mystl::my_string DD2str(Date r);
mystl::my_string JD2str(double jd);
mystl::my_string fill_str(mystl::my_string s, int n, mystl::my_string c);

mystl::my_string to_str(long in);
mystl::my_string to_str(int in);
mystl::my_string to_str(int in, uint8_t n);
mystl::my_string to_str(double in);
mystl::my_string to_str(double in, uint8_t n);
mystl::my_string to_str(double in, uint8_t n, uint8_t n2);
int astoi(mystl::my_string as);

#endif