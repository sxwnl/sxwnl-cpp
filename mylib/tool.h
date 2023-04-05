#ifndef TOOL_H
#define TOOL_H

#include <cstdlib>
#include "../mylib/mystl/string.h"
#include "../mylib/math_patch.h"

#define A2R(x) ((x) / 180 * M_PI)

inline int int2(double v) {return (int)floor(v);};


struct Date
{//儒略历结构，包含: 年 月 日 时 分 秒 (浮点型)
	int Y, M, D, h, m;
	double s;
};

//实现数字转字符
// mystl::string to_str(int in, int n, int len);
// mystl::string to_str(double in, int n, int len);

//mystl::vector<mystl::string> split(const mystl::string& src, const mystl::string& separator);
void string_replace( mystl::string &strBig, const mystl::string &strsrc, const mystl::string &strdst);
mystl::string timeStr(double jd);
mystl::string rad2strE(double d, bool tim, int ext);
mystl::string rad2str(double d, bool tim);
mystl::string rad2str2(double d);
mystl::string m2fm(double v, int fx, int fs);
double toJD(Date date);
Date setFromJD(double jd);
mystl::string DD2str(Date r);
mystl::string JD2str(double jd);
mystl::string fill_str(mystl::string s, int n, mystl::string c);

mystl::string to_str(long in);
mystl::string to_str(int in);
mystl::string to_str(int in, uint8_t n);
mystl::string to_str(double in);
mystl::string to_str(double in, uint8_t n);
mystl::string to_str(double in, uint8_t n, uint8_t n2);
int astoi(mystl::string as);

#endif