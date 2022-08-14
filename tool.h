#ifndef TOOL_H
#define TOOL_H

#define _USE_MATH_DEFINES

#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#define A2R(x) ((x) / 180 * M_PI)

inline int int2(double v) {return (int)floor(v);};

struct Date
{//儒略历结构，包含: 年 月 日 时 分 秒 (浮点型)
	int Y, M, D, h, m;
	double s;
};

//实现数字转字符
template <class T>
std::string to_str(T in, int n = 0)
{
	std::ostringstream str;
	str << (n==-1?(in>0?" ":""):"")<<std::fixed << std::setprecision(n) << in;
	return str.str();
}

std:: vector<std::string> split(const std::string& src, const std::string& separator);
void string_replace( std::string &strBig, const std::string &strsrc, const std::string &strdst);
std::string timeStr(double jd);
std::string rad2strE(double d, bool tim, int ext);
std::string rad2str(double d, bool tim);
std::string rad2str2(double d);
std::string m2fm(double v, int fx, int fs);
double toJD(Date date);
Date setFromJD(double jd);
std::string DD2str(Date r);
std::string JD2str(double jd);
#endif