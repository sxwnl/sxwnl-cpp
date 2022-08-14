#ifndef TOOL_H
#define TOOL_H

#define _USE_MATH_DEFINES

#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <cstdio>
#define A2R(x) ((x) / 180 * M_PI)

inline int int2(double v) {return (int)floor(v);};

struct Date
{//儒略历结构，包含: 年 月 日 时 分 秒 (浮点型)
	int Y, M, D, h, m;
	double s;
};

//实现数字转字符
template <class T>
std::string to_str(T in, int n = 0, int len = 0)
{
	std::ostringstream str;
	str << (in >= 0 && n != -1 ? " " : "") << std::fixed << std::setprecision(n) << in;
	std::string s = str.str();
	int slen = s.length()-(n?2:1);
	for (int i = 0; i < (len + n - slen) && len; i++)
	{
		s = " " + s;
	}
	return s;
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
std::string fill_str(std::string s, int n, std::string c);
#endif