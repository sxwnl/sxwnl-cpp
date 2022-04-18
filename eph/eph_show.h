#ifndef EPH_SHOW_H
#define EPH_SHOW_H

#include <string>
#include "../tool.h"

void rysCalc(Date d, bool is_utc, bool nasa_r);
std::string rs_search(int Y,int M,int n,bool fs);
void rs2_calc(int fs,double jd0);
void rs2_jxb();
void shengjiang(int y, int m, int d);
void shengjiang2(int y);
void shengjiang3(int y);
#endif