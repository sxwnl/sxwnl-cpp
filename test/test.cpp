/*
寿星天文历v5.05js版 C++翻译工程
API测试
2018-8-28
*/

#include <time.h>
#include <iostream>
#include <array>
#include <map>

#include "../tool.h"
#include "../lunar.h"
#include "../eph_show.h"

#define MAP_H 18
#define MAP_W  7
std::array<std::array<std::string,MAP_W>,MAP_H> strmap;

Date get_time(void)
{
	struct tm *bjs;
	time_t time0;
	time0 = time(NULL);
	bjs = localtime(&time0);
	Date dat=
	{
		bjs->tm_year + 1900,
		bjs->tm_mon+1,
		bjs->tm_mday,
		bjs->tm_hour,
		bjs->tm_min,
		(double)bjs->tm_sec
	};
	return dat;
}
void initmap(int y,int m)
{
	OB_LUN lun=yueLiCalc(y,m);
	
	std::map<std::string,std::string> str_yx={{"朔","\e[33m●"},{"上弦","\e[33m∪"},{"望","\e[33m∩"},{"下弦","\e[38;5;245m●"}};
	for (int i=0;i<18;i++)
	{
    	for (int j=0;j<7;j++)
    	{
    		if (i%3==2)
    		strmap[i][j]="\e[32m————————";
        	else 
        	strmap[i][j]="        ";
    	}
	}
	for (int i = 0;i < lun.dn;i++)
	{
		int y=lun.day[i].weeki;
		int x=lun.day[i].week;
		int j=i<9?6:4;
		strmap[y*3][x] = x==0||x==6?"\e[31m":"\e[37m";
		strmap[y*3][x] += str_rmc0[i];
		if (lun.day[i].jqmc.length())
			strmap[y*3][x]+="\e[32m◆",j--;
		if (lun.day[i].yxmc.length())
			strmap[y*3][x]+=str_yx[lun.day[i].yxmc],j--;
		for (int k = 0;k < j;k++)
			strmap[y*3][x]+=" \e[0m\e[1m";
		if (!j) 
			strmap[y*3][x]+="\e[0m\e[1m";
		if (lun.day[i].Ldi==0)
			strmap[y*3+1][x]="\e[38;5;123m",strmap[y*3+1][x]=strmap[y*3+1][x]+lun.day[i].Lleap+strmap[y*3+1][x]+lun.day[i].Lmc2+strmap[y*3+1][x]+"月"+strmap[y*3+1][x]+str_dx[lun.day[i].Ldn-29]+strmap[y*3+1][x]+(lun.day[i].Lleap.length()?"\e[0m\e[1m":"  \e[0m\e[1m");
		else
			strmap[y*3+1][x]="",strmap[y*3+1][x]=strmap[y*3+1][x]+lun.day[i].Ldc+strmap[y*3+1][x]+"    \e[0m\e[1m";
		if (lun.day[i].Ljq.length())
			strmap[y*3+1][x]="\e[34m",strmap[y*3+1][x]=strmap[y*3+1][x]+lun.day[i].Ljq+strmap[y*3+1][x]+="    \e[0m\e[1m";
			
	}
//	sprintf(strmap_head[0],"公元 %04d年%02d月%02d日 周%s %s座\n",lun.day[1].XiZ);
}
void drawmap()
{
	for (int k=0;k<MAP_W;k++)
		std::cout<<"\e[3"<<((k&&k<6)?3:1)<<";1m周"<<str_xqmc[k]<<"    ";
	std::cout<<std::endl;
	for (int k=0;k<MAP_W*4;k++)
		std::cout<<"\e[33m——";
	std::cout<<std::endl;
	for (std::array<std::string,MAP_W> i:strmap)
	{
		for (std::string j:i)
			std::cout<<j;
		std::cout<<std::endl;
	}
	std::cout<<"\e[0m";
}
int main()
{
	init_ob();
	initmap(2017,7);
	drawmap();
	std::cout<<DD2str(get_time())<<std::endl;
	/*	
	
	OB_LUN lun=yueLiCalc(2019,7);
	for (int i = 0;i < 30;i++)
	{
		std::cout<<std::endl<<lun.day[i].d<<"天:"<<lun.day[i].Lday2;
		if (lun.day[i].A.length())
		std::cout<<"A"<<lun.day[i].A;
		if (lun.day[i].B.length())
		std::cout<<"B"<<lun.day[i].B;
		if (lun.day[i].C.length())
		std::cout<<"C"<<lun.day[i].C;
	
	}*/
	rysCalc();
	std::cout<<rs_search(2008,8,200,1)<<std::endl;
	rs2_calc(5,0);
	rs2_jxb();
	return 0;
}