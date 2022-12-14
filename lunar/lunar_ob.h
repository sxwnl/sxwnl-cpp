#ifndef LUNAR_OB_H
#define LUNAR_OB_H

#include <string>
#include <vector>

struct OB_DAY
{
	/* 日的公历信息 */
	int d0;		// 2000.0起算儒略日,北京时12:00
	int di;		// 所在公历月内日序数
	int y;
	int m;
	int d; 		// 日名称(公历)
	int dn;
	int week0;
	int week;  	// 星期
	int weeki; 	// 在本月中的周序号
	int weekN; 	// 本月的总周数
	/* 日的农历信息 */
	int Ldi;	   // 日距农历月首偏移
	std::string Ldc;   // 日名称(农历),即'初一,初二等'
	
	int cur_dz; //冬至的天数
	int cur_xz; //夏至的天数
	int cur_lq; //立秋的天数
	int cur_mz; //芒种的天数
	int cur_xs; //小暑的天数
	std::string Lmc;
	std::string Lmc2;  // 月名称
	int Ldn;   		// 月大小
	std::string Lleap; 	// 闰状况 
	/* 日的农历纪年、月、日、时及星座 */
	int Lyear;		 	// 农历纪年(10进制,1984年起算,分界点可以是立春也可以是春节,在程序中选择一个)
	int Lyear0;
	std::string Lyear2;	// 干支纪年
	std::string Lyear3;	// 干支纪年(春节)
	int Lyear4;			// 干支纪年(黄帝纪元)
	int Lmonth;			// 纪月处理,1998年12月7日(大雪)开始连续进行节气计数,0为甲子
	std::string Lmonth2;   // 干支纪月
	std::string Lday2; 	// 纪日
	std::string Ltime2;	// 纪时
	std::string Ljq;	   // 节气
	std::string XiZ;   	// 星座
	/* 日的回历信息 */
	int Hyear;	 	// 年(回历)
	int Hmonth;		// 月(回历)
	int Hday;	  	// 日(回历)
	/* 日的其它信息 */
	std::string yxmc;	// 月相名称
	double yxjd;   	  // 月相时刻(儒略日)
	std::string yxsj;	// 月相时间串
	std::string jqmc;	// 节气名称
	double jqjd;	     // 节气时刻(儒略日)
	std::string jqsj;	// 节气时间串
	
	bool Fjia;
	std::string A;
	std::string B;
	std::string C;
};

struct MLBZ
{
	std::string bz_jn;
	std::string bz_jy;
	std::string bz_jr;
	std::string bz_js;
	std::string bz_JS;
	std::string bz_zty;
};

class OBA
{
public:
	static void init();
	static void getDayName(OB_DAY &r);
	static void getHuiLi(double d0,OB_DAY &r);

private:
	static std::vector<std::vector<std::string>> sFtv; //假日表,由init初始化
	static std::vector<std::string> wFtv;
};

class OBB//农历对象，气朔计算等
{
public:
	static void init();
	static std::string getNH(int y);
	static void getDayName2(OB_DAY &r);
	static void mingLiBaZi(double jd, double J, MLBZ &ob);
	static double qi_accurate(double W);
	static double so_accurate(double W);
	static double qi_accurate2(double jd);
	static double so_accurate2(double jd);
	
private:
	static std::vector<std::string> JNB;
};

extern const char *str_num[];
extern const char *str_ymc[];
extern const char *str_yxmc[];
extern const char *str_jqmc[];
extern const char *str_gan[];
extern const char *str_zhi[];
extern const char *str_sxmc[];
extern const char *str_nywx[];
extern const char *str_xqmc[];
extern const char *str_rmc[];
extern const char *str_rmc0[];
extern const char *str_xz[];
extern const char *str_dx[];
extern const char *str_ago[];
extern const char *str_fw[];
extern const char *str_sjd[];
extern const char *str_ry[];
extern const char *str_ry2[];
extern const char *str_yx[];

#endif