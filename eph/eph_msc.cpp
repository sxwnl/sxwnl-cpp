#include <cmath>

#include "eph0.h"
#include "eph.h"
#include "../tool.h"
#include "eph_msc.h"

double msc::T;//TD力学时
double msc::L;//地理纬度
double msc::fa;//地理经度
double msc::dt;//力学-世界时时差
double msc::jd;//UT世界时
double msc::dL;//黄经章动
double msc::dE;//黄纬章动
double msc::E;//交角章动
double msc::gst;//真恒星时

double msc::mHJ;//月球视黄经
double msc::mHW;//月球视黄纬
double msc::mR;//地月质心距
double msc::mCJ;//月球视赤经
double msc::mCW;//月球视赤纬
double msc::mShiJ;//月球时角

double msc::mCJ2;//时差修正后的赤道坐标
double msc::mCW2;
double msc::mR2;
double msc::mDJ;//高度角
double msc::mDW;//方位角
double msc::mPJ;//大气折射修正后的高度角
double msc::mPW;//大气折射修正后的方位角

double msc::sHJ;//太阳视黄经
double msc::sHW;//太阳视黄纬
double msc::sCJ;//太阳视赤经
double msc::sCW;//太阳视赤纬
double msc::sCJ2;//时差修正后的赤道坐标
double msc::sCW2;
double msc::sR2;
double msc::sShiJ;//太阳时角

double msc::sDJ;//高度角
double msc::sDW;//方位角
double msc::sR;
double msc::sPJ;//方位角
double msc::sPW;//高度角
double msc::sc;//时差

double msc::pty;//平恒星时
double msc::zty;//真恒星时
double msc::mRad;//月亮视半径
double msc::sRad;//太阳视半径
double msc::e_mRad;//月亮地心视半径
double msc::eShadow;//地本影在月球向径处的半径(角秒)
double msc::eShadow2;//地半影在月球向径处的半径(角秒)
double msc::mIll;//月面被照亮比例
double msc::zx_J;//中心食坐标
double msc::zx_W;

void msc::calc(double T, double L, double fa, double high)
{								//sun_moon类的成员函数。参数：T是力学时,站点经纬L,fa,海拔high(千米)
	//基本参数计算
	msc::T = T, msc::L = L, msc::fa = fa;
	msc::dt = dt_T(T);			//TD-UT
	msc::jd = T - msc::dt;	//UT
	T /= 36525.0;
	std::array<double, 2> zd = nutation2(T);
	msc::dL = zd[0];			//黄经章
	msc::dE = zd[1];			//交角章动
	msc::E = hcjj(T) + msc::dE;	//真黄赤交角
	msc::gst = pGST(msc::jd, msc::dt) + msc::dL * cos(msc::E);	//真恒星时(不考虑非多项式部分)
	std::array<double, 3> z;

	//=======月亮========
	//月亮黄道坐标
	z = m_coord(T, -1, -1, -1);	//月球坐标
	z[0] = rad2mrad(z[0] + gxc_moonLon(T) + msc::dL);
	z[1] += gxc_moonLat(T);		//补上月球光行差及章动
	msc::mHJ = z[0];
	msc::mHW = z[1];
	msc::mR = z[2];			//月球视黄经,视黄纬,地月质心距

	//月球赤道坐标
	z = llrConv(z, msc::E);	//转为赤道坐标
	msc::mCJ = z[0];
	msc::mCW = z[1];			//月球视赤经,月球赤纬

	//月亮时角计算
	msc::mShiJ = rad2mrad(msc::gst + L - z[0]);	//得到此刻天体时角
	if (msc::mShiJ > _pi)
		msc::mShiJ -= pi2;

	//修正了视差的赤道坐标
	z = parallax(z, msc::mShiJ, fa, high);	//视差修正
	msc::mCJ2 = z[0], msc::mCW2 = z[1], msc::mR2 = z[2];

	//月亮时角坐标
	z[0] += _pi / 2 - msc::gst - L;	//转到相对于地平赤道分点的赤道坐标(时角坐标)

	//月亮地平坐标
	z = llrConv(z, _pi / 2 - fa);	//转到地平坐标(只改经纬度)
	z[0] = rad2mrad(_pi / 2 - z[0]);
	msc::mDJ = z[0];
	msc::mDW = z[1];			//方位角,高度角
	if (z[1] > 0)
		z[1] += MQC(z[1]);		//大气折射修正
	msc::mPJ = z[0];
	msc::mPW = z[1];			//方位角,高度角

	//=======太阳========
	//太阳黄道坐标
	z = e_coord(T, -1, -1, -1);	//地球坐标
	z[0] = rad2mrad(z[0] + _pi + gxc_sunLon(T) + msc::dL);	//补上太阳光行差及章动
	z[1] = -z[1] + gxc_sunLat(T);	//z数组为太阳地心黄道视坐标
	msc::sHJ = z[0];
	msc::sHW = z[1];
	msc::sR = z[2];			//太阳视黄经,视黄纬,日地质心距

	//太阳赤道坐标
	z = llrConv(z, msc::E);	//转为赤道坐标
	msc::sCJ = z[0];
	msc::sCW = z[1];			//太阳视赤经,视赤纬

	//太阳时角计算
	msc::sShiJ = rad2mrad(msc::gst + L - z[0]);	//得到此刻天体时角
	if (msc::sShiJ > _pi)
		msc::sShiJ -= pi2;

	//修正了视差的赤道坐标
	z = parallax(z, msc::sShiJ, fa, high);	//视差修正
	msc::sCJ2 = z[0], msc::sCW2 = z[1], msc::sR2 = z[2];

	//太阳时角坐标
	z[0] += _pi / 2 - msc::gst - L;	//转到相对于地平赤道分点的赤道坐标

	//太阳地平坐标
	z = llrConv(z, _pi / 2 - fa);
	z[0] = rad2mrad(_pi / 2 - z[0]);
	//z[1] -= 8.794/rad/z[2]*cos(z[1]); //直接在地平坐标中视差修正(这里把地球看为球形,精度比parallax()稍差一些)
	msc::sDJ = z[0];
	msc::sDW = z[1];			//方位角,高度角

	if (z[1] > 0)
		z[1] += MQC(z[1]);		//大气折射修正
	msc::sPJ = z[0];
	msc::sPW = z[1];			//方位角,高度角

	//=======其它========
	//时差计算
	double t = T / 10, t2 = t * t, t3 = t2 * t, t4 = t3 * t, t5 = t4 * t;
	double Lon = (1753470142 + 6283319653318 * t + 529674 * t2 + 432 * t3 - 1124 * t4 - 9 * t5) / 1000000000 + _pi - 20.5 / rad;	//修正了光行差的太阳平黄经
	Lon = rad2mrad(Lon - (msc::sCJ - msc::dL * cos(msc::E)));	//(修正了光行差的平黄经)-(不含dL*cos(E)的视赤经)
	if (Lon > _pi)
		Lon -= pi2;				//得到时差,单位是弧度
	msc::sc = Lon / pi2;		//时差(单位:日)

	//真太阳与平太阳
	msc::pty = msc::jd + L / pi2;	//平太阳时
	msc::zty = msc::jd + L / pi2 + msc::sc;	//真太阳时

	//视半径
	// msc::mRad = moonRad(msc::mR,msc::mDW);  //月亮视半径(角秒)
	msc::mRad = cs_sMoon / msc::mR2;	//月亮视半径(角秒)
	msc::sRad = 959.63 / msc::sR2;	//太阳视半径(角秒)
	msc::e_mRad = cs_sMoon / msc::mR;	//月亮地心视半径(角秒)
	msc::eShadow = (cs_rEarA / msc::mR * rad - (959.63 - 8.794) / msc::sR) * 51 / 50;	//地本影在月球向径处的半径(角秒),式中51/50是大气厚度补偿
	msc::eShadow2 = (cs_rEarA / msc::mR * rad + (959.63 + 8.794) / msc::sR) * 51 / 50;	//地半影在月球向径处的半径(角秒),式中51/50是大气厚度补偿
	msc::mIll = moonIll(T);	//月亮被照面比例

	//中心食计算
	if (fabs(rad2rrad(msc::mCJ - msc::sCJ)) < 50.0 / 180.0 * _pi)
	{
		COORDP pp = lineEar({ msc::mCJ, msc::mCW, msc::mR }
							, { msc::sCJ, msc::sCW, msc::sR * cs_AU }
							, msc::gst);
		msc::zx_J = pp.J;
		msc::zx_W = pp.W;		//无解返回值是100
	}
	else
		msc::zx_J = msc::zx_W = 100;
}

std::string msc::toStr(bool fs)
{
	std::string s;
	s = "-------------------------------------------\n";
	s = s + "平太阳 " + timeStr(msc::pty) + " 真太阳 " + timeStr(msc::zty) + "\n";
	s = s + "时差 " + m2fm(msc::sc * 86400, 2, 1) + " 月亮被照亮 " + to_str(msc::mIll * 100, 2) + "% ";
	s = s + "\n";

	s = s + "-------------------------------------------\n表一       月亮            太阳\n";
	s = s + "视黄经 " + rad2str(msc::mHJ, 0) + "  " + rad2str(msc::sHJ, 0) + "\n";
	s = s + "视黄纬 " + rad2str(msc::mHW, 0) + "  " + rad2str(msc::sHW, 0) + "\n";
	s = s + "视赤经 " + rad2str(msc::mCJ, 1) + "  " + rad2str(msc::sCJ, 1) + "\n";
	s = s + "视赤纬 " + rad2str(msc::mCW, 0) + "  " + rad2str(msc::sCW, 0) + "\n";
	s = s + "距离    " + to_str(msc::mR, 2) + "千米    " + to_str(msc::sR, 8) + "AU" + "\n";

	s = s + "-------------------------------------------\n表二       月亮            太阳\n";
	s = s + "方位角 " + rad2str(msc::mPJ, 0) + "  " + rad2str(msc::sPJ, 0) + "\n";
	s = s + "高度角 " + rad2str(msc::mPW, 0) + "  " + rad2str(msc::sPW, 0) + "\n";
	s = s + "时角   " + rad2str(msc::mShiJ, 0) + "  " + rad2str(msc::sShiJ, 0) + "\n";
	s = s + "视半径   " + m2fm(msc::mRad, 2, 0) + "       " + m2fm(msc::sRad, 2, 0) + " (观测点)\n";

	if (fs)
	{
		s = s + "-------------------------------------------\n";
		s = s + "力学时" + JD2str(msc::T + J2000);
		s = s + " ΔT=" + to_str(msc::dt * 86400, 1) + "秒\n";
		s = s + "黄经章 " + to_str(msc::dL / pi2 * 360 * 3600, 2) + "\" ";
		s = s + "交角章 " + to_str(msc::dE / pi2 * 360 * 3600, 2) + "\" ";
		s = s + "ε=" + rad2str(msc::E, 0) + "\n";
		s = s + "-------------------------------------------\n";

	}
	return s;
}