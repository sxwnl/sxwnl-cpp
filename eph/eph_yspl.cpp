#include <cmath>
#include "eph0.h"
#include "eph_yspl.h"

std::array<double,7> ysPL::lT;
std::string ysPL::LX;
double ysPL::sf;

double ysPL::lineT(RE0 G,double v,double u, double r, bool n)
{//已知t1时刻星体位置、速度，求x*x+y*y=r*r时,t的值
  double b=G.y*v-G.x*u, A=u*u+v*v, B=u*b, C=b*b-r*r*v*v, D=B*B-A*C;
  if(D<0) return 0;
  D=sqrt(D); if(!n) D=-D;
  return G.t+((-B+D)/A-G.x)/v;
}
void ysPL::lecXY(double jd, RE0 &re)
{//日月黄经纬差转为日面中心直角坐标(用于月食)
  double T=jd/36525;
  std::array<double, 3> zm={}, zs={};

  //=======太阳月亮黄道坐标========
  zs = e_coord(T,-1,-1,-1);   //地球坐标
  zs[0]  = rad2mrad(zs[0]+_pi+gxc_sunLon(T));  zs[1]  =-zs[1] + gxc_sunLat(T); //补上太阳光行差
  zm = m_coord(T,-1,-1,-1); //月球坐标
  zm[0]  = rad2mrad( zm[0]+gxc_moonLon(T) );  zm[1] += gxc_moonLat(T);  //补上月球光行差就可以了

  //=======视半径=======
  re.e_mRad = cs_sMoon/zm[2]; //月亮地心视半径(角秒)
  re.eShadow = (cs_rEarA/zm[2]*rad-(959.63-8.794)/zs[2] )*51/50; //地本影在月球向径处的半径(角秒),式中51/50是大气厚度补偿
  re.eShadow2= (cs_rEarA/zm[2]*rad+(959.63+8.794)/zs[2] )*51/50; //地半影在月球向径处的半径(角秒),式中51/50是大气厚度补偿

  re.x = rad2rrad(zm[0]+_pi-zs[0]) * cos((zm[1]-zs[1])/2);
  re.y = zm[1]+zs[1];
  
  re.mr= re.e_mRad/rad,  
  re.er=re.eShadow/rad, 
  re.Er=re.eShadow2/rad;
  
  re.t = jd;
}

void ysPL::lecMax(double jd)
{ //月食的食甚计算(jd为近朔的力学时,误差几天不要紧)
  //ysPL::lT={};
  for(int i=0;i<7;i++) ysPL::lT[i]=0; //分别是:食甚,初亏,复圆,半影食始,半影食终,食既,生光
  ysPL::sf=0;
  ysPL::LX="";

  jd = MS_aLon_t2( floor((jd-4)/29.5306)*_pi*2 +_pi)*36525; //低精度的朔(误差10分钟),与食甚相差10分钟左右

  RE0 g={}, G={};
  double u,v;

  //求极值(平均误差数秒)
  u = -18461 * sin(0.057109+0.23089571958*jd)*0.23090/rad; //月日黄纬速度差
  v = (M_v(jd/36525)-E_v(jd/36525))/36525; //月日黄经速度差
  ysPL::lecXY(jd,G);
  jd -= (G.y*u+G.x*v)/(u*u+v*v); //极值时间

  //精密求极值
  double dt=60/86400.0;
  ysPL::lecXY(jd,G); ysPL::lecXY(jd+dt,g); //精密差分得速度,再求食甚
  u = (g.y-G.y)/dt;
  v = (g.x-G.x)/dt;
  dt= -(G.y*u+G.x*v)/(u*u+v*v); jd += dt; //极值时间

  //求直线到影子中心的最小值
  double x=G.x+dt*v, y=G.y+dt*u, rmin=sqrt(x*x+y*y);
  //注意,以上计算得到了极值及最小距rmin,但没有再次计算极值时刻的半径,对以下的判断造成一定的风险,必要的话可以再算一次。不过必要性不很大，因为第一次极值计算已经很准确了,误差只有几秒
  //求月球与影子的位置关系
  if(rmin<=G.mr+G.er){ //食计算
   ysPL::lT[1] = jd; //食甚
   ysPL::LX = "偏";
   ysPL::sf=(G.mr+G.er-rmin)/G.mr/2; //食分

   ysPL::lT[0] = ysPL::lineT(G,v,u, G.mr+G.er, 0); //初亏
   ysPL::lecXY(ysPL::lT[0],g);
   ysPL::lT[0] = ysPL::lineT(g,v,u, g.mr+g.er, 0); //初亏再算一次
//	std::cout<<ysPL::lT[0]<<std::endl;
   ysPL::lT[2] = ysPL::lineT(G,v,u, G.mr+G.er, 1); //复圆
   ysPL::lecXY(ysPL::lT[2],g);
   ysPL::lT[2] = ysPL::lineT(g,v,u, g.mr+g.er, 1); //复圆再算一次
  }
  if(rmin<=G.mr+G.Er){ //半影食计算
   ysPL::lT[3] = ysPL::lineT(G,v,u, G.mr+G.Er, 0); //半影食始
   ysPL::lecXY(ysPL::lT[3],g);
   ysPL::lT[3] = ysPL::lineT(g,v,u, g.mr+g.Er, 0); //半影食始再算一次

   ysPL::lT[4] = ysPL::lineT(G,v,u, G.mr+G.Er, 1); //半影食终
   ysPL::lecXY(ysPL::lT[4],g);
   ysPL::lT[4] = ysPL::lineT(g,v,u, g.mr+g.Er, 1); //半影食终再算一次
  }
  if(rmin<=G.er-G.mr){ //全食计算
   ysPL::LX = "全";
   ysPL::lT[5] = ysPL::lineT(G,v,u, G.er-G.mr, 0); //食既
   ysPL::lecXY(ysPL::lT[5],g);
   ysPL::lT[5] = ysPL::lineT(g,v,u, g.er-g.mr, 0); //食既再算一次

   ysPL::lT[6] = ysPL::lineT(G,v,u, G.er-G.mr, 1); //生光
   ysPL::lecXY(ysPL::lT[6],g);
   ysPL::lT[6] = ysPL::lineT(g,v,u, g.er-g.mr, 1); //生光再算一次
  }
}