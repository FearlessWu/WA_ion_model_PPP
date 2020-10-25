#pragma once
#include<fstream>
#include<iostream>
#include<string>
#include<cmath>
#include<cassert>
#include<sstream>
#include<set>
#include<vector>
#include<math.h>
#include<map>
#include<algorithm>
#include<exception>
#include<stdio.h>
#include<iomanip>
#include"Matrix.h"
#define LOGFILE							/*log file switch*/
#define SYS_NONE		0x00
#define SYS_GPS			0x01                /* navigation system: GPS */
#define SYS_GLO			0x02                /* navigation system: GLONASS */
#define SYS_GAL			0x04                /* navigation system: Galileo */
#define SYS_BDS			0x08                /* navigation system: BeiDou */
#define GPSMAXSAT		32
#define GLOMAXSAT		27
#define GALMAXSAT		30
#define BDSMAXSAT		35
#define MAXSAT			(GPSMAXSAT+GLOMAXSAT+GALMAXSAT+BDSMAXSAT)
#define PI				3.14159265358979323
#define D2R				(PI/180)
#define M				(1.0/5)				/*卫星数量控制常数*/
#define FE_WGS84		(1.0/298.257223563) /* earth flattening (WGS84)地球偏率 */
#define RE_WGS84		6378137.0           /* earth semimajor axis (WGS84) (m) 地球长半轴*/
#define FREQ1			1.57542E9           /* L1/E1  frequency (Hz) */
#define FREQ1_GLO		1.60200E9           /* GLONASS G1 base frequency (Hz) */
#define DFRQ1_GLO		0.56250E6           /* GLONASS G1 bias frequency (Hz/n) */
#define FREQ1_BDS		1.561098E9          /* BeiDou B1 frequency (Hz) */
#define TEC2D			40.28E16				/*TEC转换常数*/
#define SCALE			(10.0*PI/180)
#define THRESHOLD		0.005				/*收敛阈值*/
#define MinValue		-99.0
#define _CRT_SECURE_NO_WARNINGS
using std::vector;
using std::string;
using std::map;
using std::set;
using std::pair;

typedef struct {        /* time struct */
	time_t time;        /* time (s) expressed by standard time_t */
	double sec;         /* fraction of second under 1 s */
} gtime_t;
typedef struct {

	string sta;
	double coor[3];
	int ID;

}station_t;
class ObsSat_t {
public:
	gtime_t time;
	string sat;
	int prn;
	double el;
	double lat;
	double lon;
	double mf;
	double stec;
	station_t sta;
	bool operator ==(const gtime_t &time);
};
typedef struct {
	gtime_t StartTime;//建模开始时间
	gtime_t EndTime;//建模结束时间
	int sys;//建模卫星系统
	double ModelElv;//建模高度角（degree）
	int ModelInv;//建模的时间间隔（sec）
	int ModelDuration;//建模时间长度（sec）
	int nmax, mmax;//建模阶数
	int ObsInv;
	string logpath;
}opt_t;

typedef struct {
	gtime_t time;
	double y;
	double h;
	double v;
	int sat;
	int ref_sat;
}inerpre_t;

void init(opt_t opt);
bool ReadCfgFile(string cfgpath, vector<station_t> &station, string &StaPath, opt_t &opt);
void Read_all_ionf1_Files(const string StaPath, vector<station_t> &station, const opt_t &opt);
void ReadStaCoor(std::ifstream &stafid, vector<station_t>& station);
void SaveToFlie(vector<ObsSat_t> &obss, vector<station_t> &sta);
void ReadDefaultData(string &str, vector<ObsSat_t> &obss, gtime_t m_time, opt_t opt);
void PushAllDataToStr(string &str, string DefualtDataFile, string StationListCoor, string StaPath, vector<station_t> &station, opt_t opt);
void SortBySat(vector<ObsSat_t> &obs);
void ExcludeLessSat(vector<ObsSat_t> &obs);
void calcuCenterLL(vector<ObsSat_t>&obs, vector<station_t> station, double *centerLL);
map<gtime_t, map<string, vector<ObsSat_t>>> VectorSort(vector<ObsSat_t>&obs, vector<station_t> station, int obsinv);

map<int, int> SelCnstrnSat(pair<string, vector<ObsSat_t>> sta_onepoch);

int SatNum(opt_t opt);
int satsys(int sat);
bool sysexclude(int sys, opt_t opt);
inline bool operator<=(gtime_t &a, gtime_t &b) {
	return a.time <= b.time ? true : false;
}
inline bool ObsSat_t::operator ==(const gtime_t &time) {
	return this->time.time == time.time ? true : false;
}
inline gtime_t operator- (const gtime_t& a, const int inv) {
	gtime_t tmp;
	tmp.sec = a.sec;
	tmp.time = a.time - inv;
	return tmp;
}
inline gtime_t operator+ (const gtime_t& a, const int inv) {
	gtime_t tmp;
	tmp.sec = a.sec;
	tmp.time = a.time + inv;
	return tmp;
}
inline bool operator==(const gtime_t &a, const gtime_t &b) { return (a.time == b.time&&a.sec == b.sec) ? true : false; }
inline bool operator>(const gtime_t &a, const gtime_t &b) { return a.time > b.time ? true : false; }
inline bool operator<(const gtime_t &a, const gtime_t &b) { return a.time < b.time ? true : false; }
inline bool operator>=(const gtime_t &a, const gtime_t &b) { return a.time >= b.time ? true : false; }

const static double gpst0[] = { 1980,1, 6,0,0,0 }; /* gps time reference */
#ifdef __cplusplus
extern "C" {
#endif	
	gtime_t epoch2time(const double *ep);
	double timediff(gtime_t t1, gtime_t t2);
	double dot(const double *a, const double *b, int n);
	void ecef2pos(const double *r, double *pos);
	void time2epoch(gtime_t t, double *ep);
	double time2gpst(gtime_t t, int *week);
#ifdef __cplusplus
}
#endif
