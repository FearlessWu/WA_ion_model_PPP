#pragma once
#include"Model.h"

class Process {
public:
	Process(opt_t opt, vector<station_t> station);
	~Process();
	map<int, pair<int, double>> SelRefSat(pair<string, vector<ObsSat_t>> sta_onepoch);
	void ProcessModel(map<gtime_t, map<string, vector<ObsSat_t>>> &ObsSort, double *centerLL, opt_t opt, gtime_t m_time);
	void FindFirEpoch(string &str);
private:
	bool init;
	int iterflag;//是否首次迭代flag
	int inerflag;//内符合精度flag
	int cnstflag;//卫星约束flag
	gtime_t FirstEpoch;
	gtime_t LastEpoch;
	map<string, std::ofstream&> ion;
	map<string, std::ofstream&> ionfit;
	std::ofstream *fp;
	std::ofstream *fpfit;
	vector<station_t> sta;
	FILE*log;
	std::ofstream Enm_out;
	map<int, string> SYS;
	map<gtime_t, map<string, vector<ObsSat_t>>> DurData;
	opt_t opt;
private:
	
	void NorEquSup(pair<string, vector<ObsSat_t>> sta_onepoch, map<int, pair<int, double>> &ref_sat, double *centerLL, opt_t opt, gtime_t m_time, Matrix& HH, Matrix& Hy, Matrix &x, map<int, pair<int, double>> CnstrnSat);
	void AssemblyDesignMatrix(Matrix& B, Matrix &x, Matrix &h, ObsSat_t &sat, ObsSat_t &RefObs, double* centerLL, opt_t opt, gtime_t &m_time, int row);
	void OutputFunc(gtime_t m_time, Matrix &x, double *centerLL, opt_t opt);
	void WashX(Matrix &x, opt_t opt);
	void OutputInerPre(string StaName, vector<inerpre_t>& pre);
};
