#include "Model.h"
#include"Process.h"
void main(int args, char* argv[]) {
	argv[1] = "../EUstation.cfg";
	string cfgpath = argv[1];
	string StaPath;
	opt_t opt = { 0 };
	gtime_t m_time;
	vector<station_t> station;
	string DefualtDataFile = "../defualt.ionf1";
	string StationListCoor = "../StationList.coor";
	int  obsinv = 30;//建模观测值间隔
	try {
		string str;
		ReadCfgFile(cfgpath, station, StaPath, opt);
		PushAllDataToStr(str, DefualtDataFile, StationListCoor, StaPath, station, opt);
		Process Model(opt, station);
		Model.FindFirEpoch(str);
		for (m_time = opt.StartTime;m_time <= opt.EndTime;m_time.time += opt.ModelInv)
		{
			double centerLL[3] = { 0 };
			vector<ObsSat_t> obs;//存储建模观测数据
			ReadDefaultData(str, obs, m_time, opt);
			ExcludeLessSat(obs);
			calcuCenterLL(obs, station, centerLL);
			map<gtime_t, map<string, vector<ObsSat_t>>> ObsSort;
			ObsSort = VectorSort(obs, station, obsinv);
			Model.ProcessModel(ObsSort, centerLL, opt, m_time);
		}
	}
	catch (char* s) {
		std::cout << s << std::endl;
		getchar();
	}
}