#include"Model.h"
extern bool ReadCfgFile(string cfgpath, vector<station_t> &sta, string &StaPath, opt_t &opt) {
	std::ifstream fid(cfgpath);
	if (!fid.is_open()) {
		fid.close();
		throw ("ReadCfgFile error: cannot open cfg file,please check it!");
	}
	string s;
	int sysflag[4] = { 0 };
	int sys[4] = { SYS_GPS,SYS_GLO,SYS_GAL,SYS_BDS };
	double sep[6] = { 0 }, eep[6] = { 0 };
	while (getline(fid, s)) {
		if (string::size_type pos = s.find("FILE_PATH: ", 0) != string::npos)
		{
			StaPath = s.substr(11, string::npos);
		}
		if (string::size_type pos = s.find("Log_Path:	", 0) != string::npos)
		{
			opt.logpath = s.substr(10, string::npos);
		}
		if (string::size_type pos = s.find("UseGPS", 0) != string::npos) {
			sysflag[0] = atoi(s.substr(8, string::npos).c_str());
		}
		if (string::size_type pos = s.find("UseGLO", 0) != string::npos) {
			sysflag[1] = atoi(s.substr(8, string::npos).c_str());
		}
		if (string::size_type pos = s.find("UseGAL", 0) != string::npos) {
			sysflag[2] = atoi(s.substr(8, string::npos).c_str());
		}
		if (string::size_type pos = s.find("UseBDS", 0) != string::npos) {
			sysflag[3] = atoi(s.substr(8, string::npos).c_str());
		}
		if (string::size_type pos = s.find("Model_elv", 0) != string::npos) {
			opt.ModelElv = atof(s.substr(11, string::npos).c_str());
		}
		if (string::size_type pos = s.find("Model_inv", 0) != string::npos) {
			opt.ModelInv = atof(s.substr(11, string::npos).c_str());
		}
		if (string::size_type pos = s.find("nmax", 0) != string::npos) {
			opt.nmax = atoi(s.substr(6, string::npos).c_str());
		}
		if (string::size_type pos = s.find("mmax", 0) != string::npos) {
			opt.mmax = atoi(s.substr(6, string::npos).c_str());
		}
		if (string::size_type pos = s.find("obs_inv", 0) != string::npos) {
			opt.ObsInv = atoi(s.substr(11, string::npos).c_str());
		}
		if (string::size_type pos = s.find("Model_duration", 0) != string::npos) {
			opt.ModelDuration = atoi(s.substr(16, string::npos).c_str());
		}
		if (string::size_type pos = s.find("start_time", 0) != string::npos) {
			sep[0] = atof(s.substr(13, 4).c_str());
			sep[1] = atof(s.substr(18, 2).c_str());
			sep[2] = atof(s.substr(21, 2).c_str());
			sep[3] = atof(s.substr(24, 2).c_str());
			sep[4] = atof(s.substr(27, 2).c_str());
			sep[5] = atof(s.substr(30, 2).c_str());
		}
		if (string::size_type pos = s.find("end_time", 0) != string::npos) {
			eep[0] = atof(s.substr(13, 4).c_str());
			eep[1] = atof(s.substr(18, 2).c_str());
			eep[2] = atof(s.substr(21, 2).c_str());
			eep[3] = atof(s.substr(24, 2).c_str());
			eep[4] = atof(s.substr(27, 2).c_str());
			eep[5] = atof(s.substr(30, 2).c_str());
		}
		if (string::size_type pos = s.find("STATION_LIST:", 0) != string::npos)
		{
			while (fid)
			{
				string sat;
				int StaID = 0;
				string StaName;
				fid >> sat;
				if (string::size_type a = sat.find("#", 0) == string::npos)
					StaID = atof(sat.substr(0, string::npos).c_str());
				fid >> StaName;
				if (StaID != 0) {
					station_t s;
					s.sta = StaName;
					s.ID = StaID;
					sta.push_back(s);
				}
			}
		}

	}
	opt.StartTime = epoch2time(sep);
	opt.EndTime = epoch2time(eep);
	for (int i = 0;i < 4;++i) {
		if (sysflag[i] == 1)
			opt.sys |= sys[i];
	}
	fid.close();
	return true;
}
static inline bool CmpTime(ObsSat_t&a, ObsSat_t&b) {
	return (int)a.time.time< (int)b.time.time;
}
void ReadStaCoor(std::ifstream &fid, vector<station_t>& station) {
	string s;
	station.clear();
	while (getline(fid, s))
	{
		station_t sta;
		sta.sta = s.substr(0, 4);
		sta.ID = atoi(s.substr(5, 3).c_str());
		sta.coor[0] = atof(s.substr(9, 15).c_str());
		sta.coor[1] = atof(s.substr(25, 15).c_str());
		sta.coor[2] = atof(s.substr(41, 15).c_str());
		station.push_back(sta);
	}
}
extern void ReadDefaultData(string &str, vector<ObsSat_t> &obss, gtime_t m_time, opt_t opt) {
	printf("\nreading data ......\n");
	std::stringstream sstr;
	string time;
	string::size_type time_pos;
	size_t size;
	string::size_type pos;
	sstr << (int)(m_time.time - opt.ModelDuration / 2);
	time = sstr.str();
	size = str.size();
	time_pos = str.find(time, 0);
	if (time_pos == string::npos) {
		throw("ReadDefaultData error: 无此时段的建模数据！");
	}
	pos = str.find("\n", 0);
	int inv = pos + 1;
	for (int i = (time_pos - 1);i<size && (i + inv) <= size;i += inv)
	{
		string s;
		s = str.substr(i, inv);
		ObsSat_t obs = { 0 };
		obs.time.time = atoi(s.substr(0, 11).c_str());
		if (obs.time.time>(m_time.time + opt.ModelDuration / 2))
			break;
		obs.sat = s.substr(12, 3);
		obs.prn = atoi(s.substr(16, 3).c_str());
		obs.el = atof(s.substr(20, 8).c_str());
		obs.lat = atof(s.substr(29, 8).c_str());
		obs.lon = atof(s.substr(38, 9).c_str());
		obs.mf = atof(s.substr(48, 8).c_str());
		obs.stec = atof(s.substr(57, 12).c_str());
		obs.sta.sta = s.substr(70, 4);
		obs.sta.ID = atoi(s.substr(75, 2).c_str());
		if (obs.el < opt.ModelElv*D2R)//低于截止高度角的不存储
			continue;
		int sys = satsys(obs.prn);
		if (!sysexclude(sys, opt))//排除卫星系统
			continue;
		if (fabs(obs.stec)>300)
			continue;
		obss.push_back(obs);

	}
	//SortBySat(obss);
	printf("\nreading data finished\n");
}
extern void SaveToFlie(vector<ObsSat_t> &obss, vector<station_t> &sta) {
	using namespace std;
	string SavePath = "../defualt.ionf1";
	ofstream fout;
	fout.open(SavePath);
	vector<ObsSat_t>::iterator iter;
	for (iter = obss.begin();iter != obss.end();++iter)
	{
		fout << setw(11) << iter->time.time << " " << iter->sat << " " << setw(3) << iter->prn << " " << setw(8) << iter->el << " " << setw(8) << iter->lat << " " << setw(9) << iter->lon << " " << setw(8) << iter->mf << " "
			<< setw(12) << iter->stec << " " << iter->sta.sta << " " << setw(2) << iter->sta.ID << endl;
	}
	fout.close();
	string SaveStaPath = "../StationList.coor";
	ofstream outf;
	outf.open(SaveStaPath);
	outf << fixed;
	for (vector<station_t>::iterator iter = sta.begin();iter != sta.end();++iter)
	{
		outf << iter->sta << " " << setfill('0') << setw(3) << iter->ID << " " << setfill(' ') << setw(15) << (double)iter->coor[0] << " " << setw(15) << (double)iter->coor[1] << " " << setw(15) << (double)iter->coor[2] << endl;
	}
	outf.close();
}
extern void Read_all_ionf1_Files(const string StaPath, vector<station_t> &sta, const opt_t &opt) {
	clock_t read_start = clock();
	std::cout << " reading data......" << std::endl;
	vector<ObsSat_t> obss;
	vector<station_t>::iterator iter;
	int i;
	for (iter = sta.begin();iter != sta.end();++iter) {
		std::cout << " reading " << iter->sta << " station file......" << std::endl;
		i = 0;
		string StaFullPath = StaPath + iter->sta + ".ionf1";
		std::ifstream Stafile(StaFullPath);
		double ep[6];
		double coord[3];
		string s;
		string StaName;
		int satnum;
		if (!Stafile.is_open()) {
			std::cout << "\n WARNNING： CANNOT OPEN " << iter->sta << " FILE,PLEASE CHECK IT!\n" << std::endl;
			continue;
		}
		while (getline(Stafile, s)) {
			string::size_type pos;
			if (pos = s.find("STA_NAME: ", 0) != string::npos)
			{
				StaName = s.substr(10, string::npos);
			}
			if (pos = s.find("ANT_TYPE:", 0) != string::npos)
			{
				string whatever;
				Stafile >> whatever;
				Stafile >> coord[0];
				Stafile >> coord[1];
				Stafile >> coord[2];
				//coord[0] = atof(s.substr(10, 15).c_str());
				//coord[1] = atof(s.substr(26, 15).c_str());
				//coord[2] = atof(s.substr(42, 15).c_str());
			}
			if (pos = s.find("> ", 0) != string::npos) {
				ObsSat_t temp;
				temp.sta = (*iter);
				ep[0] = atof(s.substr(2, 6).c_str());
				ep[1] = atof(s.substr(7, 2).c_str());
				ep[2] = atof(s.substr(10, 2).c_str());
				ep[3] = atof(s.substr(13, 2).c_str());
				ep[4] = atof(s.substr(16, 2).c_str());
				ep[5] = atof(s.substr(19, 11).c_str());
				satnum = atof(s.substr(34, string::npos).c_str());
				temp.time = epoch2time(ep);
				for (int i = 0;i < satnum;++i) {
					Stafile >> temp.sat;
					Stafile >> temp.el;
					Stafile >> temp.lat;
					Stafile >> temp.lon;
					Stafile >> temp.mf;
					Stafile >> temp.stec;
					temp.prn = atoi(temp.sat.substr(1, string::npos).c_str());
					if (temp.sat[0] == 'R')
						temp.prn += 32;
					if (temp.sat[0] == 'E')
						temp.prn += 59;
					if (temp.sat[0] == 'C')
						temp.prn += 89;
					//if (temp.el < opt.ModelElv*D2R)//低于截止高度角的不存储
					//	continue;
					//int sys = satsys(temp.prn);
					//if (!sysexclude(sys, opt))//排除卫星系统
					//	continue;
					obss.push_back(temp);
				}

			}

		}
		for (int i = 0;i < 3;++i)
			iter->coor[i] = coord[i];
		i++;
	}
	std::cout << "\n sorting......\n";
	sort(obss.begin(), obss.end(), CmpTime);
	printf("\n saving data to defualt file......\n\n");
	SaveToFlie(obss, sta);
	clock_t read_end = clock();
	double read_dur = read_end - read_start;
	printf("\n 读文件共用时%.2f sec\n", read_dur / 1000);
}
extern void PushAllDataToStr(string &str, string DefualtDataFile, string StationListCoor, string StaPath, vector<station_t> &station, opt_t opt)
{
	std::ifstream fid(DefualtDataFile, std::ios::in);
	std::ifstream stafid(StationListCoor);
	
	if (!fid.is_open() || !stafid.is_open())
	{
		fid.close();
		Read_all_ionf1_Files(StaPath, station, opt);
		std::ifstream fid1(DefualtDataFile, std::ios::in);
		std::stringstream buf;
		buf << fid1.rdbuf();
		str = buf.str();
	}
	else
	{
		std::stringstream buf;
		buf << fid.rdbuf();
		str = buf.str();
		ReadStaCoor(stafid, station);
	}
	fid.close();
}