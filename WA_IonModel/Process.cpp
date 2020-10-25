#include"Model.h"
#include"Process.h"
static inline bool CmpSat(ObsSat_t&a, ObsSat_t&b) {
	return a.prn<b.prn ? true : false;
}
//按卫星号排序
extern void SortBySat(vector<ObsSat_t> &obs)
{
	sort(obs.begin(), obs.end(), CmpSat);

}
//排除数目过少卫星
extern void ExcludeLessSat(vector<ObsSat_t> &obs) {
	std::cout << "\nExcluding Sat ......" << std::endl;
	vector<ObsSat_t> gps, glo, gal, bds;
	vector<ObsSat_t>::iterator iter;
	for (iter = obs.begin();iter != obs.end();++iter) {
		if (satsys(iter->prn) == SYS_GPS)
			gps.push_back(*iter);
		if (satsys(iter->prn) == SYS_GLO)
			glo.push_back(*iter);
		if (satsys(iter->prn) == SYS_GAL)
			gal.push_back(*iter);
		if (satsys(iter->prn) == SYS_BDS)
			bds.push_back(*iter);
	}
	vector<vector<ObsSat_t>> gnss;
	map<int, int> sat_arr;
	map<int, int>::iterator itf;
	gnss.push_back(gps);gnss.push_back(glo);gnss.push_back(gal);gnss.push_back(bds);
	for (vector<ObsSat_t> t : gnss) {
		for (ObsSat_t s : t) {
			if (!sat_arr.count(s.prn))
				sat_arr.insert(pair<int, int>(s.prn, 0));
			else
			{
				itf = sat_arr.find(s.prn);
				itf->second++;
			}

		}
	}
	double gps_sum[3] = { 0 }, glo_sum[3] = { 0 }, gal_sum[3] = { 0 }, bds_sum[3] = { 0 };
	for (pair<int, double> s : sat_arr) {
		if (satsys(s.first) == SYS_GPS)
		{
			gps_sum[0] += s.second;
			gps_sum[1]++;
		}
		if (satsys(s.first) == SYS_GLO)
		{
			glo_sum[0] += s.second;
			glo_sum[1]++;
		}
		if (satsys(s.first) == SYS_GAL)
		{
			gal_sum[0] += s.second;
			gal_sum[1]++;
		}
		if (satsys(s.first) == SYS_BDS)
		{
			bds_sum[0] += s.second;
			bds_sum[1]++;
		}

	}

	gps_sum[2] = gps_sum[0] / gps_sum[1] * M;
	glo_sum[2] = glo_sum[0] / glo_sum[1] * M;
	gal_sum[2] = gal_sum[0] / gal_sum[1] * M;
	bds_sum[2] = bds_sum[0] / bds_sum[1] * M;
	set<int> exsat;
	for (pair<int, double> s : sat_arr) {
		if (satsys(s.first) == SYS_GPS)
		{
			if (s.second < gps_sum[2])
				exsat.insert(s.first);
		}
		if (satsys(s.first) == SYS_GLO)
		{
			if (s.second < glo_sum[2])
				exsat.insert(s.first);
		}
		if (satsys(s.first) == SYS_GAL)
		{
			if (s.second < gal_sum[2])
				exsat.insert(s.first);
		}
		if (satsys(s.first) == SYS_BDS)
		{
			if (s.second < bds_sum[2])
				exsat.insert(s.first);
		}

	}
	if (exsat.size() == 0)
		return;
	for (vector<ObsSat_t>::iterator iter = obs.begin();iter != obs.end();) {
		if (exsat.count(iter->prn))
			iter = obs.erase(iter);
		else
			iter++;
	}
}
//计算建模中心经纬度
extern void calcuCenterLL(vector<ObsSat_t>&obs, vector<station_t> station, double *centerLL) {
	using namespace std;
	std::cout << "\nCalculating center lat and lon .....\n";
	map<string, int> stalist;
	for (auto o : obs) {
		if (!stalist.count(o.sta.sta)) {
			stalist.insert(pair<string, int>(o.sta.sta, 1));
		}
		else {
			auto iter = stalist.find(o.sta.sta);
			iter->second++;
		}
	}
	double ref_xyz[3] = { 0 };
	double ref_llh[3] = { 0 };
	vector<station_t>::iterator iter;
	for (iter = station.begin();iter != station.end();++iter) {
		if (!stalist.count(iter->sta))
			continue;
		else {
			auto s = stalist.find(iter->sta);
			for (int i = 0;i < 3;++i)
			{
				//ref_xyz[i] += iter->coor[i] * s->second;
				ref_xyz[i] = iter->coor[i];

			}
			ecef2pos(ref_xyz, ref_llh);
			//cout << iter->sta << ": " << ref_llh[0] << "  " << ref_llh[1] << "  " << ref_llh[2] << endl;
			for (int i = 0;i < 3;++i)
			{
				ref_llh[i] *= s->second;
				centerLL[i] += ref_llh[i];

			}
			//std::cout<<iter->sta<<": "<<setw(4)<<cent[0]<<" "<<setw(4)<<cent[1]<<" "<<setw(4)<<cent[2]<<"\n";
		}

	}
	for (int i = 0;i < 3;++i)
		centerLL[i] /= obs.size();
	//ecef2pos(ref_xyz, centerLL);
}
//将建模时长内的数据按以下索引：历元-->站-->卫星观测值
extern map<gtime_t, map<string, vector<ObsSat_t>>> VectorSort(vector<ObsSat_t>&obs, vector<station_t> station, int obsinv) {
	map<gtime_t, map<string, vector<ObsSat_t>>> ObsSort;
	for (gtime_t epoch = obs.begin()->time;epoch.time <= obs.rbegin()->time.time;epoch.time += obsinv)
	{
		vector<ObsSat_t> Onepoch;
		vector<ObsSat_t>::iterator iter;
		int flag = 0;
		for (iter = obs.begin();iter != obs.end();++iter) {
			if (iter->time == epoch)
			{
				Onepoch.push_back(*iter);
				flag = 1;
			}
			else if (iter->time > epoch&&flag == 1)
				break;
		}
		map<string, vector<ObsSat_t>> sta_onepoch;
		for (vector<ObsSat_t>::iterator iter = Onepoch.begin();iter != Onepoch.end();++iter) {
			if (!sta_onepoch.count(iter->sta.sta))
			{
				vector<ObsSat_t> temp;
				temp.push_back(*iter);
				sta_onepoch.insert(pair<string, vector<ObsSat_t>>(iter->sta.sta, temp));
			}
			else {

				map<string, vector<ObsSat_t>>::iterator it = sta_onepoch.find(iter->sta.sta);
				it->second.push_back(*iter);
			}

		}
		ObsSort.insert(pair<gtime_t, map<string, vector<ObsSat_t>>>(epoch, sta_onepoch));


	}
	return ObsSort;
}
//构造函数
Process::Process(opt_t _opt, vector<station_t> station) {
	using namespace std;
	sta = station;
	opt = _opt;
	string logpath = opt.logpath + "log.txt";
	if (!(log = fopen(logpath.c_str(), "w")))
		throw("Log file open error: log file open failed, please close the opening log file!");

	string Enmpath = opt.logpath + "Enm.txt";
	Enm_out.open(Enmpath);
	if (!Enm_out.is_open())
		throw("Enm file open error: Enm file open failed, please close the opening Enm file!");

	fp = nullptr;
	fpfit = nullptr;
	if (station.size() != 0)
	{
		fp = new ofstream[station.size()];
		fpfit = new ofstream[station.size()];
		vector<station_t>::iterator iter = station.begin();
		for (int i = 0;i<station.size();++i)
		{
			string name = opt.logpath + "inerpre\\" + iter->sta + ".ion";
			string namefit = opt.logpath + "inerpre\\" + iter->sta + "_fit.ion";
			fp[i].open(name);
			fpfit[i].open(namefit);
			if (fp[i].is_open() && fpfit[i].is_open())
			{
				ion.insert(pair<string, ofstream&>(iter->sta, fp[i]));
				ionfit.insert(pair<string, ofstream&>(iter->sta, fpfit[i]));
			}
			else
				fprintf(log, "%s 内符合文件打开失败！\n", iter->sta.c_str());
			iter++;
		}

	}

	SYS.insert(pair<int, string>(SYS_GPS, "GPS"));
	SYS.insert(pair<int, string>(SYS_GLO, "GLO"));
	SYS.insert(pair<int, string>(SYS_GAL, "GAL"));
	SYS.insert(pair<int, string>(SYS_BDS, "BDS"));

}
//析构函数
Process::~Process() {
	if (log)
		fclose(log);
	if (Enm_out.is_open())
		Enm_out.close();
	SYS.clear();
	for (map<string, std::ofstream&>::iterator iter = ion.begin();iter != ion.end();++iter)
	{
		if (iter->second.is_open())
			iter->second.close();
	}
	if (fp != nullptr)
		delete[]fp;
	ion.clear();

	for (map<string, std::ofstream&>::iterator iter = ionfit.begin();iter != ionfit.end();++iter)
	{
		if (iter->second.is_open())
			iter->second.close();
	}
	if (fpfit != nullptr)
		delete[]fpfit;
	ionfit.clear();

	sta.clear();
}
//清洗待估参数X矩阵
void Process::WashX(Matrix &x, opt_t opt)
{
	int ntec = (opt.nmax + 1)*(opt.mmax + 1);
	double NumBase = 0;
	if (SYS_GPS&opt.sys)
	{
		for (int i = ntec;i<ntec + GPSMAXSAT;++i)
		{
			if (x[i] != 0.0)
			{
				NumBase = x[i] - 1;
				break;
			}
		}
		for (int i = ntec;i<ntec + GPSMAXSAT;++i)
		{
			if (x[i] != 0.0)
			{
				x[i] -= NumBase;
			}
		}
	}

	if (SYS_GLO&opt.sys)
	{
		for (int i = ntec + GPSMAXSAT;i<ntec + GPSMAXSAT + GLOMAXSAT;++i)
		{
			if (x[i] != 0.0)
				NumBase = x[i] - 1;
			break;
		}
		for (int i = ntec + GPSMAXSAT;i<ntec + GPSMAXSAT + GLOMAXSAT;++i)
		{
			if (x[i] != 0.0)
				x[i] -= NumBase;
		}
	}

	if (SYS_GAL&opt.sys)
	{
		for (int i = ntec + GPSMAXSAT + GLOMAXSAT;i<ntec + GPSMAXSAT + GLOMAXSAT + GALMAXSAT;++i)
		{
			if (x[i] != 0.0)
				NumBase = x[i] - 1;
			break;
		}
		for (int i = ntec + GPSMAXSAT + GLOMAXSAT;i<ntec + GPSMAXSAT + GLOMAXSAT + GALMAXSAT;++i)
		{
			if (x[i] != 0.0)
				x[i] -= NumBase;
		}
	}

	if (SYS_BDS&opt.sys)
	{
		for (int i = ntec + GPSMAXSAT + GLOMAXSAT + GALMAXSAT;i<ntec + MAXSAT;++i)
		{
			if (x[i] != 0.0)
				NumBase = x[i] - 1;
			break;
		}
		for (int i = ntec + GPSMAXSAT + GLOMAXSAT + GALMAXSAT;i<ntec + MAXSAT;++i)
		{
			if (x[i] != 0.0)
				x[i] -= NumBase;
		}
	}
}
//待估参数输出
void Process::OutputFunc(gtime_t m_time, Matrix &x, double *centerLL, opt_t opt) {
	using namespace std;
	int week, sec;
	Enm_out << fixed;
	sec = (int)time2gpst(m_time, &week);
	int ntec = (opt.mmax + 1)*(opt.nmax + 1);
	if (Enm_out.is_open())
	{
		Enm_out << "# " << setw(4) << week << " " << setw(6) << sec << " 1 " << setw(14) << centerLL[0] << " " << setw(14) << centerLL[1] << "          0.00000" << endl;
		Enm_out << "#TEC ";
		for (int i = 0;i<ntec;i++)
		{
			Enm_out << " " << setw(9) << x[i];
		}
		Enm_out << endl;

		if (opt.sys&SYS_GPS)
		{
			Enm_out << "#GPS ";
			for (int i = 0;i<GPSMAXSAT;i++)
			{
				Enm_out << " " << setw(9) << x[ntec + i];
			}
			Enm_out << endl;
		}
		if (opt.sys&SYS_GLO)
		{
			Enm_out << "#GLO ";
			for (int i = 0;i<GLOMAXSAT;i++)
			{
				Enm_out << " " << setw(9) << x[ntec + GPSMAXSAT + i];
			}
			Enm_out << endl;
		}
		if (opt.sys&SYS_GAL)
		{
			Enm_out << "#GAL ";
			for (int i = 0;i<GALMAXSAT;i++)
			{
				Enm_out << " " << setw(9) << x[ntec + GPSMAXSAT + GLOMAXSAT + i];
			}
			Enm_out << endl;
		}
		if (opt.sys&SYS_BDS)
		{
			Enm_out << "#BDS ";
			for (int i = 0;i<BDSMAXSAT;i++)
			{
				Enm_out << " " << setw(9) << x[ntec + GPSMAXSAT + GLOMAXSAT + GALMAXSAT + i];
			}
			Enm_out << endl;
		}

	}
}
//选取参考星
map<int, pair<int, double>> Process::SelRefSat(pair<string, vector<ObsSat_t>> sta_onepoch) {
	map<int, pair<int, double>> ref;
	int ref_sat[4] = { 0 };
	double max_el[4] = { 0 };
	int i;
	//选取高度角最高的为参考星
	for (auto sat : sta_onepoch.second)
	{//sat	-->取一个历元一个站的一颗卫星数据
		int sys = satsys(sat.prn);
		switch (sys) {
		case SYS_GPS:
			i = 0;
			break;
		case SYS_GLO:
			i = 1;
			break;
		case SYS_GAL:
			i = 2;
			break;
		case SYS_BDS:
			i = 3;
			break;
		default:
			throw("卫星号不正确！");
		}
		if (sat.el > max_el[i])
		{
			max_el[i] = sat.el;
			ref_sat[i] = sat.prn;
		}
	}

	ref.insert(pair<int, pair<int, double>>(SYS_GPS, pair<int, double>(ref_sat[0], max_el[0])));
	ref.insert(pair<int, pair<int, double>>(SYS_GLO, pair<int, double>(ref_sat[1], max_el[1])));
	ref.insert(pair<int, pair<int, double>>(SYS_GAL, pair<int, double>(ref_sat[2], max_el[2])));
	ref.insert(pair<int, pair<int, double>>(SYS_BDS, pair<int, double>(ref_sat[3], max_el[3])));
	return ref;
}
//输出内符合文件
void Process::OutputInerPre(string StaName, vector<inerpre_t>& pre)
{
	using namespace std;
	map<string, ofstream&>::iterator iter = ion.find(StaName);
	map<string, ofstream&>::iterator iterfit = ionfit.find(StaName);
	if (iter != ion.end() && iterfit != ionfit.end())
	{
		//ion 05 21 2018/11/04 03:00:00  0.34074
		double ep[6] = { 0 };
		time2epoch(pre.begin()->time, ep);
		for (vector<inerpre_t>::iterator itpre = pre.begin();itpre != pre.end();++itpre)
		{
			iter->second << "ion " << setfill('0') << setw(2) << itpre->sat << " " << itpre->ref_sat << " " << setw(4) << (int)ep[0] << "/" << setw(2) << (int)ep[1] << "/"
				<< setw(2) << (int)ep[2] << " " << setw(2) << (int)ep[3] << ":" << setw(2) << (int)ep[4] << ":" << setw(2) << (int)ep[5] << " " << setfill(' ') << setw(8) << itpre->y << endl;

			iterfit->second << "ion " << setfill('0') << setw(2) << itpre->sat << " " << itpre->ref_sat << " " << setw(4) << (int)ep[0] << "/" << setw(2) << (int)ep[1] << "/"
				<< setw(2) << (int)ep[2] << " " << setw(2) << (int)ep[3] << ":" << setw(2) << (int)ep[4] << ":" << setw(2) << (int)ep[5] << " " << setfill(' ') << setw(8) << itpre->h << endl;
		}
	}
}
//组装矩阵
void Process::AssemblyDesignMatrix(Matrix& B, Matrix &x, Matrix &h, ObsSat_t &sat, ObsSat_t &RefObs, double* centerLL, opt_t opt, gtime_t &m_time, int row)
{
	int nmax = opt.nmax;
	int mmax = opt.mmax;
	int ntec = (nmax + 1)*(mmax + 1);
	double m_lat = centerLL[0];
	double m_lon = centerLL[1] + (sat.time.time - m_time.time)*PI / 43200;
	int GLO_F_K[24] = { 1, -4, 5, 6, 1, -4, 5, 6, -2, -7, 0, -1, -2, -7, 0, -1, 4, -3, 3, 2, 4, -3, 3, 2 };
	for (int i = 0;i < nmax + 1;++i)
	{
		for (int j = 0;j < mmax + 1;++j)
		{
			double t1, t2, t3, t4;
			t1 = (sat.lat - m_lat) / SCALE;
			t2 = (sat.lon - m_lon) / SCALE;
			t3 = (RefObs.lat - m_lat) / SCALE;
			t4 = (RefObs.lon - m_lon) / SCALE;
			if ((satsys(sat.prn) == SYS_GPS && satsys(RefObs.prn) == SYS_GPS) || (satsys(sat.prn) == SYS_GAL && satsys(RefObs.prn) == SYS_GAL))
			{
				B(row, i*(mmax + 1) + j) = (sat.mf*pow(t1, i)*pow(t2, j) - RefObs.mf*pow(t3, i)*pow(t4, j))*TEC2D / (FREQ1*FREQ1);
			}
			else if (satsys(sat.prn) == SYS_GLO && satsys(RefObs.prn) == SYS_GLO)
			{
				double satfreq = GLO_F_K[sat.prn - GPSMAXSAT - 1] * DFRQ1_GLO + FREQ1_GLO;
				double refreq = GLO_F_K[RefObs.prn - GPSMAXSAT - 1] * DFRQ1_GLO + FREQ1_GLO;
				B(row, i*(mmax + 1) + j) = sat.mf*pow(t1, i)*pow(t2, j)*TEC2D / (satfreq*satfreq) - RefObs.mf*pow(t3, i)*pow(t4, j)*TEC2D / (refreq*refreq);
			}
			else if (satsys(sat.prn) == SYS_BDS && satsys(RefObs.prn) == SYS_BDS)
			{
				B(row, i*(mmax + 1) + j) = (sat.mf*pow(t1, i)*pow(t2, j) - RefObs.mf*pow(t3, i)*pow(t4, j))*TEC2D / (FREQ1_BDS*FREQ1_BDS);
			}
			else
			{
				throw("AssemblyDesignMatrix error: 卫星号有误");
			}
		}
	}
	B(row, ntec + sat.prn - 1) = 1;
	B(row, ntec + RefObs.prn - 1) = -1;
	Matrix Bi(1, B.getCol());
	for (int i = 0;i<B.getCol();++i)
		Bi[i] = B(row, i);
	Matrix t = Bi*x;
	h[row] = t[0];


}
//法方程叠加
void Process::NorEquSup(pair<string, vector<ObsSat_t>> sta_onepoch, map<int, pair<int, double>> &ref_sat, double *centerLL, opt_t opt, gtime_t m_time, Matrix& HH, Matrix& Hy, Matrix &x, map<int, pair<int, double>> CnstrnSat)
{

#ifdef LOGFILE 
	if(inerflag!=0)
		fprintf(log, "%s:\n", sta_onepoch.first.c_str());
#endif
	
	string StaName = sta_onepoch.first;
	//sta_onepoch	--> <站名,数据>
	int GLO_F_K[24] = { 1, -4, 5, 6, 1, -4, 5, 6, -2, -7, 0, -1, -2, -7, 0, -1, 4, -3, 3, 2, 4, -3, 3, 2 };
	map<int, vector<ObsSat_t>>sys_onepoch;
	//sys_onepoch	--> map<sys,数据>
	int sys;
	for (vector<ObsSat_t>::iterator iter = sta_onepoch.second.begin();iter != sta_onepoch.second.end();++iter)
	{
		sys = satsys(iter->prn);
		if (!sys_onepoch.count(sys))
		{
			vector<ObsSat_t> temp;
			temp.push_back(*iter);
			sys_onepoch.insert(pair<int, vector<ObsSat_t>>(sys, temp));
		}
		else
		{
			map<int, vector<ObsSat_t>>::iterator temp;
			temp = sys_onepoch.find(sys);
			temp->second.push_back(*iter);
		}
	}
	//int nsat=SatNum(opt);
	int npara = (opt.nmax + 1)*(opt.mmax + 1) + MAXSAT;

	for (auto SingleObs : sys_onepoch)
	{	//SingleObs	--> <sys,数据>
		auto strsys = SYS.find(SingleObs.first);
	
#ifdef LOGFILE
		if(inerflag!=0)
			fprintf(log, "%s:\n", strsys->second.c_str());
#endif
		
		pair<int, double> ref;

		ref = ref_sat.find(SingleObs.first)->second;
		vector<ObsSat_t> single = SingleObs.second;//-->测站单个系统的数据
		ObsSat_t RefObs;
		for (vector<ObsSat_t>::iterator iter = single.begin();iter != single.end();++iter)
		{
			if (iter->prn == ref.first)
			{
				RefObs = *iter;//参考卫星数据
				break;
			}
		}


		int nobs = single.size();
		if (nobs <= 1)
			continue;
		Matrix B(nobs, npara);
		Matrix v(nobs, 1);
		Matrix W(nobs, nobs);
		Matrix h(nobs, 1);
		Matrix y(nobs, 1);
		vector<double> sig;
		int i = 0;
		vector<inerpre_t> pre;
		for (auto sat : single)
		{
			if (sat.prn == RefObs.prn)
				continue;
			if ((satsys(sat.prn) == SYS_GPS && satsys(RefObs.prn) == SYS_GPS) || (satsys(sat.prn) == SYS_GAL && satsys(RefObs.prn) == SYS_GAL))
				y[i] = (sat.stec - RefObs.stec)*TEC2D / (FREQ1*FREQ1);
			else if (satsys(sat.prn) == SYS_GLO && satsys(RefObs.prn) == SYS_GLO)
			{
				double satfreq = GLO_F_K[sat.prn - GPSMAXSAT - 1] * DFRQ1_GLO + FREQ1_GLO;
				double refreq = GLO_F_K[RefObs.prn - GPSMAXSAT - 1] * DFRQ1_GLO + FREQ1_GLO;
				y[i] = sat.stec*TEC2D / (satfreq*satfreq) - RefObs.stec*TEC2D / (refreq*refreq);
			}
			else if (satsys(sat.prn) == SYS_BDS && satsys(RefObs.prn) == SYS_BDS)
			{
				y[i] = (sat.stec - RefObs.stec)*TEC2D / (FREQ1_BDS*FREQ1_BDS);
			}
			else
			{
				throw("NorEquSup error: 卫星PRN有误");
			}
			AssemblyDesignMatrix(B, x, h, sat, RefObs, centerLL, opt, m_time, i);

			v[i] = y[i] - h[i];
			if (iterflag == 0)
				sig.push_back(1 / (2.5*2.5));
			else
			{
				sig.push_back(1 / (v[i] * v[i]));
			}
			if (inerflag != 0)
			{
				inerpre_t tmp;
				tmp.time = sat.time;
				tmp.y = y[i];
				tmp.h = h[i];
				tmp.v = v[i];
				tmp.sat = sat.prn;
				tmp.ref_sat = RefObs.prn;
				pre.push_back(tmp);
			}

			i++;
		};

		/*constraint sat dcb*/
		int cnstrnsat, ntec;
		cnstrnsat = CnstrnSat.find(SingleObs.first)->second.first;
		ntec = (opt.mmax + 1)*(opt.nmax + 1);
		B(nobs - 1, ntec + cnstrnsat - 1) = 1;
		if (iterflag == 0)
			v[nobs - 1] = 1;
		else
			v[nobs - 1] = 0;
		sig.push_back(1E8);
		CreateDiag(W, sig);

		if (inerflag != 0)
		{
#ifdef LOGFILE 
			fprintf(log, "constraint sat: %02d\n", cnstrnsat);
			fprintf(log, "reference  sat: %02d\n", RefObs.prn);
			fprintf(log, " B Matrix:\n");
			B.fileprint(log);
			//fprintf(log," x Matrix:\n");
			//x.fileprint(log);
			fprintf(log, " W Matrix:\n");
			W.fileprint(log);
			fprintf(log, " v Matrix:\n");
			v.fileprint(log);
#endif
			OutputInerPre(StaName, pre);
		}
		else
		{
			Hy = Hy + B.trans()*W*v;
			HH = HH + B.trans()*W*B;
		}
	}
}
//计算卫星个数
extern int SatNum(opt_t opt)
{
	int nsat = 0;
	switch (opt.sys)
	{
	case SYS_GPS:
		nsat = GPSMAXSAT;
		break;
	case SYS_GLO:
		nsat = GLOMAXSAT;
		break;
	case SYS_GAL:
		nsat = GALMAXSAT;
		break;
	case SYS_BDS:
		nsat = BDSMAXSAT;
		break;
	case (SYS_GPS | SYS_GLO) :
		nsat = GPSMAXSAT + GLOMAXSAT;
		break;
	case (SYS_GPS | SYS_GAL) :
		nsat = GPSMAXSAT + GALMAXSAT;
		break;
	case (SYS_GPS | SYS_BDS) :
		nsat = GPSMAXSAT + BDSMAXSAT;
		break;
	case (SYS_GPS | SYS_GLO | SYS_GAL) :
		nsat = GPSMAXSAT + GLOMAXSAT + GALMAXSAT;
		break;
	case (SYS_GPS | SYS_GLO | SYS_GAL | SYS_BDS) :
		nsat = MAXSAT;
		break;
	}
	return nsat;
}
//建模核心处理函数
void Process::ProcessModel(map<gtime_t, map<string, vector<ObsSat_t>>> &ObsSort, double *centerLL, opt_t opt, gtime_t m_time) {
	int iter_time = 1;
	double ep[6];
	double lastmax = 99;
	time2epoch(m_time, ep);
	printf("\n%4d %02d %02d %02d:%02d:%02d Modeling .....\n", (int)ep[0], (int)ep[1], (int)ep[2], (int)ep[3], (int)ep[4], (int)ep[5]);

	iterflag = 0;//是否首次迭代flag
	cnstflag = 0;//卫星约束flag
	inerflag = 0;//内符合精度flag

	int npara = (opt.nmax + 1)*(opt.mmax + 1) + MAXSAT;
	Matrix HH(npara, npara);
	Matrix HY(npara, 1);
	Matrix x(npara, 1);

#ifdef LOGFILE
		fprintf(log, "建模时刻：%4d %02d %02d %02d:%02d:%02d\n", (int)ep[0], (int)ep[1], (int)ep[2], (int)ep[3], (int)ep[4], (int)ep[5]);
#endif
	
	map<int, pair<int, double>> CnstrnSat;
	while (1)
	{
		if (inerflag != 0)
		{
#ifdef LOGFILE 
			//fprintf(log, "\n第%02d次迭代\n", iter_time);
#endif
		}
		for (auto onepoch : ObsSort)
		{//onepoch	-->取一个历元所有数据
#ifdef LOGFILE
			double epp[6];
			time2epoch(onepoch.first, epp);
			if(inerflag!=0)
				fprintf(log, "\n历元时刻：%4d %02d %02d %02d:%02d:%02d\n", (int)epp[0], (int)epp[1], (int)epp[2], (int)epp[3], (int)epp[4], (int)epp[5]);
#endif
			for (auto sta_onepoch : onepoch.second)
			{//sta_onepoch	-->取一个历元一个站的数据
				map<int, pair<int, double>> ref_sat;
				ref_sat = SelRefSat(sta_onepoch);
				if (cnstflag == 0)
					CnstrnSat = ref_sat;
				NorEquSup(sta_onepoch, ref_sat, centerLL, opt, m_time, HH, HY, x, CnstrnSat);
				cnstflag = 1;
			}
		}
		if (inerflag != 0)
		{
			std::cout << "\nComplished,going to the next time ......\n";
			break;
		}
		printf("第%d次迭代\n", iter_time);

		iterflag++;
		iter_time++;

		set<int> erow, ecol, rc;
		Matrix HB, Hy, Hh;
		erow = EraseZerosRow(HH, HB);
		ecol = EraseZerosCol(HB, Hh);
		rc = EraseZerosRow(HY, Hy);
		if (rc.size() != erow.size() || erow.size() != ecol.size())
			throw("Matrix error: HH、Hy矩阵不正确！");
		Matrix dx(Hy.getRow(), 1);
		dx = Hh.inv()*Hy;
	
		//x.print();
		//WashX(x,opt);

		double max = MinValue;
		for (int i = 0;i<(opt.nmax + 1)*(opt.mmax + 1);++i)
		{
			if (max<fabs(dx[i]))
				max = fabs(dx[i]);
		}
		if (max == MinValue)
			throw("ProcessModel error: max=MinValue ");

		printf("max= %.5f\n", max);
		if (max <= THRESHOLD||max>=lastmax)
		{
			OutputFunc(m_time, x, centerLL, opt);
			inerflag++;
		}
		if (inerflag == 0)
		{
			int j = 0;
			for (int i = 0;i < x.getRow();++i)
			{
				if (erow.count(i))
					continue;
				x[i] += dx[j];
				++j;
			}
		}
		lastmax = max;
	}

}
void Process::FindFirEpoch(string &str)
{
	string::size_type size,pos,time_pos;
	string s;
	gtime_t stime = { 0 }, etime={ 0 },time;

	time_pos = 0;
	size = str.size();
	pos = str.find("\n", 0);
	int inv = pos + 1;
	s = str.substr(0, inv);
	stime.time = atoi(s.substr(0, 11).c_str());
	s = str.substr(str.size() - inv, string::npos);
	etime.time = atoi(s.substr(0, 11).c_str());
	
	for (time = stime;time <= etime;time = time+ opt.ObsInv)
	{
		if ((time >= opt.StartTime) && (time <= opt.EndTime))
		{
			FirstEpoch = time;
			if (etime >= opt.EndTime)
			{
				int i = (int)(opt.EndTime.time-time.time)/opt.ObsInv;
				LastEpoch = FirstEpoch + i*opt.ObsInv;
			}
			break;
		}
	}
	
}