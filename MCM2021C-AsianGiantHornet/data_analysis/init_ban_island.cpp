/*
*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
* 
* @author: Xingjian Bai 
* @date: 2021-02-05 18:20:28
* @description: 
*  
* 
* @notes: 
* g++ -O2 -fsanitize=address -ftrapv 
*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  */
#include <bits/stdc++.h>
#define F first
#define S second
#define MP make_pair
#define TIME (double)clock()/CLOCKS_PER_SEC
#define PI 3.14159265358980
using namespace std;
typedef long long ll;
typedef long double ld;
typedef pair <int, int> pii;
typedef pair <ld, ld> pll;
const int mod = 1000000007; 
const ll oo = 1e18;
const ld eps = 1e-8;
#define debug(x) cerr << "(debug mod) " << #x << " = " << x << endl

const ld Sla = 364000, Sln = 288200;
const pll la = MP(45.088689 * Sla, 49.948004 * Sla);	//MP(45.488689 * Sla, 49.548004 * Sla);
const pll ln = MP(-125.165014 * Sln, -116.373687 * Sln);//ln = MP(-124.665014 * Sln, -116.873687 * Sln);
const ld l = 32808.4 * 1.2; //l feet = 12 km

pll epic_ct;
pii epic_id;
vector<pair<pll, ld> > all_sut;
vector<pll> all_data;

ld M(ld x, ld y) {return (x + y) * 0.5;}
ld D(pll a) {return fabs(a.F - a.S);}
ld Rand() {return (ld)rand()/RAND_MAX;}
bool within (pll r, ld x) {return r.F <= x + eps && x <= r.S + eps;}
ld Gaussian (pll p) {return exp(-(p.F*p.F+p.S*p.S)*8)/sqrt(2*PI);}
bool if_winter (int x) {return x % 365 <= 183;}

struct R {
	//each region
	pll x, y;
	int rt; // the larger region it belongs to
	ld popu, popi, sut; int sut_num; //population, population index, sutability, index

	R () {x = y = MP(oo, -oo); rt = -1;}
	R (pll _x, pll _y) {x = _x, y = _y; rt = -1;}

	pll center()		{return MP(M(x.F, x.S), M(y.F, y.S));}
	ld area ()			{return D(x)*D(y);}
	bool if_in(pll p)	{return within (x, p.F) && within(y, p.S);}
	void count_population () {
		for (int i = 0; i < all_data.size(); i ++)
			popu += Gaussian(MP(all_data[i].F - center().F / Sla, all_data[i].S - center().S / Sln));
		popi = -log(popu);
	}
	void popu_justify() {if (sut == 0) popu = oo, popi = -100;}
	void count_sutability () {
		sut_num = sut = 0;
		for (int i = 0; i < all_sut.size(); i ++) {
			pll p = all_sut[i].F;
			ld v = all_sut[i].S;
			ld var = l/2;
			if (within(x, p.F) && within(y, p.S))
				sut += v, sut_num += 1;
			// else if ((within(x, p.F - var) || within(x, p.F + var)) && (within(y, p.S - var) || within(y, p.S + var)))
			// 	sut += v*0.9, sut_num += 1;
			// else if ((within(x, p.F - var*2) || within(x, p.F) || within(x, p.F + var*2)) && (within(y, p.S - var*2) || within(y, p.S) || within(y, p.S + var*2)))
			// 	sut += v*0.8, sut_num += 1;
		}
		if (sut_num == 0)	sut = 0;
		else				sut /= sut_num;
		// cerr << sut << " " << sut_num << endl;
	}
	ld popi_to_sut()	{return (popi + 8) / 30.0;}
	void output ()		{cout << "grid: [" << x.F  / Sla << ", " << x.S  / Sla << "], [" << y.F / Sln << ", " << y.S / Sln << "], s = " << x.S / Sla - x.F / Sla << "\n";}
	void output_ct ()	{cout << center().F / Sla << ", " << center().S / Sln << '\n';}
	void output_ct_id (){cout << center().F / Sla << ", " << center().S / Sln << ", " << sut << '\n';}
	void output_pop ()	{cout << center().F / Sla << ", " << center().S / Sln << ", " << popu << ", " << sut << '\n';}
	void output_popi ()	{cout << center().F / Sla << ", " << center().S / Sln << ", " << popu << ", " << popi << '\n';}
	void output_sut ()	{cout << center().F / Sla << ", " << center().S / Sln << ", " << (sut>eps) << '\n';}
};

void decide_epic_ct () {
	epic_ct.F = 48.994 * Sla;
	epic_ct.S = -122.702 * Sln;
}

const int MAXN = 100;
int n, m;
vector <vector <R> > grid;


pii divide_to_grids(ld l) {
	pii epic_id;
	n = ceil(D(la) / l);
	m = ceil(D(ln) / l);
	ld lx, ly;	int i, j;
	for (i = 0, lx = la.F; i < n; i ++, lx += l) {
		vector <R> row(m);
		for (j = 0, ly = ln.F; j < m; j ++, ly += l) {
			row[j] = R(MP(lx, lx + l), MP(ly, ly + l));
			row[j].rt = i * m + j;
			if (row[j].if_in(epic_ct))
				epic_id = MP(i, j);
		}
		grid.push_back(row);
	}
	return epic_id;
}

void output_grid () {
	cerr << "tot grids: " << n << " * " << m << " = " << n * m << endl;
	ld mx = -oo, mn = oo;
	for (int i = 0; i < grid.size(); i ++)
		for (int j = 0; j < grid[i].size(); j ++) {
			grid[i][j].output_sut();
			mx = max (mx, grid[i][j].popu);
			mn = min (mn, grid[i][j].popu);
			// grid[i][j].output_ct();
		}
	cerr << mx << " -- " << mn << endl;
}

// ld adjmp[3200][3200];
// void output_adjmap() {

// 	for (int i = 0; i < n; i ++)
// 		for (int j = 0; j < m; j ++) {
// 			int cnt = (i != 0) + (j != 0) + (i < n - 1) + (j < m - 1);
// 			if (i)	adjmp[i*m+j][i*m+j-m] = 1.0 / cnt;
// 			if (j)	adjmp[i*m+j][i*m+j-1] = 1.0 / cnt;
// 			if (j<m-1)	adjmp[i*m+j][i*m+j+1] = 1.0 / cnt;
// 			if (i<n-1)	adjmp[i*m+j][i*m+j+m] = 1.0 / cnt;
// 		}
// 	for (int i = 0; i < n * m; i ++) {
// 		for (int j = 0; j < n * m; j ++) {
// 			cout << adjmp[i][j];
// 			if (j != n * m - 1)
// 				cout << ',';
// 		}
// 		cout << '\n';
// 	}
// }

ld hab[MAXN][MAXN];
#define MAXL 1000 //the starting number of bees
#define MINL 0
struct State {
	ld lev[MAXN][MAXN];
	State () {memset(lev, 0, sizeof lev);}
	int tid;
	ld tot () {ld sum = 0;for (int i = 0; i < n; i ++)for (int j = 0; j < m; j ++)sum += lev[i][j];return sum;}
	ld av () {return tot() / n / m;}
	int non_emp () {int sum = 0;for (int i = 0; i < n; i ++)for (int j = 0; j < m; j ++)sum += (lev[i][j] > eps);return sum;}
	ld mx () {ld mx_v = 0;for (int i = 0; i < n; i ++)for (int j = 0; j < m; j ++)mx_v = max(lev[i][j], mx_v);return mx_v;}
	ld top10av () {
		priority_queue <ld> q; ld sum = 0;
		for (int i = 0; i < n; i ++) for (int j = 0; j < m; j ++) q.push(lev[i][j]);
		for (int i = 0; i < 10; i ++) {sum += q.top(); q.pop();}
		return sum / 10;
	}
	void init () {
		memset(lev, 0, sizeof lev);
		lev[epic_id.F][epic_id.S] = MAXL;
	}
	void output () {
		cerr << "data at t = " << tid << ": tot = " << tot() << ", non-emp = " << non_emp() << ", mx = " << mx() << ", top10av = " << top10av() << "\n";
		// for (int i = 0; i < n; i ++)
		// 	for (int j = 0; j < m; j ++)
		// 		cout << grid[i][j].center().F << ", " << grid[i][j].center().S << ", " << lev[i][j] << '\n';
	}
	void output2 () {
		cerr << tid << ", " << tot() << '\n';
		// for (int i = 0; i < n; i ++)
		// 	for (int j = 0; j < m; j ++)
		// 		cout << grid[i][j].center().F << ", " << grid[i][j].center().S << ", " << lev[i][j] << '\n';
	}
	void print(string info) {
		stringstream ss;
		ss << "frames/frame" << tid << ".csv";
		// ss << "frames/frame" << tid << "-" << tid/365 << "-" << info << ".csv";
		freopen(ss.str().c_str(), "w", stdout);
		cout << "latitude, longitude, count\n";
		for (int i = 0; i < n; i ++)
			for (int j = 0; j < m; j ++)
				cout << grid[i][j].center().F / Sla << ", " << grid[i][j].center().S / Sln << ", " << lev[i][j] << endl;
		fclose(stdout);
	}
	void print2() {
		stringstream ss;
		ss << "frames/frame_mat" << tid << ".csv";
		freopen(ss.str().c_str(), "w", stdout);
		for (int i = 0; i < n; i ++) {
			for (int j = 0; j < m; j ++)
				cout << lev[i][j] << ' ';
			cout << endl;
		}
		fclose(stdout);
	}
} stt[365 * 100];

//popi
const ld EQ = 2;
ld calc_prob_pop (int px, int py, int cx, int cy) {
	if (cx < 0 || cx >= n || cy >= m || cy < 0) //cy < 0 || 
		return EQ;
	// if (cy < 0)
	// 	return 0;
	// cerr << grid[cx][cy].popi << " " << grid[px][py].popi << endl;
	return exp(grid[cx][cy].popi - grid[px][py].popi) + EQ;
}

void propagate (State pre, bool if_winter) {
	//propagate State(t) to State(t + 1)
	ld p_rem = (if_winter ? 0.9 : 0.7);
	for (int i = 0; i < n; i ++)
		for (int j = 0; j < m; j ++) {
			ld rep = (if_winter ? 1 : 1 + 0.02 * grid[i][j].popi_to_sut());
			// if (Rand() < p_rem) {
			// 	stt[pre.tid + 1].lev[i][j] += pre.lev[i][j];
			// 	continue ;
			// }
			pre.lev[i][j] *= rep;
			stt[pre.tid + 1].lev[i][j] += pre.lev[i][j] * p_rem;
			pre.lev[i][j] *= (1 - p_rem);

			//equal chance to escape
			ld flee[5];
			flee[0] = calc_prob_pop (i, j, i - 1, j);
			flee[1] = calc_prob_pop (i, j, i, j - 1);
			flee[2] = calc_prob_pop (i, j, i + 1, j);
			flee[3] = calc_prob_pop (i, j, i, j + 1);
			flee[4] = calc_prob_pop (i, j, i, j);
			ld exps = flee[0] + flee[1] + flee[2] + flee[3] + flee[4];
			// if (Rand() > 0.999)
				// cerr << "pp1: " << flee[0] << " " << flee[1] << " " <<  flee[2] << " " <<  flee[3] << " " << flee[4] << endl;
			if (i > 0)		stt[pre.tid + 1].lev[i - 1][j] += pre.lev[i][j] * flee[0] / exps;
			if (j > 0)		stt[pre.tid + 1].lev[i][j - 1] += pre.lev[i][j] * flee[1] / exps;
			if (i < n - 1)	stt[pre.tid + 1].lev[i + 1][j] += pre.lev[i][j] * flee[2] / exps;
			if (j < m - 1)	stt[pre.tid + 1].lev[i][j + 1] += pre.lev[i][j] * flee[3] / exps;
							stt[pre.tid + 1].lev[i][j] += pre.lev[i][j] * flee[4] / exps;
		}
	stt[pre.tid + 1].tid = pre.tid + 1;
}


void count_popu () {
	for (int i = 0; i < n; i ++)
		for (int j = 0; j < m; j ++)
			grid[i][j].count_population();
}
void sut_to_popu() {
	for (int i = 0; i < n; i ++)
		for (int j = 0; j < m; j ++)
			grid[i][j].popu_justify();
}
void count_sut () {
	for (int i = 0; i < n; i ++)
		for (int j = 0; j < m; j ++)
			grid[i][j].count_sutability();
	for (int i = 0; i < n; i ++)
		for (int j = 0; j < m; j ++)
			if (grid[i][j].sut == 0 && i && i < n - 1 && grid[i - 1][j].sut && grid[i + 1][j].sut)
				grid[i][j].sut = 1, grid[i][j].sut_num = 1;
}

int months[13] = {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
ld coef[4] = {7.07317231, 0.53359878, -0.18054375, 0.04970395};
struct Rep {
	int date, type;
	pll loc; pii pos;
	ld q_time;
	void find_pos() {
		for (int i = 0; i < n; i ++)
			for (int j = 0; j < m; j ++)
				if(grid[i][j].if_in(loc))
					pos = MP(i, j);
	}
	void input (int m, int d, string tp, ld lx, ld ly) {
		date = d;
		for (int i = 1; i < m; i ++)
			date += months[i];
		q_time = max((ld)0, coef[2] * sin(coef[1] * (date / 365.0 * 12 + 1) + coef[0]) + coef[3]);
		type = (tp == "Unverified" || tp == "Unprocessed" ? 0 : (tp == "Negative" ? -1 : 1));
		loc = MP(lx * Sla, ly * Sln);
		find_pos();
	}
	void output () {
		cout << pos.F << " " << pos.S << " " << date << " " << type << " " << q_time << endl;
	}
}; vector <Rep> reps;

//sorting of requests
int time_stamp = 0;
bool ord (const Rep &a, const Rep &b) {
	return a.q_time * stt[time_stamp].lev[a.pos.F][a.pos.S] > b.q_time * stt[time_stamp].lev[b.pos.F][b.pos.S];
}

void input_all_data () {
	freopen("all_data.csv", "r", stdin);
	ld x, y, z;
	while (cin >> x >> y) {
		if (x == -1)
			break;
		all_data.push_back(MP(x, y));
	}
	while (cin >> x >> y >> z) {
		if (x == -1)
			break;
		if (within(la, y * Sla) && within(ln, x * Sln))
			all_sut.push_back(MP(MP(y * Sla, x * Sln), z));
	}
	ld lx, ly; string tp; int m, d;
	while (cin >> m >> d >> tp >> lx >> ly) {
		Rep ne; ne.input(m, d, tp, lx, ly);
		reps.push_back(ne);
		// if (ne.type == 1) {
		// 	cerr << "find pos: " << m << " " << d << " " << tp << " " << lx << " " << ly << endl;
		// 	ne.output();
		// }
	}
	fclose(stdin);
}

void evolve() {
	stt[0].init();
    for (int T = 0; T <= 365*12; T ++) {
    	// cerr << " : " << T << endl;
    	propagate(stt[T], if_winter(T));
    	if (T % 10 == 0)
	    	stt[T].output();
	    if ((T <= 300 && T % 50 == 0) || (T % 200 == 0)) {
    		// cerr << (T % 365 ? "summer" : "winter") << " comes! " << endl;
    		stt[T].print(T % 365 ? " " : " ");
    		// stt[T].print2();
    	}
    }
}

int main() {
	srand(time(0));
	cout << fixed << setprecision(5);
	cerr << fixed << setprecision(5);
    ios::sync_with_stdio(false);
    decide_epic_ct();
    
    epic_id = divide_to_grids(l);
    input_all_data();
    count_sut ();
    count_popu ();
    sut_to_popu();
    // output_grid();
    // output_adjmap();


    evolve();
    return 0;
}