/*
*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
* 
* @author: Xingjian Bai 
* @date: 2021-02-07 13:38:16
* @description: 
*  
* 
* @notes: 
* g++ -fsanitize=address -ftrapv 
*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  */
#include <bits/stdc++.h>
#define F first
#define S second
#define MP make_pair
#define TIME (double)clock()/CLOCKS_PER_SEC
using namespace std;
typedef long long ll;
typedef long double ld;
typedef pair <int, int> pii;
const int mod = 1000000007; 
const ll oo = 1e18;
const ld eps = 1e-8;
#define debug(x) cerr << "(debug mod) " << #x << " = " << x << endl

int mn[13] = {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
int bar[1010];

int main() {
    ios::sync_with_stdio(false);
    freopen("time2.txt", "r", stdin);
    freopen("time4.txt", "w", stdout);
    int m, d;
    while (cin >> m >> d) {
    	int ans = 0;
    	// for (int i = 1; i < m; i ++)
    		// ans += mn[i];
    	bar[m] ++;
    }
  //   for (int i = 1; i <= 365; i += 7) {
		// int sum = 0;
  //   	for (int j = i; j < i + 7; j ++)
  //   		sum += bar[j];
  //   	cout << sum << endl;
  //   }
    for (int i = 1; i <= 12; i ++)
        cout << bar[i] << endl;
    return 0;
}