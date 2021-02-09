// bool create (pii x, pii y) {
// 	x.F = max (0, x.F); x.S = min (n - 1, x.S);
// 	y.F = max (0, y.F); y.S = min (m - 1, y.S);
// 	if (x.F > x.S || y.F > y.S)	return false;
// 	rid ++;
// 	R rt;
// 	for (int i = x.F; i <= x.S; i ++)
// 		for (int j = y.F; j <= y.S; j ++) {
// 			grid[i][j].rt = rid;
// 			rt.x.F = min (rt.x.F, grid[i][j].x.F), rt.x.S = max (rt.x.S, grid[i][j].x.S);
// 			rt.y.F = min (rt.y.F, grid[i][j].y.F), rt.y.S = max (rt.y.S, grid[i][j].y.S);
// 		}
// 	rt.bounded();
// 	pts.push_back(rt);
// 	return true;
// }

// void divide (pii x, pii y, int l) {
// 	if (l <= 9) {
// 		bool ret = false;
// 		for (int dx = -1; dx <= 1; dx ++)
// 			for (int dy = -1; dy <= 1; dy ++)
// 				if (dx || dy)
// 					ret |= create(MP(x.F + dx * l, x.S + dx * l), MP(y.F + dy * l, y.S + dy * l));
// 		if (!ret)
// 			return ;
// 		divide (MP(x.F - l, x.S + l), MP(y.F - l, y.S + l), l * 3);
// 	}
// 	else {
// 		pii ax = x;
// 		while (ax.S > 0) {
// 			ax.F -= l, ax.S -= l; create (ax, y);
// 			pii ay = y;	while (ay.S > 0) { ay.F -= l, ay.S -= l; create (ax, ay);}
// 			ay = y;		while (ay.F < m) { ay.F += l, ay.S += l; create (ax, ay);}
// 		}
// 		ax = x;
// 		pii ay = y;	while (ay.S > 0) { ay.F -= l, ay.S -= l; create (x, ay);}
// 		ay = y;		while (ay.F < m) { ay.F += l, ay.S += l; create (x, ay);}
// 		while (ax.F < n) {
// 			ax.F += l, ax.S += l; create (ax, y);
// 			pii ay = y;	while (ay.S > 0) { ay.F -= l, ay.S -= l; create (ax, ay);}
// 			ay = y;		while (ay.F < m) { ay.F += l, ay.S += l; create (ax, ay);}
// 		}
// 	}
// }

// void decide_epic_ct () {
// 	epic_ct.F = 48.994 * Sla;
// 	epic_ct.S = -122.702 * Sln;
// 	// freopen("positive_position.txt", "r", stdin);
// 	// pll pt, sum = MP(0.0, 0.0); int num = 0;
// 	// while (cin >> pt.F >> pt.S)
// 	// 	sum = MP(sum.F + pt.F, sum.S + pt.S), num ++;
// 	// epic_ct = MP(sum.F / num, sum.S / num);
// }