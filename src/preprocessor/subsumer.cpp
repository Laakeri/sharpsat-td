#include "subsumer.hpp"

#include <vector>
#include <algorithm>
#include <iostream>

#include "utils.hpp"

#define getD(e,i) data[(e.B) + (i)]

namespace sspp {

void Subsumer::assumeSize(size_t size) {
	if (data.size() < size) {
		data.resize(size);
	}
}

bool Subsumer::isPrefix(const vecP a, const vecP b) {
	for (size_t i = 0; i < a.size(); i++) {
		if (getD(a, i) != getD(b, i)) return false;
	}
	return true;
}

bool Subsumer::CSO1(const vector<vecP>& D, size_t b, size_t e, const vecP S, 
					   size_t j, size_t d, const vector<int>& d0Index) {
	while (j < S.size() && getD(S, j) < getD(D[b], d)) j++;
	if (j >= S.size()) return false;
	
	if (getD(S, j) == getD(D[b], d)) {
		size_t ee = b;
		if (d == 0) {
			ee = d0Index[getD(S, j) + 1] - 1;
		}
		else if (e - b < BSconst) {
			while (ee + 1 <= e && getD(D[ee + 1], d) == getD(S, j)) ee++;
		}
		else {
			int mi = b;
			int ma = e - 1;
			while (mi <= ma) {
				int mid = (mi+ma)/2;
				if (getD(D[mid + 1], d) == getD(S, j)) {
					ee = mid + 1;
					mi = mid + 1;
				}
				else {
					ma = mid - 1;
				}
			}
		}
		if (S.size() > d + 1 && D[b].size() == d + 1) {
			return true;
		}
		if (j + 1 <= S.size()) {
			if (CSO1(D, b, ee, S, j + 1, d + 1, d0Index)) {
				return true;
			}
		}
		b = ee + 1;
	}
	else {
		if (d == 0) {
			b = d0Index[getD(S, j)];
		}
		else if (e - b < BSconst) {
			while (b <= e && getD(D[b], d) < getD(S, j)) b++;
		}
		else {
			int mi = b;
			int ma = e;
			while (mi <= ma) {
				int mid = (mi+ma)/2;
				if (getD(D[mid], d) < getD(S, j)) {
					b = mid + 1;
					mi = mid + 1;
				}
				else {
					ma = mid - 1;
				}
			}
		}
	}
	if (b <= e) return CSO1(D, b, e, S, j, d, d0Index);
	return false;
}

vector<vector<Lit>> Subsumer::Subsume(const vector<vector<Lit>>& clauses) {
	vector<vecP> D(clauses.size());
	vector<int> d0Index;
	ALIt++;
	int pid = 0;
	int maxFreq = 0;
	size_t dataP = 0;
	for (size_t i = 0; i < clauses.size(); i++) {
		if (clauses[i].empty()) {
			return {{}};
		}
		for (size_t j = 0; j < clauses[i].size(); j++) {
			data.push_back(clauses[i][j]);
		}
		D[i].B = dataP;
		D[i].E = dataP + clauses[i].size();
		D[i].O = i;
		dataP += D[i].size();
		
		for (size_t j = D[i].B; j < D[i].E; j++) {
			if ((int)ALU.size() <= data[j]) {
				ALU.resize(data[j] + 1);
				ALI.resize(data[j] + 1);
			}
			if (ALU[data[j]] != ALIt) {
				ALU[data[j]] = ALIt;
				if (pid >= (int)ALF.size()) ALF.resize(pid + 1);
				ALF[pid] = 0;
				ALI[data[j]] = pid++;
			}
			data[j] = ALI[data[j]];
			ALF[data[j]]++;
			maxFreq = max(maxFreq, ALF[data[j]]);
		}
	}
	vector<vector<int> > cSort(maxFreq + 1);
	for (int i = 0; i < pid; i++) {
		cSort[ALF[i]].push_back(i);
	}
	vector<size_t> perm(pid);
	size_t i2 = 0;
	for (int i = 0; i <= maxFreq; i++) {
		for (int t : cSort[i]) {
			perm[t] = i2++;
		}
	}
	for (size_t i = 0; i < D.size(); i++) {
		for (size_t j = D[i].B; j < D[i].E; j++) {
			data[j] = perm[data[j]];
		}
		std::sort(data.data() + D[i].B, data.data() + D[i].E);
	}
	auto cmp = [&](vecP a, vecP b) {
		if (getD(a, 0) < getD(b, 0)) return true;
		else if (getD(a, 0) > getD(b, 0)) return false;
		size_t sz = min(a.size(), b.size());
		for (size_t i = 1; i < sz; i++) {
			if (getD(a, i) < getD(b, i)) return true;
			else if (getD(a, i) > getD(b, i)) return false;
		}
		return a.size() < b.size();
	};
	vector<vector<vecP> > cSort2(pid);
	for (size_t i = 0; i < D.size(); i++) {
		cSort2[getD(D[i], 0)].push_back(D[i]);
	}
	size_t di = 0;
	for (int i = 0; i < pid; i++) {
		sort(cSort2[i].begin(), cSort2[i].end(), cmp);
		for (size_t j = 0; j < cSort2[i].size(); j++) {
			D[di++] = cSort2[i][j];
		}
	}
	d0Index.resize(pid + 1);
	int idx = 0;
	for (size_t i = 0; i < D.size(); i++) {
		if (i == 0 || getD(D[i], 0) > getD(D[i - 1], 0)) {
			while (idx <= getD(D[i], 0)) {
				d0Index[idx++] = i;
			}
		}
	}
	while (idx <= pid) {
		d0Index[idx++] = (int)D.size();
	}
	vector<int> subsumed(D.size());
	size_t S = 0;
	for (size_t i = 1; i < D.size(); i++) {
		if (D[S].size() <= D[i].size() && isPrefix(D[S], D[i])) {
			subsumed[i] = true;
		}
		else {
			S = i;
		}
	}
	for (size_t i = 0; i + 1 < D.size(); i++) {
		if (!subsumed[i] && CSO1(D, i + 1, D.size() - 1, D[i], 0, 0, d0Index)) {
			subsumed[i] = true;
		}
	}
	vector<vector<Lit>> ret;
	for (size_t i = 0; i < D.size(); i++) {
		if (!subsumed[i]) {
			ret.push_back(clauses[D[i].O]);
		}
	}
	return ret;
}
} // namespace sspp