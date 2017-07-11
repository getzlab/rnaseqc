#ifndef INTERVAL_H
#define INTERVAL_H

#include "GTF.h"
#include <vector>
#include <utility>
#include <set>

//struct Interval {
//	coord start;
//	coord end;
//};

typedef Feature Interval;

typedef std::set<Interval, bool(*) (const Interval&, const Interval&)> resultSet;
//typedef std::unordered_set<Interval> resultSet;

bool compIntervalStart(const Interval&, const Interval&);
bool compIntervalEnd(const Interval&, const Interval&);
bool intersectPoint(const Interval&, const coord);
bool intersectInterval(const Interval&, const Interval&);
bool operator==(const Interval&, const Interval&);

class IntervalTree {
	std::vector<Interval> starts;
	std::vector<Interval> ends;
	coord center;
	IntervalTree* left;
	IntervalTree* right;
	std::vector<std::pair<coord, Interval*> > points;
    
    void walkPoints(resultSet*, Interval&);

public:
    IntervalTree(std::vector<Interval>&);
    //IntervalTree(std::vector<Interval> &intervals) : IntervalTree(intervals, true) {};
	~IntervalTree();
	std::vector<Interval>* queryPoint(coord);
    std::vector<Interval>* queryPoint(coord, bool);
    resultSet* queryInterval(Interval&);
};

IntervalTree* constructIntervalTree(std::vector<Interval>&);
#endif //INTERVAL_H
