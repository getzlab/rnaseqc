#include "Interval.h"

#include <algorithm>
#include <iterator>
using std::sort;
using std::copy;
using std::back_insert_iterator;
using std::vector;
using std::pair;
using std::make_pair;
using std::inserter;


bool compIntervalStart(const Interval &a, const Interval &b)
{
	return a.start < b.start;
}

bool compIntervalEnd(const Interval &a, const Interval &b)
{
	return a.end < b.end;
}

bool intersectPoint(const Interval &a, const coord x)
{
	return (x >= a.start) && (x <= a.end);
}

bool intersectInterval(const Interval &a, const Interval &b)
{
	return intersectPoint(a, b.start) || intersectPoint(a, b.end) || intersectPoint(b, a.start);
}

//bool operator==(const Interval &a, const Interval &b)
//{
//	return (a.start == b.start) && (a.end == b.end);
//}

IntervalTree::IntervalTree(vector<Interval> &intervals) : starts(vector<Interval>()), ends(vector<Interval>()), center(0ll), left(nullptr), right(nullptr)
{
	//pick a midpoint that will (roughly) balance the left and right sides
	this->center = (intervals.begin()->start + intervals.rbegin()->end)/2;
	vector<Interval> left, right;
	this->starts = vector<Interval>();
	auto curr = intervals.begin();
	auto end = intervals.end();
	//iterate over the intervals, sorting them into the left, center, and right bins
	while (curr != end)
	{
        
		if (intersectPoint(*curr, this->center)) this->starts.push_back(*curr);
		else if (curr->end < this->center) left.push_back(*curr);
		else right.push_back(*curr);
		++curr;
	}
    intervals.clear();
	if (left.size()) //if any intervals were put in the left bin, generate a tree on the left
	{
		this->left = new IntervalTree(left);
		left.clear();
	}
	if (right.size()) //same thing for the right bin
	{
		this->right = new IntervalTree(right);
		right.clear();
	}
	//now copy the contents from the starts bin (which were already sorted by start point) into the ends list, and sort by ends
	this->ends = vector<Interval>();
	copy(this->starts.begin(), this->starts.end(), back_insert_iterator<vector<Interval> >(this->ends));
	sort(this->ends.begin(), this->ends.end(), compIntervalEnd);

	//now create a list of <start/end point, interval> pairs
	//Used for querying intersections against intervals
	auto head = this->starts.begin();
	auto tail = this->starts.end();
	this->points = vector<pair<coord, Interval*> >();
	while (head != tail)
	{
		//for each interval in this tree, add a pair for it's start and end points
		this->points.push_back(make_pair(head->start, &(*head)));
		this->points.push_back(make_pair(head->end, &(*head)));
		++head;
	}
	sort(this->points.begin(), this->points.end()); //by default, it will sort by the first element

}

IntervalTree::~IntervalTree()
{
	delete this->left;
	this->left = nullptr;
	delete this->right;
	this->right = nullptr;
}

vector<Interval>* IntervalTree::queryPoint(coord x)
{
    return this->queryPoint(x, true);
}

vector<Interval>* IntervalTree::queryPoint(coord x, bool recursive)
{
	vector<Interval>* output = new vector<Interval>();
	if (x <= this->center)
	{ //if the point is left of center, query the left tree, and add intervals from this tree intersecting the point
		if (this->left != nullptr && recursive)
		{
			//copy in results from the left tree iff there is a left tree
			vector<Interval>* left = this->left->queryPoint(x);
			copy(left->begin(), left->end(), back_insert_iterator<vector<Interval> >(*output));
			delete left;
			left = nullptr;
		}
		auto curr = this->starts.begin();
		auto end = this->starts.end();
		while (curr != end && curr->start <= x)
		{
			//since this point is left of the center and any interval in this tree intersects the center
			//then, as long as an interval starts left of x, we know it intersects
			output->push_back(*curr);
			++curr;
		}
	}
	else
	{
		//same idea as above, but applied to interval end points
		if (this->right != nullptr && recursive)
		{
			vector<Interval>* right = this->right->queryPoint(x);
			copy(right->begin(), right->end(), back_insert_iterator<vector<Interval> >(*output));
			delete right;
			right = nullptr;
		}
		auto curr = this->ends.rbegin();
		auto end = this->ends.rend();
		while (curr != end && curr->end >= x)
		{
			output->push_back(*curr);
			++curr;
		}
	}
	return output;
}

void IntervalTree::walkPoints(resultSet *results, Interval &other)
{
    if (this->left != nullptr && other.start < this->center) this->left->walkPoints(results, other);
    if (this->right != nullptr && other.end > this->center) this->right->walkPoints(results, other);
    auto lower = this->points.begin();
    auto upper = this->points.end();
    //binary search to find the first point >= the start of the interval
    while (lower != upper && lower->first < other.start)
    {
        auto target = lower + (upper-lower);
        if (target == upper) break;
        lower = target;
    }
    //now add intervals into the results set as long as they have at least one point (start or end) <= the end of the interval
    while (lower != upper && lower->first <= other.end)
    {
        results->insert(*(lower->second));
        ++lower;
    }
}

resultSet* IntervalTree::queryInterval(Interval &other)
{
	resultSet *results = new resultSet(compIntervalStart);
    //resultSet *results = new resultSet();
    this->walkPoints(results, other);
	//now check for intervals which entirely encompass this interval by doing a point query on this interval's start
	//really, any point inside this interval would do just as fine
	vector<Interval>* encompassing = this->queryPoint(other.start, false);
	//copy these results into the set (to avoid duplicate intervals)
	copy(encompassing->begin(), encompassing->end(), inserter(*results, results->end()));
    
//	delete encompassing;
//	encompassing = nullptr;
//	vector<Interval>* output = new vector<Interval>();
//	//copy back out to a vector
//	copy(results.begin(), results.end(), back_insert_iterator<vector<Interval> >(*output));
//	return output;
    return results;
}

IntervalTree* constructIntervalTree(vector<Interval> &intervals)
{
	//sort intervals first, so we can do our best to keep the tree roughly balanced during construction
	sort(intervals.begin(), intervals.end(), compIntervalStart);
	return new IntervalTree(intervals);
}
