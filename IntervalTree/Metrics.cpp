//
//  Metrics.cpp
//  IntervalTree
//
//  Created by Aaron Graubert on 7/5/17.
//  Copyright © 2017 Aaron Graubert. All rights reserved.
//

#include "Metrics.h"


/*void increment(std::string);
 
 friend std::ofstream& operator<<(std::ofstream&, Metrics&);
 };
 
 
*/

void Metrics::increment(std::string key)
{
    this->counter[key]++;
}

void Metrics::increment(std::string key, int n)
{
    this->counter[key] += n;
}

unsigned long Metrics::get(std::string key)
{
    return this->counter[key];
}

double Metrics::frac(std::string a, std::string b)
{
    return (double) this->get(a) / this->get(b);
}

std::ofstream& operator<<(std::ofstream &stream, const Metrics &counter)
{
    auto beg = counter.counter.begin();
    auto end = counter.counter.end();
    while (beg != end)
    {
        stream << beg->first << "\t" << beg->second << std::endl;
        ++beg;
    }
    return stream;
}
