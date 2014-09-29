/*
Model integration example 2: csv.cpp
	
	  Copyright 2014 Matthew V Talluto, Isabelle Boulangeat, Dominique Gravel
	
	  This program is free software; you can redistribute it and/or modify
	  it under the terms of the GNU General Public License as published by
	  the Free Software Foundation; either version 3 of the License, or (at
	  your option) any later version.
	  
	  This program is distributed in the hope that it will be useful, but
	  WITHOUT ANY WARRANTY; without even the implied warranty of
	  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
	  General Public License for more details.
	
	  You should have received a copy of the GNU General Public License
	  along with this program; if not, write to the Free Software
	  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
	
	


	Rudimentary CSV parser to assist reading input data for the sampler
	See the documentation in csv.hpp for information on using this class

*/


#include "csv.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>


using std::string;
using std::vector;

CSV::CSV (const char * filename, size_t header)
{
	// parse the input file
	std::ifstream inputFile;
	inputFile.open(filename);
	if(!inputFile.is_open()) {
		std::stringstream err;
		err << "Failed to open file <" << filename << ">\n";
		throw std::runtime_error(err.str());
	}
	string line;
	
	if(header > 0) {
		for(size_t i = 0; i < header && getline(inputFile, line); i++)
			continue;
	}

	while(std::getline(inputFile, line)) {
		vector<string> row = split_line(line);
		_data.push_back(str_to_double(row));
	}
	inputFile.close();
}


CSV::CSV()
{ }


vector<vector<double> > CSV::data() const
{
	return _data;
}


vector<double> CSV::columns(size_t colNum) const
{
	vector<double> result;
	for(vector<vector<double> >::const_iterator row = _data.begin(); row != _data.end(); row++)
		result.push_back((*row).at(colNum));
	return result;
}


vector<vector<double> > CSV::columns(size_t start, size_t end) const
{
	vector<vector<double> > result;
	for(vector<vector<double> >::const_iterator row = _data.begin(); row != _data.end(); row++) {
		vector<double> cols;
		for(size_t i = start; i < end; i++) {
			cols.push_back((*row).at(i));
		}
		result.push_back(cols);
	}
	
	return result;
}


vector<string> CSV::split_line(const string &str)
{
    std::stringstream lineStream(str);
    string cell;
    vector<string> result;

    while(std::getline(lineStream,cell,','))
    {
        result.push_back(cell);
    }
    return result;
}


vector<double> CSV::str_to_double(const vector<string> &dat)
{
	vector<double> result;
	for(vector<string>::const_iterator cellDat = dat.begin(); cellDat != dat.end(); cellDat++) {
		std::stringstream strDat(*cellDat);
		double val;
		strDat >> val;
		result.push_back(val);
	}
	return result;
}
