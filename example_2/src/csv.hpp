#ifndef CSV_H
#define CSV_H

#include <vector>
#include <string>
#include <iostream>

class CSV {
  public:
  	std::vector<std::vector<double> > data() const;
  	std::vector<double> columns(size_t colNum) const ;
  	std::vector<std::vector<double> > columns(size_t start, size_t end) const;

	/*
		CSV object constructor
		filename: name of the csv file to open
		header: number of header rows to skip
		
		throws a runtime error if opening file fails
	*/
  	CSV(const char * filename, size_t header = 0);
  	CSV();

  private:
	std::vector<std::vector<double> > _data;
  	
	std::vector<double> str_to_double(const std::vector<std::string> &dat);
	std::vector<std::string> split_line(const std::string &str);

};

#endif