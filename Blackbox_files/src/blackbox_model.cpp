// blackbox_model.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <math.h> /* pow */
#include <string.h>
#include <sstream>
#include "program_options.hpp"

// #include <boost/program_options.hpp>
using std::string;
using namespace std;

int main(int argc, char** argv)
{

	// Parse command line arguments
	try {
		program_options::parse(argc, argv);
	} catch (const std::exception &x) {
		std::cerr << "bb: " << x.what() << '\n';
		std::cerr << "usage: bb [-p|--parameters] <parameter float> [-x|--variables] <variable float> [-v|--verbose] [-g|--gradients] ...\n";
		return EXIT_FAILURE;
	}


	std::vector<double> x = program_options::variables();
	std::vector<double> p = program_options::parameters();
	std::vector<double> g, d_obj, d_g1, d_g2, d_g3;
	double obj;

	// display inputs
	if (program_options::verbose()) {
		cout << "x = [ ";
		for (size_t i = 0; i < x.size(); ++i)
		{
			cout << x[i] << ", ";
		}
		cout << "]\n";
		cout << "p = [ ";
			for (size_t i = 0; i < p.size(); ++i)
		{
			cout << p[i] << ", ";
		}
		cout << "]\n";
	}

	obj = pow(pow(x[0], 2) + x[1] - 11, 2) + pow(x[0] + pow(x[1], 2) - 7,2);
	d_obj.push_back(4 * x[0] * (pow(x[0], 2) + x[1] - 11) + 2 * (x[0] + pow(x[1], 2) - 7));
	d_obj.push_back(2 * (pow(x[0], 2) + x[1] - 11) + 4 * x[1] * (x[0] + pow(x[1], 2) - 7));

	g.push_back(-(pow(x[0] - p[0], 2) - (x[1] - p[1]) - 5));
	d_g1.push_back(2*p[0]-2*x[0]);
	d_g1.push_back(1);

	g.push_back(-(pow(x[0] - p[2], 2) + pow(x[1] - p[3], 2) - 4));
	d_g2.push_back(2*p[2]-2*x[0]);
	d_g2.push_back(2*p[3]-2*x[1]);

	g.push_back(-(pow(p[6]*(x[0] - p[4]), 2) + pow(p[7]*(x[1] - p[5]), 2) - 8));
	d_g3.push_back(-p[6]*(-2*p[4] + 2*x[0]));
	d_g3.push_back(-p[7]*(-2*p[5] + 2*x[1]));

	if (program_options::verbose()) {
		cout << "f(x): " << obj << ", ";
		cout << "g(x): [ " << g[0] << " " << g[1] << " " << g[2] << " ]\n";
		if (program_options::return_gradients()) {
			cout << "df(x): [ " << d_obj[0] << " " << d_obj[1] << " ]\n";
			cout << "dg1(x): [ " << d_g1[0] << " " << d_g1[1] << " ]\n";
			cout << "dg2(x): [ " << d_g2[0] << " " << d_g2[1] << " ]\n";
			cout << "dg3(x): [ " << d_g3[0] << " " << d_g3[1] << " ]\n";
		}
	}

	ofstream outfile("output.txt");
	outfile.precision(7);

	if (outfile.is_open())
	{
		outfile << obj << "\n";
		outfile << g[0] << " " << g[1] << " " << g[2] << "\n";
		if (program_options::return_gradients()) {
			outfile << d_obj[0] << " " << d_obj[1] << "\n";
			outfile << d_g1[0] << " " << d_g1[1] << "\n";
			outfile << d_g2[0] << " " << d_g2[1] << "\n";
			outfile << d_g3[0] << " " << d_g3[1];
		}
		outfile.close();
	}
	else cout << "Unable to open file\n";
	return 0;
}