#include <iostream>
#include <vector>
using namespace std;
class apple
{
public:
	apple(int num);
	apple();
	~apple();
	
	int setWeight(double );
	double  getWeight();

private:
	vector<double> vec;
	double weight; 
};